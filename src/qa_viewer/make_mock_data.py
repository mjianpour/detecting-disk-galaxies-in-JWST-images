"""Generate mock pipeline outputs so the QA app runs without real data.

Creates, under mock_data/:
  - catalog.csv                mock SourceCatalog + statmorph metrics
  - mosaic.fits                fake full pointing (SCI HDU) with injected sources
  - cutouts/{source_id}.fits   per-source stamps: SCI HDU + SEG HDU (mask)

Sources are Sersic2D profiles: low-n ("disky") and high-n ("spheroid-like")
mixed, with statmorph-style metrics drawn to correlate with the input n so
the disk classification cuts (n < 2.5 and G < -0.14*M20 + 0.33) do
something sensible in the app.

Usage: python make_mock_data.py [output_dir]
"""

import os
import sys

import numpy as np
import pandas as pd
from astropy.io import fits

RNG = np.random.default_rng(42)

MOSAIC_SHAPE = (900, 900)
N_SOURCES = 18
CUTOUT_HALF = 32          # cutout is (2*CUTOUT_HALF)^2 pixels
NOISE_SIGMA = 0.02        # background noise, arbitrary surface-brightness units
PIXSCALE = 0.03           # arcsec/pix, mock WCS-ish numbers for ra/dec
RA0, DEC0 = 214.90, 52.90  # roughly EGS field


def make_sources():
    """Random source parameters, kept away from mosaic edges."""
    margin = CUTOUT_HALF + 10
    sources = []
    for i in range(N_SOURCES):
        disky = RNG.random() < 0.55
        n = RNG.uniform(0.7, 1.8) if disky else RNG.uniform(2.8, 4.5)
        sources.append(dict(
            source_id=i + 1,
            x=RNG.uniform(margin, MOSAIC_SHAPE[1] - margin),
            y=RNG.uniform(margin, MOSAIC_SHAPE[0] - margin),
            amplitude=RNG.uniform(0.3, 1.5),
            r_eff=RNG.uniform(3.0, 9.0),
            n=n,
            ellip=RNG.uniform(0.05, 0.6) if disky else RNG.uniform(0.0, 0.3),
            theta=RNG.uniform(0, np.pi),
        ))
    return sources


def sersic_profile(xx, yy, s):
    """Elliptical Sersic surface brightness (same parametrization as
    astropy's Sersic2D, but with the Ciotti & Bertin (1999) approximation
    for b_n to avoid a scipy dependency)."""
    n = s["n"]
    bn = 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n**2)
    cos_t, sin_t = np.cos(s["theta"]), np.sin(s["theta"])
    x_maj = (xx - s["x"]) * cos_t + (yy - s["y"]) * sin_t
    x_min = -(xx - s["x"]) * sin_t + (yy - s["y"]) * cos_t
    r = np.sqrt(x_maj**2 + (x_min / (1 - s["ellip"]))**2)
    return s["amplitude"] * np.exp(-bn * ((r / s["r_eff"])**(1 / n) - 1))


def render_mosaic(sources):
    yy, xx = np.mgrid[0:MOSAIC_SHAPE[0], 0:MOSAIC_SHAPE[1]]
    image = RNG.normal(0.0, NOISE_SIGMA, MOSAIC_SHAPE)
    for s in sources:
        # only evaluate in a local box for speed
        h = CUTOUT_HALF + 10
        x0, x1 = int(s["x"]) - h, int(s["x"]) + h
        y0, y1 = int(s["y"]) - h, int(s["y"]) + h
        image[y0:y1, x0:x1] += sersic_profile(xx[y0:y1, x0:x1], yy[y0:y1, x0:x1], s)
    return image


def mock_metrics(s):
    """statmorph-style metrics loosely correlated with the input profile."""
    n_meas = s["n"] * RNG.normal(1.0, 0.12)
    if s["n"] < 2.5:  # disky inputs: low Gini, more negative-ish M20 relation
        gini = RNG.normal(0.42, 0.04)
        m20 = RNG.normal(-1.55, 0.20)
        conc = RNG.normal(2.7, 0.25)
    else:
        gini = RNG.normal(0.58, 0.04)
        m20 = RNG.normal(-2.05, 0.20)
        conc = RNG.normal(4.0, 0.3)
    asym = np.clip(RNG.normal(0.08, 0.06), 0, None)
    # a few sources get bad-fit flags to exercise the QA display
    flag = int(RNG.random() < 0.12)
    flag_sersic = int(RNG.random() < 0.15)
    is_disk = (n_meas < 2.5) and (gini < -0.14 * m20 + 0.33)
    return dict(sersic_n=round(n_meas, 3), gini=round(gini, 4), m20=round(m20, 4),
                concentration=round(conc, 3), asymmetry=round(asym, 4),
                is_disk=is_disk, flag=flag, flag_sersic=flag_sersic)


def make_seg_mask(cut, source_id):
    """Simple threshold-based segmentation label map for the cutout."""
    seg = np.zeros(cut.shape, dtype=np.int32)
    seg[cut > 4 * NOISE_SIGMA] = source_id
    return seg


def main(out_dir):
    cutout_dir = os.path.join(out_dir, "cutouts")
    os.makedirs(cutout_dir, exist_ok=True)

    sources = make_sources()
    mosaic = render_mosaic(sources)

    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(mosaic.astype(np.float32), name="SCI"),
    ]).writeto(os.path.join(out_dir, "mosaic.fits"), overwrite=True)

    rows = []
    for s in sources:
        sid = s["source_id"]
        xi, yi = int(round(s["x"])), int(round(s["y"]))
        cut = mosaic[yi - CUTOUT_HALF:yi + CUTOUT_HALF,
                     xi - CUTOUT_HALF:xi + CUTOUT_HALF]
        seg = make_seg_mask(cut, sid)

        sci_hdu = fits.ImageHDU(cut.astype(np.float32), name="SCI")
        sci_hdu.header["SRCID"] = sid
        fits.HDUList([
            fits.PrimaryHDU(),
            sci_hdu,
            fits.ImageHDU(seg, name="SEG"),
        ]).writeto(os.path.join(cutout_dir, f"{sid}.fits"), overwrite=True)

        rows.append(dict(
            source_id=sid, x=round(s["x"], 2), y=round(s["y"], 2),
            ra=round(RA0 + (s["x"] - MOSAIC_SHAPE[1] / 2) * PIXSCALE / 3600, 6),
            dec=round(DEC0 + (s["y"] - MOSAIC_SHAPE[0] / 2) * PIXSCALE / 3600, 6),
            **mock_metrics(s),
        ))

    catalog = pd.DataFrame(rows)
    catalog.to_csv(os.path.join(out_dir, "catalog.csv"), index=False)
    print(f"Wrote {len(catalog)} sources to {out_dir}/")
    print(catalog[["source_id", "sersic_n", "gini", "m20", "is_disk"]].to_string(index=False))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "mock_data")
    main(out)
