# JWST/CEERS Disk Galaxy Detection & Bar Dynamics Pipeline

A research pipeline for detecting disk-type galaxy morphology in JWST/NIRCam
imaging from the CEERS survey, as a precursor to identifying barred disks and
characterizing bar dynamics.

**Status: Stage 1, early development.** The repository currently contains
raw CEERS mosaics and lightweight utilities for inspecting them. The
detection/classification pipeline described below is planned but not yet
implemented (see [Roadmap](#roadmap--status)).

---

## Scientific Background

### Data Source

This project uses imaging from the **Cosmic Evolution Early Release Science
Survey (CEERS)**, a JWST Director's Discretionary Early Release Science
program targeting the Extended Groth Strip (EGS) extragalactic field. CEERS
data is released in epochs and data-release versions (this project uses
**DR0.6**), combining imaging from two of JWST's instruments:

- **NIRCam** (Near Infrared Camera) — near-infrared imaging, ~0.6–5.0 μm,
  split into a short-wavelength (SW) and long-wavelength (LW) channel
  captured simultaneously. Higher angular resolution (PSF FWHM ~0.03–0.15″
  depending on filter), large field of view, HgCdTe detectors. This is the
  primary imaging dataset for this project.
- **MIRI** (Mid-Infrared Instrument) — mid-infrared imaging, ~5–25 μm,
  Si:As detectors, coarser resolution (PSF FWHM ~0.2–0.9″), smaller field
  of view. Traces cooler dust and PAH emission features; not used for
  morphological analysis in this project due to insufficient resolution
  for structural decomposition.

CEERS observations were scheduled using JWST's **coordinated parallel
observing** mode: at any given pointing, one instrument is designated
*prime* (its science requirements set the telescope's exact pointing, roll
angle, and dither pattern), while a second instrument observes *in
parallel*, capturing whatever falls within its own offset field of view.
In Epoch 1 (June 2022), MIRI was prime with NIRCam in parallel; in Epoch 2
(December 2022), NIRSpec (and NIRCam WFSS) was prime with NIRCam/MIRI
imaging in parallel, observed at a 180° position-angle rotation relative to
Epoch 1 to help fill coverage gaps caused by the instruments' differing
fields of view. As a result, NIRCam and MIRI coverage across CEERS
pointings is not uniform — not every pointing has every filter.

### NIRCam Filters Used

| Filter | ~Wavelength | Channel | Notes |
|---|---|---|---|
| F115W | 1.15 μm | SW | Rest-frame UV at z~1–3; traces young stars/star formation |
| F150W | 1.50 μm | SW | Rest-frame UV/blue-optical |
| F200W | 2.00 μm | SW | Transitional band |
| F277W | 2.77 μm | LW | Rest-frame optical at z~1–2; common detection band |
| F356W | 3.56 μm | LW | Rest-frame optical/near-IR; paired with F277W for source detection |
| F410M | 4.10 μm | LW (medium) | Narrow bandpass, lower per-pixel S/N |
| **F444W** | 4.44 μm | LW | **Reddest NIRCam broadband; primary science filter for this project.** Rest-frame optical/near-IR at z~1–4, least affected by dust extinction, best traces the underlying (old) stellar mass distribution rather than star-forming clumps. Coarsest PSF among NIRCam filters (~0.145″ FWHM). |

**Filter choice rationale:** This project uses **F444W alone**. Bars and
disks are structural features of the evolved stellar population, best
traced in rest-frame optical/near-IR light rather than rest-frame UV
(which is dominated by young star clusters and dust, and can obscure the
smooth underlying stellar distribution). For the target redshift range
(z~1–4), F444W samples closer to this rest-frame optical/near-IR regime
than the bluer NIRCam filters, making it the physically appropriate
single-band choice for both disk detection and, in later stages, bar
identification. A multi-filter detection stack (e.g. F277W+F356W, as used
in several CEERS reference catalogs) would improve completeness for faint
sources, but is not required for morphological measurement once a source
is detected, and is not used at this stage of the project.

### Redshift (z)

Redshift quantifies the stretching of a galaxy's emitted light toward
longer wavelengths due to the expansion of the universe:

```
z = (λ_observed − λ_emitted) / λ_emitted
```

- **Spectroscopic redshift** is measured by identifying known atomic
  transition wavelengths (e.g. Hα, Lyman-α) in a resolved spectrum and
  comparing their observed position to laboratory rest-frame values.
- **Photometric redshift ("photo-z")**, used for most sources in broadband
  imaging surveys like CEERS, is estimated by fitting template galaxy
  spectra (redshifted across a grid of trial z values) to the observed
  fluxes in each broadband filter, and selecting the best-fit
  template/z combination (e.g. via EAZY, BPZ, LePhare). This is the likely
  origin of redshift estimates for sources in the working catalog (~1000
  galaxies) unless spectroscopic (NIRSpec) redshifts are available for a
  given source.

Redshift matters for filter interpretation because a fixed *observed*
wavelength band (e.g. F444W) samples different *rest-frame* light
depending on a galaxy's z — directly affecting which physical structures
(old stellar disk vs. young star-forming clumps) dominate the observed
morphology.

### Pipeline Goal (Staged)

1. **Stage 1 (current):** Detect galaxies and classify disk-type
   morphology in F444W imaging.
2. **Stage 2 (planned):** Identify barred disk galaxies among detected
   disks.
3. **Stage 3 (planned):** Characterize bar dynamics (e.g. pattern speed,
   corotation radius) for identified barred galaxies.

### Stage 1 Method

- **Source detection & deblending:** `photutils.segmentation`
- **Morphological measurement:** `statmorph` — Gini, M20, concentration,
  asymmetry, Sérsic index, per detected source
- **Disk classification (non-parametric cut):**
  - Sérsic index: `n < 2.5`
  - Gini–M20 locus: `G < -0.14·M20 + 0.33`

This approach was chosen over a CNN/GCNN classification route because it
does not require labeled training data or pretrained weights, at the cost
of being a coarser, threshold-based classification rather than a learned
one.

### Known Caveats

- PSF resolution in F444W (~0.145″ FWHM) is the coarsest among NIRCam
  filters used in CEERS; this is not expected to significantly affect
  Stage 1 disk/non-disk classification, but may affect the precision of
  Stage 2/3 bar identification and dynamics for compact or distant
  sources, since bar detection (ellipse fitting, Fourier m=2
  decomposition) is more sensitive to PSF smoothing than global
  Sérsic/Gini-M20 classification.
- Detecting sources in F444W alone (rather than a multi-band stack) may
  reduce completeness for faint sources relative to CEERS reference
  catalogs, but does not affect the validity of morphology measurements
  for sources that are detected.
- Not all CEERS pointings/filters have equal or complete coverage due to
  the prime/parallel observing scheme described above; this should be
  accounted for when combining catalogs across multiple pointings.

---

## Repository Layout

```
.
├── src/
│   ├── open_fits_images.py   # interactive FITS viewer (ZScale, NaN-aware)
│   └── high-res.py           # export a FITS science array to a full-resolution PNG
├── images/                   # CEERS DR0.6 F444W mosaics (gitignored — not tracked in git)
│   └── compressed-images/    # gzip-compressed copies of the mosaics above
└── requirements.txt
```

`images/` is listed in `.gitignore` because the mosaics are large
(2–3.5 GB each); this repository ships code, not data. See
[Getting the data](#getting-the-data) below for where the pointings used
here come from.

## Getting the Data

This project currently works with CEERS DR0.6 NIRCam F444W mosaics for
**pointings 4, 5, 7, 8, 9, and 10**. Mosaics are `hlsp_ceers_jwst_nircam_*`
FITS files, one per pointing, obtained from the CEERS public data release
(MAST / CEERS HLSP). Place downloaded mosaics in `images/` following the
existing naming convention:

```
images/hlsp_ceers_jwst_nircam_nircam<N>_f444w_dr0.6_i2d.fits
```

Each FITS file is a multi-extension file; the scripts in this repo read
the `SCI` (science image) extension. `ERR` (error) and `WHT`
(weight/inverse-variance) extensions are also present in these files and
will be used by later pipeline stages.

There is no automated download script yet — see [Roadmap](#roadmap--status).

## Installation

Requires Python 3.12+.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

The two scripts in `src/` are standalone command-line utilities for
inspecting a single FITS mosaic. Both scripts:

- Read the `SCI` extension only.
- Treat non-finite pixels (NaNs in blank padding regions around the
  detector footprint) as masked, and compute display scaling
  (`ZScaleInterval`) from **finite pixels only** — this avoids the
  all-black / washed-out images that result from including NaN padding in
  the stretch calculation.

### View a mosaic interactively

```bash
python3 src/open_fits_images.py images/hlsp_ceers_jwst_nircam_nircam4_f444w_dr0.6_i2d.fits
```

Opens a matplotlib window showing the image with ZScale contrast, prints
basic diagnostics (finite-pixel fraction, min/max/median, ZScale limits) to
the terminal, and displays a colorbar in surface-brightness units
(MJy/sr).

### Export a full-resolution PNG

```bash
python3 src/high-res.py images/hlsp_ceers_jwst_nircam_nircam4_f444w_dr0.6_i2d.fits images/nircam4_highres.png
```

Writes a pixel-for-pixel PNG of the science array (no resampling, no
axes/labels) — the sharpest possible rendering of the mosaic, at the cost
of file size (mosaics can be tens to hundreds of megapixels, so output
files may run to hundreds of MB). Useful for zooming into individual
sources by eye ahead of automated detection. An alternative, labeled
high-DPI export path (with axes, title, and colorbar, capped by DPI ×
figure size rather than true native resolution) is sketched in a commented
block at the bottom of the script for presentation-quality figures.

## Roadmap / Status

Stage 1 (galaxy detection + disk classification in F444W) is in progress.
Stages 2–3 (bar identification, bar dynamics) are out of scope until Stage
1 produces a validated disk catalog.

- [x] FITS viewer with NaN-aware display scaling (`open_fits_images.py`)
- [x] Full-resolution PNG export utility (`high-res.py`)
- [ ] Script to download/organize CEERS DR0.6 F444W mosaics for pointings
      4, 5, 7, 8, 9, 10
- [ ] FITS I/O utilities: load `SCI`/`ERR`/`WHT` extensions consistently
      across scripts, correctly handling NaN/blank padding
- [ ] (Optional) Reproject/coadd multiple pointings into a single mosaic
      per filter (`reproject.mosaicking.reproject_and_coadd`), for cases
      where a continuous image is needed (visualization, sources
      straddling tile boundaries) — per-pointing catalogs are otherwise
      preferred over a full combined mosaic
- [ ] Background estimation and source detection/deblending
      (`photutils.segmentation`) on F444W mosaics
- [ ] Per-source morphology measurement (`statmorph`): Gini, M20,
      concentration, asymmetry, Sérsic index
- [ ] Disk classification cut: `n < 2.5` AND `G < -0.14*M20 + 0.33`
- [ ] Structured output catalog (CSV/FITS table): source ID, RA/Dec (from
      WCS), pointing ID, morphological parameters, disk classification
      flag
- [ ] Visual inspection utility: cutout thumbnails of classified disk
      candidates for manual sanity-checking
- [ ] Sanity statistics: sources detected per pointing, fraction
      classified as disks, distribution plots of Sérsic index and
      Gini-M20 with the classification boundary overlaid

Explicitly **out of scope for now**: Stage 2/3 bar identification and
dynamics, and CNN/GCNN-based classification (noted as a possible future
alternative to the current threshold-based approach, not currently
planned).

## Requirements

See `requirements.txt`. Core dependencies: `astropy`, `numpy`,
`matplotlib`. Planned additions for later Stage 1 work: `photutils`,
`statmorph`, `reproject`.
