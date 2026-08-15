"""Pipeline output I/O — no Streamlit imports here.

This module is the boundary between the QA app and the detection/morphology
pipeline outputs. To point the app at real pipeline products, only the paths
(and possibly column names, see CATALOG_COLUMNS) need to change; the UI code
in app.py never touches FITS files or catalog parsing directly.

Conventions supported for segmentation-mask cutouts (see MaskConvention):
  - "same_file":     cutouts/{source_id}.fits with SCI in HDU 'SCI' (or 0/1)
                     and the mask in HDU 'SEG' (or the following HDU).
  - "parallel_dir":  cutouts/{source_id}.fits for science,
                     seg_cutouts/{source_id}.fits for the mask (HDU 0 or 'SEG').
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

# Columns the app knows how to display/filter. Only source_id is required;
# everything else degrades gracefully if absent.
CATALOG_COLUMNS = [
    "source_id", "x", "y", "ra", "dec",
    "sersic_n", "gini", "m20", "concentration", "asymmetry",
    "is_disk", "flag", "flag_sersic",
]

MASK_CONVENTIONS = ("same_file", "parallel_dir")


@dataclass(frozen=True)
class DataConfig:
    """Where the pipeline outputs live."""
    catalog_path: str
    cutouts_dir: str
    mask_convention: str = "same_file"   # one of MASK_CONVENTIONS
    seg_cutouts_dir: Optional[str] = None  # used when mask_convention == "parallel_dir"
    mosaic_path: Optional[str] = None      # full-pointing image (optional)


def load_catalog(catalog_path: str) -> pd.DataFrame:
    """Load the source catalog from CSV or a FITS table into a DataFrame."""
    ext = os.path.splitext(catalog_path)[1].lower()
    if ext in (".fits", ".fit", ".fits.gz"):
        df = Table.read(catalog_path).to_pandas()
    else:
        df = pd.read_csv(catalog_path)
    if "source_id" not in df.columns:
        raise ValueError(
            f"Catalog {catalog_path} has no 'source_id' column "
            f"(found: {list(df.columns)})"
        )
    if "is_disk" in df.columns:
        df["is_disk"] = df["is_disk"].astype(bool)
    return df


def _first_image_hdu(hdul: fits.HDUList, preferred_name: str) -> np.ndarray:
    """Return data from the named HDU if present, else the first HDU with data."""
    if preferred_name in hdul:
        return hdul[preferred_name].data
    for hdu in hdul:
        if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
            return hdu.data
    raise ValueError("No 2D image HDU found")


def load_cutout(config: DataConfig, source_id) -> tuple[np.ndarray, Optional[np.ndarray], fits.Header]:
    """Load (science, segmentation_mask, science_header) for one source.

    The mask is returned as-is (integer labels); None if not found.
    """
    sci_path = os.path.join(config.cutouts_dir, f"{source_id}.fits")
    with fits.open(sci_path) as hdul:
        sci = _first_image_hdu(hdul, "SCI").astype(np.float64)
        header = hdul["SCI"].header if "SCI" in hdul else hdul[0].header
        seg = None
        if config.mask_convention == "same_file":
            if "SEG" in hdul:
                seg = np.asarray(hdul["SEG"].data)
            else:
                # fall back: any integer-typed image HDU after the science one
                for hdu in hdul[1:]:
                    if hdu.data is not None and np.issubdtype(hdu.data.dtype, np.integer):
                        seg = np.asarray(hdu.data)
                        break

    if config.mask_convention == "parallel_dir" and config.seg_cutouts_dir:
        seg_path = os.path.join(config.seg_cutouts_dir, f"{source_id}.fits")
        if os.path.exists(seg_path):
            with fits.open(seg_path) as hdul:
                seg = np.asarray(_first_image_hdu(hdul, "SEG"))

    return sci, seg, header


def load_mosaic(mosaic_path: str, max_size: int = 1500) -> tuple[np.ndarray, int]:
    """Load the full-pointing science image, downsampled by simple striding
    so its largest dimension is <= max_size. Returns (image, stride)."""
    with fits.open(mosaic_path) as hdul:
        data = _first_image_hdu(hdul, "SCI").astype(np.float64)
    stride = max(1, int(np.ceil(max(data.shape) / max_size)))
    return data[::stride, ::stride], stride
