from astropy.io import fits
from astropy.visualization import ZScaleInterval
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

# python3 export_high_res.py <path_to_fits_file> [output_path]

if len(sys.argv) < 2:
    print("Usage: python3 export_high_res.py <path to fits file> [output path]")
    sys.exit(1)

file_path = sys.argv[1]
if len(sys.argv) == 3:
    output_path = sys.argv[2]

with fits.open(file_path) as hdul:
    hdul.info()
    data = hdul['SCI'].data.astype(np.float64)
    header = hdul['SCI'].header

print(f"Image shape: {data.shape}  ({data.shape[1]} x {data.shape[0]} pixels)")

# --- Robust scaling, ignoring NaNs/blank padding ---
finite = np.isfinite(data)
finite_data = data[finite]
vmin, vmax = ZScaleInterval().get_limits(finite_data)
print(f"ZScale vmin/vmax: {vmin:.4g} / {vmax:.4g}")

# --- Normalize to 0-1, clip, replace NaNs with 0 (black) ---
norm = np.clip((data - vmin) / (vmax - vmin), 0, 1)
norm = np.nan_to_num(norm, nan=0.0)

# ============================================================
# OPTION 1: Pure pixel-for-pixel export (RECOMMENDED for max detail)
# One array pixel = one image pixel, no axes/labels, no resampling.
# This is the sharpest possible output and will be large (100s of MB - ~1GB
# depending on the array size, since it's basically lossless).
# ============================================================
plt.imsave(output_path, norm, cmap='gray', origin='lower', vmin=0, vmax=1)
print(f"Saved native-resolution image to: {output_path}")

# ============================================================
# OPTION 2: High-DPI figure WITH axes, ticks, colorbar (for presentations)
# Uncomment below if you want labeled axes instead of a bare pixel image.
# Note: this is capped by dpi x figsize, so it won't match true native
# resolution unless you push dpi very high (e.g. 600-1200) and a large figsize.
# ============================================================
# fig_output_path = output_path.replace('.png', '_labeled.png')
# fig_w = data.shape[1] / 300   # figure size in inches at target dpi
# fig_h = data.shape[0] / 300
# fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
# im = ax.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
# ax.set_title(file_path.split('/')[-1])
# cbar = fig.colorbar(im, ax=ax, label='MJy/sr')
# fig.savefig(fig_output_path, dpi=300, bbox_inches='tight')
# print(f"Saved labeled high-DPI image to: {fig_output_path}")