from astropy.io import fits
from astropy.visualization import ZScaleInterval
import numpy as np
import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Usage: python3 open_fits_images.py <path to the file>")
    sys.exit(1)

file_path = sys.argv[1]

with fits.open(file_path) as hdul:
    hdul.info()
    data = hdul['SCI'].data.astype(np.float64)
    header = hdul['SCI'].header

# --- Diagnostics: what's actually in this array? ---
finite = np.isfinite(data)
print(f"Total pixels: {data.size}")
print(f"Finite pixels: {finite.sum()} ({100*finite.sum()/data.size:.1f}%)")
print(f"Min/Max (finite only): {np.nanmin(data):.4g} / {np.nanmax(data):.4g}")
print(f"Median (finite only): {np.nanmedian(data):.4g}")

# --- Robust scaling: only use finite pixels ---
finite_data = data[finite]
vmin, vmax = ZScaleInterval().get_limits(finite_data)
print(f"ZScale vmin/vmax: {vmin:.4g} / {vmax:.4g}")

# Mask NaNs to a neutral value for display only (doesn't affect saved data)
display_data = np.where(finite, data, np.nan)

plt.figure(figsize=(8, 8))
plt.imshow(display_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
plt.colorbar(label='MJy/sr')
plt.title(file_path.split('/')[-1])
plt.show()