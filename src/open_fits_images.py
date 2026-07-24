from astropy.io import fits
from astropy.visualization import ZScaleInterval
import numpy as np
import sys
import matplotlib.pyplot as plt

# python3 open_fits_images.py your_f444w_image.fits <path to the file>

if len(sys.argv) != 2:
    print("Usage: python3 open_fits_images.py <path to the file>")
    sys.exit(1)

file_path = sys.argv[1]

with fits.open(file_path) as hdul:
    hdul.info()  # inspect extensions - JWST files often have SCI, ERR, DQ, etc.
    data = hdul['SCI'].data.astype(np.float64)   # science image
    header = hdul['SCI'].header
    wht = hdul['WHT'].data if 'WHT' in hdul else None  # weight/inverse-variance, if present
    err = hdul['ERR'].data if 'ERR' in hdul else None  # error array, if present

# Display the image using ZScale for better contrast
vmin, vmax = ZScaleInterval().get_limits(data)
plt.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
plt.show()