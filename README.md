# Introduction

This pipeline has implemented to take a large CEERS mosaic and identify barred galaxies from all detected sources in that large image. 

Expected inputs could be downlaoded [from here](https://ceers.github.io/dr06.html) and expected filter is F444W, since the stellar bar of the galaxy has formed better in this filter (if there is any!). 

Pipleline will create some files in the `/src` path of porject after the job is done. 

```text
/src
|____ceers_barred_galaxies_report.pdf
|____/barred_galaxies/
     |____/strongly_barred/...
     |____/weakly_barred/...
|____/cutouts_physical/... 
```

After the job is done ... 

- A `.pdf` file named `ceers_barred_galaxies_report.pdf` is created to present a report of detected barred galaxies and their properties such as **bar eccentricity, PA, bar length, RA, declanation, etc.** 
- Every detected galaxies are stored in `/cutouts_physical` folder as individual `.fits` files. 
- Barred galaxies are stored in the folder `/barred_galaxies` in two separate folders `/strongly_barred` and `/weakly_barred` as individual `.fits` files, each file standing for barred galaxy. 

## JWST Properties & Characteristics of the Filter F444W 

- Each pixel corresponds to $0.0300\ \text {arcsec}$ on the sky
    > Source: https://ceers.github.io/dr06.html#:~:text=The%20mosaics%20are%20on%20pixel%20scales%20of%200.03%22/pix%20and%20are%20pixel%2Daligned%20across%20all%20filters.%20The%20science%20images%20are%20in%20units%20of%20MegaJansky%20per%20steradian%20(MJy/sr)%2C%20and%20so%20the%20zeropoint%20is%2028.08652%20AB%20mag.
- $FWHM$: $4.70\ \text {pix}$  
- $PSF$ Area: $17.3$ ${\text {pix}}^2$ 
# Contents

## How to Use this Repository.
Note first that the original end to end pipeline has implemented in `/src/full_test.ipynb`. 
### 1. Dependencies & Installation 
Visit the file `requirements.txt` in the root of project to review the requirements of project. You can isntall all the dependencies of this project by running this command in terminal while being in the root of the porject. 

```txt
pip install -r requirements.py 
```
### 2. Input CEERS mosaic
Enter the path to your `.fits` file in the third cell of notebook, this way: 

```python
image = '../path/to/your/file.fits' # e.g. image = '../images/hlsp_ceers_jwst_nircam_nircam8_f444w_dr0.6_i2d.fits'
file = fits.open(image)
```
### 3. Run the Notebook
Run the notebook and note this might take a while. When the job is done, last notebook will log this message: 

```text
Saved: ../detecting-disk-galaxies-in-JWST images/src/ceers_barred_galaxies_report.pdf (XX barred candidates / XXX detected galaxies)
```
### 4. Better Displaying
You can display an individual exported `.fits` galaxy using the notebook `better displaying.ipynb`. Just enter the path to your file, run the cell and you're good to go. 

## Noise Reduction Methods

### Pre-Reduced Noises
First thing to note is that `hlsp_ceers_jwst_nircam_nircam8_f444w_dr0.6_i2d.fits` and similar files, downloaded [from here](https://ceers.github.io/dr06.html) have gone through a noise reduction pipeline before being published. Noises such as, snowballs, wisps, $\frac{1}{f}$ noise have subtracted earlier by CEERS staff. 

> Source: https://ceers.github.io/dr06.html#:~:text=WFC3%20F160W-,DATA%20REDUCTION,Pipeline%20parameter%20files%20are%20available%20on%20GitHub%20at%20ceers%2Dnircam.,-EPOCH%201%262
## Source Detection & Filtering

### Sky Threshold
We could use the header `BKG_SUB` of the original `hlsp_ceers_jwst_nircam_nircam8_f444w_dr0.6_i2d.fits` file from the beginning to gain access to the sources (objects in the sky), with no sky pixel included. But, for education purposes we decided to do subtract the sky pixels by ourselves. We found and subtracted the sky pixels in this way: 

1. Defining a global sky threshold
2. Defining a global noise spread value
3. Any pixel above the sky and noise spread will be known as a source. 
4. Replacing the global sky threshold and noise spread by corresponding values at each $\sim 200 * 200\ \text{pix}^2$ tile. 
5. Subtractnig again from the raw file (no sky and noise pixel subtracted!)
6. Every remaining pixels will be known as a source at this stage. 

Following cells do these jobs: 

```python
h, w = sci.shape
tile = 200
n_rows, n_cols = h // tile, w // tile
medians = np.full((n_rows, n_cols), np.nan)
spreads = np.full((n_rows, n_cols), np.nan)

for i in range(n_rows):
    for j in range(n_cols):
        box      = sci[i*tile:(i+1)*tile, j*tile:(j+1)*tile]
        box_mask = nodata[i*tile:(i+1)*tile, j*tile:(j+1)*tile]
        real = box[~box_mask]
        if real.size < 0.5 * box.size:
            continue
        mean, med, stddev = sigma_clipped_stats(real)
        medians[i, j] = med
        spreads[i, j] = stddev
```

```python
sky   = np.nanmedian(medians)
noise = np.nanmedian(spreads)
rms   = np.sqrt(np.nanmedian(spreads**2))
n = 5
threshold = sky + n * noise
suspects = (sci > threshold) & ~nodata
print(suspects.mean())
```

```python
# Use the entire mosaic rather than a 2000 x 2000 testing tile.
y0, y1 = 0, sci.shape[0]
x0, x1 = 0, sci.shape[1]

cut = sci[y0:y1, x0:x1].astype(np.float64)
cut_mask = nodata[y0:y1, x0:x1]

print(f"Processing full mosaic: {cut.shape[1]} x {cut.shape[0]} pixels")
print(f"Invalid-pixel fraction: {cut_mask.mean():.4f}")
```

```python
from scipy.ndimage import zoom

# Replace empty sky tiles with the global estimates before interpolation.
medians_filled = np.where(np.isfinite(medians), medians, sky)
spreads_filled = np.where(np.isfinite(spreads), spreads, rms)

# Interpolate the tile-level background and noise maps to the full mosaic size.
zoom_factors = (
    cut.shape[0] / medians_filled.shape[0],
    cut.shape[1] / medians_filled.shape[1],
)

sky_cut = zoom(medians_filled, zoom_factors, order=1, mode="nearest")
rms_cut = zoom(spreads_filled, zoom_factors, order=1, mode="nearest")

# Ensure the interpolated arrays exactly match the mosaic dimensions.
sky_cut = sky_cut[:cut.shape[0], :cut.shape[1]]
rms_cut = rms_cut[:cut.shape[0], :cut.shape[1]]

if sky_cut.shape != cut.shape or rms_cut.shape != cut.shape:
    raise ValueError(
        f"Background-map shape mismatch: image={cut.shape}, "
        f"sky={sky_cut.shape}, noise={rms_cut.shape}"
    )

cut_sub = cut - sky_cut
cut_sub[~np.isfinite(cut_sub)] = 0.0

print(f"Full-mosaic background maps created: {sky_cut.shape}")
```
### Convolution & Deblending
We convolve the remaining pixels (sources) to make sure no pixel of a source is outsanding from the rest of the source pixels. So, the entire pixels of the source together known as a single source. It also helps not confusing a noise pixel with a source.

Following cell convolves the large image of sky: 

```python
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel

kernel = make_2dgaussian_kernel(5.0, size=11)
conv = convolve(cut_sub, kernel, mask=cut_mask)
segm_conv = detect_sources(conv, threshold=n * rms_cut, n_pixels=5, mask=cut_mask)
print(segm_conv.n_labels)
```
### Cataloging the Sky
Then we catalog the detected sources by running the follwing cell:

```python
segm_local = detect_sources(cut_sub, threshold=n * rms_cut, n_pixels=5, mask=cut_mask)
print(segm_local.n_labels)
```
### Filtering the Detected Souces
Detected sources are passed through a set of procedural filters and each filter will drop some of the sources. Those filters are: 

1. Negative flux sources (`flag_badflux`)
2. `Semi-minor axis of the ellipse < MIN_SEMIMINOR` or `ellipse eccentricity > maximum eccentricity` (`np.mean(ecc) + std(ecc)`) (`flag_unresolved`)
3. $\text {area of the source} < ~ 17.3\ \text {pix}^2$ (`psf_area_pix`) (`flag_small`)
4. Small signal to noise ratio `(SNR) < MIN_SNR`(`flag_lowsnr`)
5. Source too close to image edge (box within `EDGE_BUFFER` pixels of border) (`flag_edge`)
6. Too many bad pixels in source footprint (`badfrac > MAX_BADFRAC = 0.05`) (`flag_baddata`)
7. Mismatch between area and SNR normalized values (`|area_norm - snr_norm| > 3.0`) (`flag_snr_area_mismatch`)

remaining sources are all known as a galaxy. 
## Indentifying Barred Galaxies (First level Esitmation)
## Methodology
There are different ways to identify the stellar bar in an individual galaxy. The nucleus/bar of a galaxy is the brightest area of that galaxy. We take advantage of this fact to ... 
1. Identify the central brightest area of a given galaxy
2. Fit an ellipse to that area
3. We define $q = \frac{\text{minor-axis of ellipse}}{\text{major-axis of ellipse}}$  
4. If $q < 0.40$, galaxy will be known as **strongly barred**. 
5. If $0.40 < q < 0.60$, galaxy will be known as **weakly barred**. 
## Code Cells
## Bar Frcation
## Results
## Processing Extracted Barred Galaxies Precisely
This part is still in progress ... 
