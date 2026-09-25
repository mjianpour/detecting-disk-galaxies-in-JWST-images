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

```python
from scipy.ndimage import (
    binary_closing,
    generate_binary_structure,
    label as nd_label,
    center_of_mass,
    sum as nd_sum,
)

def is_barred(result: dict):
    """Return 'weakly_barred', 'strongly_barred', or None for a galaxy cutout."""
    image_path = f"cutouts_physical/phys_src_{int(result['label']):05d}.fits"

    with fits.open(image_path) as hdul:
        sci_image = hdul["SCI"].data.astype(float)

    _, _, _ = sigma_clipped_stats(sci_image, sigma=3.0)

    # Smooth lightly to suppress pixel noise while retaining bar-like structure.
    kernel_light = make_2dgaussian_kernel(1.5, size=5)
    conv_image = convolve(sci_image, kernel_light)

    _, _, sky_std_conv = sigma_clipped_stats(conv_image, sigma=3.0)
    _, sky_median_conv, _ = sigma_clipped_stats(conv_image, sigma=3.0)
    conv_sub = conv_image - sky_median_conv

    peak_conv = np.nanmax(conv_sub)
    if not np.isfinite(peak_conv) or peak_conv <= 0:
        return None

    threshold = max(3.0 * sky_std_conv, 0.53 * peak_conv)
    binary = conv_sub > threshold

    structure = generate_binary_structure(2, 2)
    binary = binary_closing(binary, structure=structure, iterations=1)

    labels, n_labels = nd_label(binary, structure=structure)
    if n_labels == 0:
        return None

    if n_labels > 1:
        sizes = nd_sum(binary, labels, range(1, n_labels + 1))
        largest = int(np.argmax(sizes)) + 1
        cy, cx = center_of_mass(binary, labels, largest)
        r_keep = 1.5 * np.sqrt(sizes[largest - 1] / np.pi)

        keep = [largest]
        for label_id in range(1, n_labels + 1):
            if label_id == largest or sizes[label_id - 1] < 3:
                continue
            yi, xi = center_of_mass(binary, labels, label_id)
            if np.hypot(yi - cy, xi - cx) < r_keep:
                keep.append(label_id)

        binary = np.isin(labels, keep)

    n_mask_pix = int(binary.sum())
    if n_mask_pix == 0:
        return None

    r_eq = float(np.sqrt(n_mask_pix / np.pi))
    if r_eq < 4.0:
        return None

    ys, xs = np.where(binary)
    weights = np.clip(conv_sub[binary], 0.0, None)
    total_weight = weights.sum()
    if not np.isfinite(total_weight) or total_weight <= 0:
        return None

    x_center = float((weights * xs).sum() / total_weight)
    y_center = float((weights * ys).sum() / total_weight)

    dx = xs - x_center
    dy = ys - y_center
    mu_xx = float((weights * dx * dx).sum() / total_weight)
    mu_yy = float((weights * dy * dy).sum() / total_weight)
    mu_xy = float((weights * dx * dy).sum() / total_weight)

    covariance = np.array([[mu_xx, mu_xy], [mu_xy, mu_yy]])
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)

    if eigenvalues[1] <= 0:
        return None

    axis_ratio = float(np.sqrt(max(eigenvalues[0], 0.0) / eigenvalues[1]))
    if not np.isfinite(axis_ratio) or axis_ratio <= 0:
        return None

    a_fit = r_eq / np.sqrt(axis_ratio)

    # Broad, preliminary candidate filter; retain plausible cases for the profile test.
    if not (0.0 < axis_ratio < 0.75 and n_mask_pix >= 40 and a_fit >= 3.0):
        return None
        
    # Temporary bins based on mask elongation, not final bar strength.
    if axis_ratio < 0.50:
        return "strongly_barred"
    return "weakly_barred"
```

```python
from pathlib import Path
import shutil

barred_root = Path("barred_galaxies")
weak_dir = barred_root / "weakly_barred"
strong_dir = barred_root / "strongly_barred"

weak_dir.mkdir(parents=True, exist_ok=True)
strong_dir.mkdir(parents=True, exist_ok=True)

barred_counts = {
    "weakly_barred": 0,
    "strongly_barred": 0,
}

for result in physical_results:
    if not result["is_galaxy"]:
        continue

    category = is_barred(result)
    if category is None:
        continue

    label = int(result["label"])
    source_path = Path("cutouts_physical") / f"phys_src_{label:05d}.fits"
    destination_dir = weak_dir if category == "weakly_barred" else strong_dir

    if not source_path.exists():
        print(f"Missing cutout, skipped: {source_path}")
        continue

    shutil.copy2(source_path, destination_dir / source_path.name)
    barred_counts[category] += 1
    print(f"{category}: {source_path.name}")

print(
    f"Suspected barred galaxies copied: "
    f"{barred_counts['weakly_barred']} weakly barred, "
    f"{barred_counts['strongly_barred']} strongly barred"
)
```
## Results
These are a few galaxies detected as barred galaxies. 

<img width="1287" height="391" alt="image" src="https://github.com/user-attachments/assets/ec335f16-2e27-4154-9b8b-f16bfd5f2513" />

## Problems With Current Approach

1. Some individual galaxies are detected as siblings
2. Some bars are over / underestimated
3. Stellar bar area in multi-galactical sources are not well detected
4. Some galaxies are visually verified as an unbarred type, but still in the exports
5. We have no idea how many candidates are dropped in earlier filters and we can't announce an exact number for that. 

Note. Some of these problems will resolved by a more precise processing over candidates! 
## Processing Extracted Barred Galaxies Precisely
This part is still in progress ... 

---


