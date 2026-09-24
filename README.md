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

## JWST Properties & Characteristics
# Contents

## How to Use this Repository.
Note first that the original end to end pipeline has implemented in `/src/full_test.ipynb`. 
### 1. Installation 
You can isntall all the dependencies of this project by running this command in terminal while being in the root of the porject. 

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
### Convolution & Deblending
### Source Detection
### Filtering the Detected Souces
## Indentifying Barred Galaxies (First level Esitmation)
## Processing Extracted Barred Galaxies Precisely

