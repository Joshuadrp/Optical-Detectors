# Optical Detectors: Stellar Photometry Assignment

## Overview
This project implements a stellar photometry pipeline to process FITS images that were taken 
with two optical filters F336W and F555W, the first one providing images via UV light and the second
one via visible light, we know this due to its wavelengths 336nm and 555nm. The workflow is divided in four tasks,
given by the assignment briefing. Each of them is explained through the README, as well as their 
individual objectives, methods and results. The code is modular, utilizing functions efficiently and clearly. 

## Project Structure

```
Optical-Detectors/
├── data/
│   ├── F336W/          # UV filter FITS files (336 nm)
│   └── F555W/          # Visible filter FITS files (555 nm)
├── OpticalDetectors.py
└── README.md
```

## Usage

Run the complete pipeline:
```bash
python your-code.py
```
The script executes all four tasks sequentially and generates output files and visualizations.
It executes automatically, but once a plot is visualized, it needs to be closed for code to continue. 

## Description

### Task 1: Cosmic Ray Removal

Combines multiple FITS exposures using median stacking to remove cosmic ray hits and hot pixels, producing clean combined images for each filter.

**Method:**
- Reads all FITS files from each filter directory
- Computes median across exposures
- Outputs combined FITS files

**Outputs:**
- `combined_F336W.fits` - Cleaned UV filter image
- `combined_F555W.fits` - Cleaned visible filter image  
- `combined_images.png` - Side-by-side visualization

### Task 2: Star Detection, Cross-Matching & Star Catalog Created

Detects stellar sources, hopefully starts, in both filters using brightness 
thresholding with sigma-clipped background estimation, then cross-matches sources between filters by position, after this
a catalog containing stars within the two images is created as a text file.

**Method:**
- **Background estimation:** 3-iteration sigma-clipping to remove outlier pixels
- **Detection threshold:** σ × background_std (default σ = 1)
- **Size filtering:** 3-200 pixel area to reject noise and extended sources
- **Edge masking:** 50-pixel buffer from image borders
- **Cross-matching:** 2-pixel radius nearest-neighbor matching between filters
- **Position averaging:** Mean coordinates from both detections

**Outputs:**
- `matched_stars.txt` - Catalog of cross-matched source positions
- `star_detections.png` - Visualization showing detections in each filter and matched sources in a different image.

### Task 3: Aperture Photometry

Performs circular aperture photometry on all matched sources with local background 
subtraction to measure stellar fluxes and compute calibrated magnitudes.

**Method:**
- **Aperture:** 5-pixel radius circle centered on each source
- **Background annulus:** 8-12 pixel radii for local sky estimation
- **Background subtraction:** Median background × aperture area
- **Flux measurement:** Sum of background-subtracted aperture pixels
- **Zero-point calibration:** Brightest star set to magnitude 15.0
- **Magnitude formula:** mag = -2.5 × log₁₀(flux) + zero_point

**Outputs:**
- `photometry_catalog.txt` - Complete catalog with positions, aperture radius, magnitudes (F336W, F555W), and color

### Task 4: Hertzsprung-Russell Diagram

Creates a color-magnitude diagram plotting stellar color (F336W - F555W) versus magnitude 
to reveal the main sequence and other stellar populations.

**Method:**
- **X-axis:** Color = mag_F336W - mag_F555W (temperature proxy)
- **Y-axis:** Magnitude F336W (inverted, brightness indicator)
- **Plot format:** Scatter plot with inverted y-axis (bright at top)

**Outputs:**
- `hr_diagram.png` - HR diagram visualization

## Results

**Typical Performance:**
- ~1400-1600 matched sources detected
- Clear main sequence visible in HR diagram
- Color range: -6 to +8 mag (hot blue to cool red stars)
- Magnitude range: 14 to 30 (bright to faint)

## Key Features

- **Robust background estimation** using iterative sigma-clipping
- **Size filtering** to reject cosmic rays, hot pixels, and extended sources
- **Local background subtraction** via annulus for accurate photometry
- **Relative magnitude calibration** anchored to brightest star
- **Clean visualization** of all processing steps

## Technical Notes

**Filters:**
- **F336W (336 nm):** Near-UV, sensitive to hot blue stars
- **F555W (555 nm):** Visible green, approximates V-band

**Limitations:**
- Magnitude system is relative (not absolute calibration)

**Algorithm Choices:**
- Median combination removes cosmic rays more effectively than mean
- Sigma-clipping prevents bright stars from biasing background estimate
- Size filtering removes most non-stellar artifacts
- Local annulus background accounts for spatial variations in sky level

## Output Files Summary

| File | Description                     |
|------|---------------------------------|
| `combined_F336W.fits` | Median-stacked UV image         |
| `combined_F555W.fits` | Median-stacked visible image    |
| `combined_images.png` | Combined images visualization   |
| `matched_stars.txt` | Cross-matched stars catalogue   |
| `star_detections.png` | Detection visualization         |
| `photometry_catalog.txt` | Full photometry with magnitudes |
| `hr_diagram.png` | Color-magnitude diagram         |

## References

### Code Implementation Resources

1. SciPy ndimage module: https://docs.scipy.org/doc/scipy/reference/ndimage.html
2. Astropy FITS I/O: https://docs.astropy.org/en/stable/io/fits/
3. NumPy statistical functions: https://numpy.org/doc/stable/reference/routines.statistics.html
4. Matplotlib visualization: https://matplotlib.org/stable/gallery/index.html

## Author
**Joshua Rodriguez** for:
Space Detector Laboratory - Optical Detectors Assignment

*Data provided by Professor Andrew Kirwan *