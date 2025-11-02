import os
from astropy.io import fits
from astropy.visualization import simple_norm
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage

"""
TASK 1
COSMIC RAY REMOVAL
"""

def cosmic_ray_removal(directory, output_file=None):
    """
    This function combines our fits images into a single image,
    which as seen in class can remove most of the cosmic rays.
    We will have two images, one for each filter.
    """

    fits_files = [
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.lower().endswith('.fits')
    ]

    image_data = []
    for file in fits_files:
        with fits.open(file) as hdul:
            for hdu in hdul:
                if hdu.data is not None:
                    data = hdu.data.astype(np.float32)
                    image_data.append(data)
                    break
            else:
                raise ValueError(f"No image data found in {file}")

    combined_image = np.median(image_data, axis=0)

    if output_file is None:
        folder_name = os.path.basename(os.path.normpath(directory))
        output_file = f"combined_{folder_name}.fits"

    # Creates the combined image.
    fits.writeto(output_file, combined_image, overwrite=True)
    print(f"Combined image saved as {output_file}.")

    return combined_image

f336w = cosmic_ray_removal("data/F336W")
f555w = cosmic_ray_removal("data/F555W")

# """
# Visualisation TASK 1
# """
# fig, axes = plt.subplots(1, 2, figsize=(10, 5))
#
# axes[0].imshow(f336w, cmap='gray', origin='lower', vmin=np.percentile(f336w, 1), vmax=np.percentile(f336w, 99))
# axes[1].imshow(f555w, cmap='gray', origin='lower', vmin=np.percentile(f555w, 1), vmax=np.percentile(f555w, 99))
# axes[0].set_title("Combined F336W")
# axes[1].set_title("Combined F555W")
#
# for ax in axes:
#     ax.set_xticks([])
#     ax.set_yticks([])
#
# plt.tight_layout()
# plt.show()

"""
TASK 2
STAR FINDING and STAR CATALOGUE
"""

def star_finding(image, threshold_sigma=0.5, edge_buffer=50):
    """
    Star detection based on brightness for now, avoiding
    detection on the edge of the image.

    -Background Substraction
    -Filter by brightness
    -Avoid left and bottom edges
    """
    background = np.median(image)
    background_std = np.std(image)
    background_subtracted = image - background

    #Filter by brightness with subtracted background
    threshold = threshold_sigma * background_std
    source_mask = background_subtracted > threshold

    labeled_image, num_sources = ndimage.label(source_mask)
    print(f"Initial detections: {num_sources}")

    # Extract sources and exclude left and bottom edge
    sources = []

    for source_id in range(1, num_sources + 1):
        source_pixels = labeled_image == source_id
        area = np.sum(source_pixels)

        y_coords, x_coords = np.where(source_pixels)
        x_center = np.mean(x_coords)
        y_center = np.mean(y_coords)

        if (x_center > edge_buffer and
                edge_buffer < y_center):
            sources.append({
                'id': len(sources) + 1,
                'x': x_center,
                'y': y_center,
                'area': area
            })

    print(f"Sources after edge filtering: {len(sources)}")
    return sources

print("DETECTING SOURCES IN F336W")
sources_f336w = star_finding(f336w)

print("DETECTING SOURCES IN F555W")
sources_f555w = star_finding(f555w)

def cross_match_sources(src1,src2, match_radius=2.5):
    """
    Cross-match sources between F336W and F555W.

    -src1: list of sources from fist filter
    -src2: list of sources from second filter
    -match_radius: radius of cross-match
    """
    matched_stars = []
    for s1 in src1:
        # Find the closest source in filter 2
        min_distance = float('inf')
        closest_match = None

        for s2 in src2:
            # Calculate distance between sources
            distance = np.sqrt((s1['x'] - s2['x']) ** 2 + (s1['y'] - s2['y']) ** 2)

            if distance < min_distance:
                min_distance = distance
                closest_match = s2

        # If closest match is within radius, is a match.
        if closest_match is not None and min_distance <= match_radius:
            matched_stars.append({
                'filter1': s1,
                'filter2': closest_match,
                'distance': min_distance
            })

    print(f"  Matched sources: {len(matched_stars)}")
    return matched_stars

# Matched stars from F336W and F555W
matched = cross_match_sources(sources_f336w, sources_f555w, match_radius=3.0)


def star_catalog(matched_sources, output_file='star_catalog.txt'):
    """
    Create a final catalog from matched sources.

    Parameters:
    - matched_sources: list of matched stars
    - output_file: filename to save catalog
    """
    catalog = []

    for i, match in enumerate(matched_sources):
        # Average the coordinates from both filters
        x_avg = (match['filter1']['x'] + match['filter2']['x']) / 2.0
        y_avg = (match['filter1']['y'] + match['filter2']['y']) / 2.0

        catalog.append({
            'Object ID': i + 1,
            'x-center': x_avg,
            'y-center': y_avg
        })

    # Save to file
    with open(output_file, 'w') as f:
        f.write("Object ID\tx-center\ty-center\n")

        # Write each source
        for s in catalog:
            f.write(f"{s['Object ID']}\t{s['x-center']:.2f}\t{s['y-center']:.2f}\n")

    return catalog

# Create the final catalog
final_catalog = star_catalog(matched)

"""
TASK 3
PHOTOMETRY CATALOG
"""


def aperture_photometry(image, x, y, aperture_radius=5, inner_annulus=8, outer_annulus=12):
    """
    Perform circular aperture photometry on a source.

    Parameters:
    - image: 2D array of the image
    - x, y: center coordinates of the source
    - aperture_radius: radius of circular aperture for source flux
    - inner_annulus: inner radius for background annulus
    - outer_annulus: outer radius for background annulus

    Returns:
    - flux: background-subtracted flux
    - background: estimated local background per pixel
    """
    height, width = image.shape

    # Create coordinate grids
    y_grid, x_grid = np.ogrid[:height, :width]

    # Calculate distance from source center
    distance = np.sqrt((x_grid - x) ** 2 + (y_grid - y) ** 2)

    # Define aperture and annulus masks
    aperture_mask = distance <= aperture_radius
    annulus_mask = (distance >= inner_annulus) & (distance <= outer_annulus)

    # Measure background in annulus
    if np.sum(annulus_mask) > 0:
        background_pixels = image[annulus_mask]
        background_per_pixel = np.median(background_pixels)
    else:
        background_per_pixel = 0

    # Measure flux in aperture
    aperture_pixels = image[aperture_mask]
    total_flux = np.sum(aperture_pixels)

    # Subtract background contribution
    n_aperture_pixels = np.sum(aperture_mask)
    background_contribution = background_per_pixel * n_aperture_pixels
    flux = total_flux - background_contribution

    return flux, background_per_pixel


def flux_to_magnitude(flux, zero_point):
    """
    Convert flux to magnitude with zero-point calibration.
    Handles negative/zero flux by returning NaN.

    Parameters:
    - flux: measured flux
    - zero_point: calibration constant (typically 25-26 for HST)
    """
    if flux <= 0:
        return np.nan
    return -2.5 * np.log10(flux) + zero_point


def perform_photometry_catalog(catalog, image_f336w, image_f555w, aperture_radius=5, reference_mag=15.0):
    """
    Perform photometry on all sources in catalog for both filters.

    Parameters:
    - catalog: list of sources with x-center, y-center
    - image_f336w: F336W image array
    - image_f555w: F555W image array
    - aperture_radius: radius for aperture photometry
    - reference_mag: magnitude to assign to the brightest star (default 15.0)

    Returns:
    - photometry_catalog: catalog with magnitudes added
    """
    photometry_catalog = []
    fluxes_f336w = []
    fluxes_f555w = []

    # First pass: measure all fluxes
    for source in catalog:
        x = source['x-center']
        y = source['y-center']

        # Perform photometry in both filters
        flux_f336w, bg_f336w = aperture_photometry(image_f336w, x, y, aperture_radius)
        flux_f555w, bg_f555w = aperture_photometry(image_f555w, x, y, aperture_radius)

        fluxes_f336w.append(flux_f336w)
        fluxes_f555w.append(flux_f555w)

    # Find the brightest star (maximum flux) in each filter
    max_flux_f336w = max(fluxes_f336w)
    max_flux_f555w = max(fluxes_f555w)

    # Calculate zero-points based on brightest star
    zero_point_f336w = reference_mag + 2.5 * np.log10(max_flux_f336w)
    zero_point_f555w = reference_mag + 2.5 * np.log10(max_flux_f555w)

    print(f"Brightest star flux F336W: {max_flux_f336w:.2f} → zero-point: {zero_point_f336w:.2f}")
    print(f"Brightest star flux F555W: {max_flux_f555w:.2f} → zero-point: {zero_point_f555w:.2f}")

    # Second pass: convert all fluxes to magnitudes
    for i, source in enumerate(catalog):
        mag_f336w = flux_to_magnitude(fluxes_f336w[i], zero_point_f336w)
        mag_f555w = flux_to_magnitude(fluxes_f555w[i], zero_point_f555w)

        photometry_catalog.append({
            'Object ID': source['Object ID'],
            'x-center': source['x-center'],
            'y-center': source['y-center'],
            'aperture_radius': aperture_radius,
            'mag_F336W': mag_f336w,
            'mag_F555W': mag_f555w
        })

    return photometry_catalog


# Perform photometry on final catalog
print(f"\n{'=' * 50}")
print("PERFORMING APERTURE PHOTOMETRY")
print(f"{'=' * 50}")

# Use background-subtracted images for photometry
f336w_bgsub = f336w - np.median(f336w)
f555w_bgsub = f555w - np.median(f555w)

photometry_catalog = perform_photometry_catalog(final_catalog, f336w_bgsub, f555w_bgsub)

# """
# debug
# """
# Debug: Check a few flux values
# print("\nDEBUG - First 5 sources:")
# for i, source in enumerate(photometry_catalog[:5]):
#     x, y = source['x-center'], source['y-center']
#     flux_f336w, bg = aperture_photometry(f336w_bgsub, x, y, aperture_radius=5)
#     print(f"Source {i+1}: flux_F336W = {flux_f336w:.2f}, mag = {source['mag_F336W']:.3f}")

# Save photometry catalog
with open('photometry_catalog.txt', 'w') as f:
    f.write("Object_ID\tx-center\ty-center\taperture_radius\tmag_F336W\tmag_F555W\n")
    for source in photometry_catalog:
        f.write(f"{source['Object ID']}\t{source['x-center']:.2f}\t{source['y-center']:.2f}\t"
                f"{source['aperture_radius']:.1f}\t{source['mag_F336W']:.3f}\t{source['mag_F555W']:.3f}\n")

print(f"Photometry catalog saved: photometry_catalog.txt")
print(f"Total sources with photometry: {len(photometry_catalog)}")


# """
# Visualization of matched stars between filters.
# """
#
# fig, axes = plt.subplots(1, 2, figsize=(16, 7))
#
# # F336W with final catalog sources
# axes[0].imshow(f336w, cmap='gray', origin='lower',
#                vmin=np.percentile(f336w, 1), vmax=np.percentile(f336w, 99))
# for source in final_catalog:
#     axes[0].plot(source['x-center'], source['y-center'], 'r+', markersize=4, markeredgewidth=0.8)
# axes[0].set_title(f"F336W - Final Catalog: {len(final_catalog)} matched sources")
# axes[0].set_xlabel("X pixel")
# axes[0].set_ylabel("Y pixel")
#
# # F555W with final catalog sources
# axes[1].imshow(f555w, cmap='gray', origin='lower',
#                vmin=np.percentile(f555w, 1), vmax=np.percentile(f555w, 99))
# for source in final_catalog:
#     axes[1].plot(source['x-center'], source['y-center'], 'r+', markersize=4, markeredgewidth=0.8)
# axes[1].set_title(f"F555W - Final Catalog: {len(final_catalog)} matched sources")
# axes[1].set_xlabel("X pixel")
# axes[1].set_ylabel("Y pixel")
#
# plt.tight_layout()
# plt.show()






# """
# VISUALIZATION
# """
# fig, axes = plt.subplots(1, 2, figsize=(14, 6))
#
# # F336W with detections
# axes[0].imshow(f336w, cmap='gray', origin='lower',
#                vmin=np.percentile(f336w, 1), vmax=np.percentile(f336w, 99))
# for source in sources_f336w:
#     axes[0].plot(source['x'], source['y'], 'r+', markersize=5, markeredgewidth=0.5)
# axes[0].set_title(f"F336W - {len(sources_f336w)} sources detected")
# axes[0].set_xlabel("X pixel")
# axes[0].set_ylabel("Y pixel")
#
# # F555W with detections
# axes[1].imshow(f555w, cmap='gray', origin='lower',
#                vmin=np.percentile(f555w, 1), vmax=np.percentile(f555w, 99))
# for source in sources_f555w:
#     axes[1].plot(source['x'], source['y'], 'r+', markersize=5, markeredgewidth=0.5)
# axes[1].set_title(f"F555W - {len(sources_f555w)} sources detected")
# axes[1].set_xlabel("X pixel")
# axes[1].set_ylabel("Y pixel")
#
# plt.tight_layout()
# plt.show()


