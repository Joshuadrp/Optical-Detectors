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

def star_finding(image, threshold_sigma=1, edge_buffer=50):
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


# Test on both filters
print("=" * 50)
print("DETECTING SOURCES IN F336W")
print("=" * 50)
sources_f336w = star_finding(f336w)

print("\n" + "=" * 50)
print("DETECTING SOURCES IN F555W")
print("=" * 50)
sources_f555w = star_finding(f555w)

def cross_match_sources(src1,src2, match_radius=3):
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

    print(f"\nCross-matching results (radius={match_radius} pixels):")
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
        for source in catalog:
            f.write(f"{source['Object ID']}\t{source['x-center']:.2f}\t{source['y-center']:.2f}\n")

    return catalog

# Create the final catalog
final_catalog = star_catalog(matched)


"""
VISUALIZATION
"""
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


