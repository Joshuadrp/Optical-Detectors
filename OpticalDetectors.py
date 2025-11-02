import os
from astropy.io import fits
from astropy.visualization import simple_norm
import numpy as np
import matplotlib.pyplot as plt

"""
TASK 1
COSMIC RAY REMOVAL
"""

def cosmic_ray_removal(directory, output_file=None):
    """
    This function combines our fits images into a single image.
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

"""
Visualisation
"""

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

axes[0].imshow(f336w, cmap='gray', origin='lower', vmin=np.percentile(f336w, 1), vmax=np.percentile(f336w, 99))
axes[1].imshow(f555w, cmap='gray', origin='lower', vmin=np.percentile(f555w, 1), vmax=np.percentile(f555w, 99))
axes[0].set_title("Combined F336W")
axes[1].set_title("Combined F555W")

for ax in axes:
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.show()