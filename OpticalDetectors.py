import os
import numpy as np
from astropy.io import fits
from scipy import ndimage
import matplotlib.pyplot as plt

# ============================================================
# TASK 1: COSMIC RAY REMOVAL
# ============================================================

def combine_fits_images(directory, output_file=None):
    """Combine FITS images using median to remove cosmic rays."""
    fits_files = [os.path.join(directory, f) for f in os.listdir(directory)
                  if f.lower().endswith('.fits')]

    image_data = []
    for file in fits_files:
        with fits.open(file) as hdul:
            for hdu in hdul:
                if hdu.data is not None:
                    image_data.append(hdu.data.astype(np.float32))
                    break

    combined = np.median(image_data, axis=0)

    if output_file is None:
        folder_name = os.path.basename(os.path.normpath(directory))
        output_file = f"combined_{folder_name}.fits"

    fits.writeto(output_file, combined, overwrite=True)
    print("OK")
    return combined


def show_combined_images(img1, img2):
    """Plot two combined images side by side."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    imgs = [img1, img2]
    titles = ["F336W (UV)", "F555W (Visible)"]

    for ax, img, title in zip(axes, imgs, titles):
        vmin, vmax = np.percentile(img, [1, 99])
        ax.imshow(img, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis('off')

    plt.tight_layout()
    plt.savefig("combined_images.png", dpi=150)
    plt.show()


# ============================================================
# TASK 2: DETECT STARS, CROSS MATCHING AND CATALOGUE CREATION
# ============================================================
def find_sources(image, sigma=1, min_area=3, max_area=200, edge=50):
    """Detect stars in the merged FITS image using a simple threshold method."""

    #Estimate background with sigma clipping
    flat = image.flatten()
    for _ in range(3):  # repeat to remove outliers
        med = np.median(flat)
        std = np.std(flat)
        flat = flat[np.abs(flat - med) < 3 * std]

    bkg = np.median(flat)
    bkg_std = np.std(flat)
    thresh = sigma * bkg_std

    #Create mask for bright pixels
    mask = (image - bkg) > thresh
    labeled, n = ndimage.label(mask)

    sources = []
    for i in range(1, n + 1):
        pix = (labeled == i)
        area = np.sum(pix)
        if not (min_area <= area <= max_area):
            continue  # skip too small or too large

        y, x = np.where(pix)
        xc, yc = np.mean(x), np.mean(y)

        # skip sources near edges
        if xc > edge and yc > edge:
            sources.append({'x': xc, 'y': yc, 'area': area})

    return sources


def cross_match(src1, src2, radius=2.0):
    """Cross-match sources between two filters (e.g. UV and visible)."""
    matches = []
    for s1 in src1:
        #compute distances from this source to all in the other list
        dists = [np.hypot(s1['x'] - s2['x'], s1['y'] - s2['y']) for s2 in src2]
        if not dists:
            continue

        dmin = min(dists)  # find nearest neighbor distance
        if dmin <= radius:  # if close enough, count as a match
            s2 = src2[dists.index(dmin)]
            # average positions from both images
            matches.append({'x': (s1['x'] + s2['x']) / 2,
                            'y': (s1['y'] + s2['y']) / 2})

    print(f"  Matched {len(matches)} sources")
    return matches


def show_detections(img1, img2, src1, src2, matched):
    """Display detections and matched stars."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    data = [
        (img1, src1, f"F336W detections ({len(src1)})"),
        (img2, src2, f"F555W detections ({len(src2)})"),
        (img1, matched, f"Matched stars ({len(matched)})")
    ]

    for ax, (img, src, title) in zip(axes, data):
        vmin, vmax = np.percentile(img, [1, 99])
        ax.imshow(img, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        for s in src:
            ax.plot(s['x'], s['y'], 'r+', markersize=4)
        ax.set_title(title)
        ax.axis('off')

    plt.tight_layout()
    plt.savefig("star_detections.png", dpi=150)
    plt.show()


def save_catalog(sources, filename, cols=None):
    """Save source list to a tab-separated text file."""
    if cols is None:
        cols = ["ID", "x", "y"]
    with open(filename, "w") as f:
        f.write("\t".join(cols) + "\n")
        for i, s in enumerate(sources, 1):
            row = [f"{i}", f"{s['x']:.2f}", f"{s['y']:.2f}"]
            if "mag_F336W" in s:
                row += [f"{s['ap_radius']:.1f}",
                        f"{s['mag_F336W']:.3f}",
                        f"{s['mag_F555W']:.3f}",
                        f"{s['color']:.3f}"]
            f.write("\t".join(row) + "\n")
    print(f"[OK] Catalog saved: {filename}")


# ============================================================
# TASK 3: APERTURE PHOTOMETRY AND CATALOG CREATION
# ============================================================

def aperture_flux(image, x, y, r=5, r_in=8, r_out=12):
    """Return aperture flux after local background subtraction."""
    h, w = image.shape
    yy, xx = np.ogrid[:h, :w]
    dist = np.sqrt((xx - x) ** 2 + (yy - y) ** 2)

    # Define circular regions
    ap_mask = dist <= r
    ann_mask = (dist >= r_in) & (dist <= r_out)

    # Measure background in annulus, subtract from aperture
    bkg = np.median(image[ann_mask]) if np.any(ann_mask) else 0
    flux = np.sum(image[ap_mask]) - bkg * np.sum(ap_mask)
    return flux


def measure_mags(sources, img1, img2, ref_mag=15.0):
    """Compute magnitudes for matched stars."""
    f1 = [aperture_flux(img1, s['x'], s['y']) for s in sources]
    f2 = [aperture_flux(img2, s['x'], s['y']) for s in sources]

    # Calibrate using brightest star from filters
    zp1 = ref_mag + 2.5 * np.log10(max(f1))
    zp2 = ref_mag + 2.5 * np.log10(max(f2))
    print(f"Zero-points: F336W={zp1:.2f}, F555W={zp2:.2f}")

    catalog = []
    for s, flux1, flux2 in zip(sources, f1, f2):
        if flux1 > 0 and flux2 > 0:
            mag1 = -2.5 * np.log10(flux1) + zp1
            mag2 = -2.5 * np.log10(flux2) + zp2
            catalog.append({
                'x': s['x'],
                'y': s['y'],
                'ap_radius': 5,
                'mag_F336W': mag1,
                'mag_F555W': mag2,
                'color': mag1 - mag2
            })

    return catalog


# ============================================================
# TASK 4: HERTZSPRUNG–RUSSELL DIAGRAM
# ============================================================

def plot_hr(catalog):
    """Plot the HR diagram (color–magnitude diagram)."""
    color = [s["color"] for s in catalog]
    mag = [s["mag_F336W"] for s in catalog]

    plt.figure(figsize=(8, 10))
    plt.scatter(color, mag, s=15, c='gold', edgecolors='black', alpha=0.7)
    plt.gca().invert_yaxis()

    plt.xlabel("Color (F336W - F555W)")
    plt.ylabel("Magnitude (F336W)")
    plt.title("Hertzsprung–Russell Diagram")

    plt.grid(alpha=0.3, linestyle='--')
    plt.text(0.05, 0.95, f"N = {len(catalog)}", transform=plt.gca().transAxes,
             ha='left', va='top', fontsize=9, bbox=dict(facecolor='wheat', alpha=0.5))
    plt.tight_layout()
    plt.savefig("hr_diagram.png", dpi=150)
    plt.show()


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("\n --- COSMIC RAY REMOVAL --- ")
    f336w = combine_fits_images("data/F336W")
    f555w = combine_fits_images("data/F555W")
    show_combined_images(f336w, f555w)

    print("\n --- STAR DETECTION & MATCHING --- ")
    src336 = find_sources(f336w)
    src555 = find_sources(f555w)
    print("Matching...")
    matched = cross_match(src336, src555)
    show_detections(f336w, f555w, src336, src555, matched)
    save_catalog(matched, "matched_stars.txt")

    print("\n --- APERTURE PHOTOMETRY --- ")
    f336w -= np.median(f336w)
    f555w -= np.median(f555w)
    catalog = measure_mags(matched, f336w, f555w)
    save_catalog(catalog, "photometry_catalog.txt",
                 cols=["ID", "x", "y", "ap_radius", "mag_F336W", "mag_F555W", "color"])

    print("\n --- HERTZSPRUNG–RUSSELL DIAGRAM --- ")
    plot_hr(catalog)

    print("\n Complete ✅")
