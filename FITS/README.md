# FITS Files - Astronomy and Space Science Data

FITS (Flexible Image Transport System) is the standard format for astronomical data. Developed in the 1970s, it's used by NASA, ESA, and observatories worldwide for storing images, spectra, and tables from telescopes and space missions.

## Available Datasets

### hubble_wfpc2_sample.fits
**Size:** 684 KB
**Instrument:** Wide Field Planetary Camera 2 (WFPC2)
**Telescope:** Hubble Space Telescope
**Source:** NASA FITS samples

Real astronomical image from Hubble's WFPC2 camera. WFPC2 was Hubble's workhorse camera from 1993-2009, capturing iconic images like the Pillars of Creation and deep field observations.

### hubble_foc_sample.fits
**Size:** 4.1 MB
**Instrument:** Faint Object Camera (FOC)
**Telescope:** Hubble Space Telescope
**Source:** NASA FITS samples

Astronomical observations from Hubble's FOC, one of the original instruments. FOC was designed for high-resolution imaging of faint objects and operated from 1990-2002.

## About FITS Format

FITS files consist of:
- **Header:** ASCII metadata with keywords (SIMPLE, BITPIX, NAXIS, etc.)
- **Data:** N-dimensional arrays (images, spectra, cubes)
- **Extensions:** Multiple HDUs (Header-Data Units) in one file

**Key Features:**
- Platform-independent binary format
- Self-describing with extensive metadata
- Support for images, tables, and multi-dimensional data
- World Coordinate System (WCS) for sky coordinates
- FITS standard maintained by IAU

## Installation

```bash
# Install astropy (recommended)
pip install astropy

# Or install with extra tools
pip install astropy[all]

# Install fitsio (faster I/O)
pip install fitsio

# Command-line tools
pip install astropy  # includes fitsheader, fitsinfo
```

## Usage Examples

### Python (Astropy - Recommended)
```python
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# Open FITS file
hdul = fits.open('hubble_wfpc2_sample.fits')

# Print file information
hdul.info()

# Access header (first HDU)
header = hdul[0].header
print("\nHeader keywords:")
for key in ['TELESCOP', 'INSTRUME', 'DATE-OBS', 'EXPTIME']:
    if key in header:
        print(f"{key}: {header[key]}")

# Access image data
data = hdul[0].data

if data is not None:
    print(f"\nImage shape: {data.shape}")
    print(f"Data type: {data.dtype}")
    print(f"Min/Max: {np.min(data):.2f} / {np.max(data):.2f}")
    print(f"Mean: {np.mean(data):.2f}")

    # Display image
    plt.figure(figsize=(10, 10))
    plt.imshow(data, cmap='gray', origin='lower', vmin=np.percentile(data, 1), vmax=np.percentile(data, 99))
    plt.colorbar(label='Pixel Value')
    plt.title(f"{header.get('OBJECT', 'Hubble Image')}")
    plt.xlabel('X Pixel')
    plt.ylabel('Y Pixel')
    plt.savefig('hubble_image.png', dpi=150)

hdul.close()
```

### World Coordinate System (WCS)
```python
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

# Open FITS with WCS
hdul = fits.open('hubble_wfpc2_sample.fits')
wcs = WCS(hdul[0].header)
data = hdul[0].data

# Plot with sky coordinates
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection=wcs)

ax.imshow(data, cmap='viridis', origin='lower')
ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
ax.set_title('Hubble Image with WCS')
ax.grid(color='white', ls='dotted')

plt.savefig('hubble_wcs.png', dpi=150)
```

### Read/Write FITS Files
```python
from astropy.io import fits
import numpy as np

# Create new FITS file
data = np.random.random((100, 100))

hdu = fits.PrimaryHDU(data)
hdu.header['OBJECT'] = 'Test Image'
hdu.header['TELESCOP'] = 'Simulated'
hdu.header['INSTRUME'] = 'Python'

hdu.writeto('new_image.fits', overwrite=True)

# Append extension
new_data = np.random.random((50, 50))
hdu2 = fits.ImageHDU(new_data, name='EXTENSION1')

hdul = fits.open('new_image.fits', mode='append')
hdul.append(hdu2)
hdul.flush()
hdul.close()
```

### Extract Specific Extension
```python
from astropy.io import fits

# Open and access specific extension
with fits.open('hubble_foc_sample.fits') as hdul:
    # Get specific extension
    if len(hdul) > 1:
        ext1 = hdul[1]
        print(f"Extension name: {ext1.name}")
        print(f"Extension shape: {ext1.data.shape if ext1.data is not None else 'No data'}")

    # Or by name
    try:
        sci = hdul['SCI']  # Science data
        err = hdul['ERR']  # Error array
        dq = hdul['DQ']    # Data quality
    except KeyError:
        print("Extension not found")
```

### FITS Tables
```python
from astropy.io import fits
from astropy.table import Table

# Read FITS table
table = Table.read('table.fits')
print(table)

# Or with fits module
with fits.open('table.fits') as hdul:
    data = hdul[1].data  # Tables usually in extension 1
    print(f"Columns: {data.columns.names}")
    print(f"Rows: {len(data)}")

    # Access columns
    col1 = data['COLUMN_NAME']
```

### Python (fitsio - Faster I/O)
```python
import fitsio

# Read data
data = fitsio.read('hubble_wfpc2_sample.fits')
print(f"Shape: {data.shape}")

# Read header
header = fitsio.read_header('hubble_wfpc2_sample.fits')
print(header)

# Read specific HDU
data_ext = fitsio.read('hubble_foc_sample.fits', ext=1)

# Write FITS
fitsio.write('output.fits', data, header=header, clobber=True)
```

### Command Line Tools
```bash
# Print file info
fitsinfo hubble_wfpc2_sample.fits

# Print header
fitsheader hubble_wfpc2_sample.fits

# Print all headers
fitsheader hubble_wfpc2_sample.fits -e 0 -e 1 -e 2

# Convert to image
# (requires additional tools like fits2png)
```

## Image Processing Examples

### Stacking Multiple Images
```python
from astropy.io import fits
import numpy as np

# Load multiple images
images = []
for filename in ['image1.fits', 'image2.fits', 'image3.fits']:
    with fits.open(filename) as hdul:
        images.append(hdul[0].data)

# Stack (median combine)
stacked = np.median(images, axis=0)

# Save result
fits.writeto('stacked.fits', stacked, overwrite=True)
```

### Photometry
```python
from astropy.io import fits
from photutils import CircularAperture, aperture_photometry
import numpy as np

# Load image
with fits.open('hubble_wfpc2_sample.fits') as hdul:
    data = hdul[0].data

# Define apertures (x, y positions)
positions = [(100, 100), (200, 200)]
apertures = CircularAperture(positions, r=10)

# Perform aperture photometry
phot_table = aperture_photometry(data, apertures)
print(phot_table)
```

### Background Subtraction
```python
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np

with fits.open('hubble_wfpc2_sample.fits') as hdul:
    data = hdul[0].data

# Estimate background
mean, median, std = sigma_clipped_stats(data, sigma=3.0)
print(f"Background: {median:.2f} +/- {std:.2f}")

# Subtract background
data_sub = data - median

# Save
fits.writeto('background_subtracted.fits', data_sub, overwrite=True)
```

## Convert FITS to Other Formats

```python
from astropy.io import fits
import numpy as np
from PIL import Image

# FITS to PNG
with fits.open('hubble_wfpc2_sample.fits') as hdul:
    data = hdul[0].data

    # Normalize to 0-255
    data_norm = (data - np.min(data)) / (np.max(data) - np.min(data))
    data_8bit = (data_norm * 255).astype(np.uint8)

    # Save as PNG
    img = Image.fromarray(data_8bit)
    img.save('hubble.png')

# FITS to HDF5
import h5py
with fits.open('hubble_wfpc2_sample.fits') as hdul:
    with h5py.File('hubble.h5', 'w') as h5f:
        h5f.create_dataset('image', data=hdul[0].data)

        # Save header as attributes
        for key, value in hdul[0].header.items():
            try:
                h5f['image'].attrs[key] = value
            except:
                pass  # Skip problematic keys
```

## Common FITS Keywords

- **SIMPLE:** True for standard FITS
- **BITPIX:** Bits per pixel (-64=double, -32=float, 16=short, 32=int)
- **NAXIS:** Number of dimensions
- **NAXIS1, NAXIS2:** Size of dimensions
- **TELESCOP:** Telescope name
- **INSTRUME:** Instrument name
- **OBJECT:** Object observed
- **DATE-OBS:** Observation date
- **EXPTIME:** Exposure time (seconds)
- **CRPIX1, CRPIX2:** Reference pixel for WCS
- **CRVAL1, CRVAL2:** Sky coordinates at reference pixel
- **CD1_1, CD2_2:** Pixel scale (degrees/pixel)

## Applications

FITS is used for:
- **Optical astronomy:** Hubble, VLT, Keck, JWST
- **Radio astronomy:** VLA, ALMA, LOFAR
- **X-ray/Gamma-ray:** Chandra, XMM-Newton, Fermi
- **Planetary science:** Mars rovers, planetary missions
- **Simulation data:** Astrophysical simulations
- **Time-series:** Light curves, spectra

## Major Telescopes Using FITS

- Hubble Space Telescope (HST)
- James Webb Space Telescope (JWST)
- Chandra X-ray Observatory
- Spitzer Space Telescope
- Very Large Telescope (VLT)
- Keck Observatory
- ALMA
- And virtually all professional observatories

## Resources

- **FITS Standard:** [fits.gsfc.nasa.gov](https://fits.gsfc.nasa.gov/)
- **Astropy FITS:** [docs.astropy.org/en/stable/io/fits](https://docs.astropy.org/en/stable/io/fits/)
- **FITS Primer:** [fits.gsfc.nasa.gov/fits_primer.html](https://fits.gsfc.nasa.gov/fits_primer.html)
- **NASA FITS Support:** [heasarc.gsfc.nasa.gov/fits](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html)
