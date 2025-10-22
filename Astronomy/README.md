# Astronomy Datasets

This directory contains astronomical datasets in FITS and HDF5 formats for telescope data processing and analysis.

## Datasets Overview

Total size: **2.7 MB**

| Dataset | Format | Size | Description |
|---------|--------|------|-------------|
| Galaxy Image | FITS | 2.1 MB | 512×512 CCD image of elliptical galaxy |
| Stellar Spectrum | FITS | 57 KB | Optical spectrum (4000-7000 Å) |
| Variable Star Light Curves | HDF5 | 489 KB | 10 stars × 3 photometric bands |
| Exoplanet Transits | HDF5 | 83 KB | 5 planetary transit observations |

---

## 1. Galaxy Image (FITS)

**File:** `galaxy_ngc1234_r_band.fits` (2.1 MB)

Simulated CCD observation of an elliptical galaxy in R-band filter.

**Image Properties:**
- **Dimensions:** 512 × 512 pixels
- **Pixel scale:** 0.396 arcsec/pixel
- **Exposure time:** 300 seconds
- **Filter:** R (red)
- **Units:** ADU (Analog-to-Digital Units)

**Galaxy Properties:**
- **Type:** Elliptical galaxy (Sersic n=4 profile)
- **Effective radius:** 40 pixels (~16 arcseconds)
- **Axis ratio:** 0.7 (b/a)
- **Position angle:** 30°
- **RA/Dec:** 123.456°, +45.678°

**Additional Features:**
- 50 point sources (stars) with Gaussian PSF
- Poisson sky background (~100 ADU)
- Read noise (σ = 5 ADU)
- PSF FWHM: ~2 pixels (Gaussian σ = 2.0)

**FITS Header:**
```
TELESCOP = 'Example Observatory'
INSTRUME = 'CCD Camera'
OBJECT   = 'NGC1234'
FILTER   = 'R'
EXPTIME  = 300.0
AIRMASS  = 1.15
GAIN     = 1.5
RDNOISE  = 5.0
PIXSCALE = 0.396
```

**Usage Example (Python):**
```python
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# Read FITS file
hdu = fits.open('galaxy_ngc1234_r_band.fits')
data = hdu[0].data
header = hdu[0].header

# Display image
plt.figure(figsize=(10, 10))
plt.imshow(data, cmap='gray', vmin=90, vmax=200, origin='lower')
plt.colorbar(label='Counts (ADU)')
plt.title(f"{header['OBJECT']} - {header['FILTER']} band")
plt.xlabel('X (pixels)')
plt.ylabel('Y (pixels)')
plt.show()

# Compute image statistics
print(f"Mean: {np.mean(data):.1f} ADU")
print(f"Std dev: {np.std(data):.1f} ADU")
print(f"Max: {np.max(data):.1f} ADU")
```

**Usage Example (DS9 visualization):**
```bash
# Open in DS9
ds9 galaxy_ngc1234_r_band.fits &

# Or use fv (FITS viewer)
fv galaxy_ngc1234_r_band.fits &
```

**Photometry Example (Python):**
```python
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry
from astropy.stats import sigma_clipped_stats

# Background estimation
mean, median, std = sigma_clipped_stats(data, sigma=3.0)

# Source detection
daofind = DAOStarFinder(fwhm=3.0, threshold=5*std)
sources = daofind(data - median)
print(f"Found {len(sources)} sources")

# Aperture photometry
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=5.0)
phot_table = aperture_photometry(data - median, apertures)
print(phot_table)
```

---

## 2. Stellar Spectrum (FITS)

**File:** `star_hd12345_spectrum.fits` (57 KB)

Optical spectrum of a solar-type star (T_eff ~ 6000 K).

**Spectrum Properties:**
- **Wavelength range:** 4000 - 7000 Angstroms
- **Spectral pixels:** 2048
- **Dispersion:** ~1.46 Å/pixel
- **Resolving power:** R ~ 5000
- **SNR:** ~100

**Prominent Features:**
- H-alpha absorption (6563 Å)
- H-beta absorption (4861 Å)
- Mg I triplet (~5175 Å)
- Na D doublet (5890, 5896 Å)
- Blackbody-like continuum

**FITS Structure:**
Binary table with 3 columns:
- `WAVELENGTH` (Angstrom): Wavelength array
- `FLUX` (erg/s/cm²/Å): Flux values
- `ERROR` (erg/s/cm²/Å): 1-sigma uncertainties

**Usage Example (Python):**
```python
from astropy.io import fits
import matplotlib.pyplot as plt

# Read spectrum
hdu = fits.open('star_hd12345_spectrum.fits')
data = hdu[1].data

wavelength = data['WAVELENGTH']
flux = data['FLUX']
error = data['ERROR']

# Plot spectrum
plt.figure(figsize=(12, 6))
plt.plot(wavelength, flux, 'k-', linewidth=0.5, label='Observed')
plt.fill_between(wavelength, flux-error, flux+error, alpha=0.3, label='1σ error')
plt.xlabel('Wavelength (Å)')
plt.ylabel('Flux (erg/s/cm²/Å)')
plt.title('HD12345 Optical Spectrum')
plt.legend()
plt.grid(alpha=0.3)

# Mark spectral lines
lines = {
    'H-alpha': 6563,
    'H-beta': 4861,
    'Na D': 5893,
    'Mg I': 5175
}

for name, wave in lines.items():
    plt.axvline(wave, color='r', linestyle='--', alpha=0.5)
    plt.text(wave, plt.ylim()[1]*0.9, name, rotation=90, va='bottom')

plt.show()
```

**Radial Velocity Measurement:**
```python
from astropy.modeling import models, fitting
import numpy as np

# Find H-alpha line
ha_region = (wavelength > 6555) & (wavelength < 6570)
ha_wave = wavelength[ha_region]
ha_flux = flux[ha_region]

# Fit Gaussian absorption
g_init = models.Gaussian1D(amplitude=-200, mean=6563, stddev=2)
fit_g = fitting.LevMarOptimizer()
g_fitted = fit_g(g_init, ha_wave, ha_flux)

# Compute velocity shift
c = 299792.458  # km/s
rest_wavelength = 6562.8  # Rest-frame H-alpha
observed = g_fitted.mean.value
velocity = c * (observed - rest_wavelength) / rest_wavelength
print(f"Radial velocity: {velocity:.2f} km/s")
```

---

## 3. Variable Star Light Curves (HDF5)

**File:** `variable_star_lightcurves.h5` (489 KB)

Multi-band photometric time series for 10 variable stars.

**Contents:**
- **10 variable stars:**
  - Cepheid variables (P = 5-20 days)
  - RR Lyrae stars (P = 0.3-0.9 days)
  - Eclipsing binaries (P = 1-10 days)
- **3 photometric bands:** g, r, i
- **200-500 epochs** per star
- **Time span:** 1 year

**HDF5 Structure:**
```
/star_001/ through /star_010/
  ├── attributes: ra, dec, var_type, period, amplitude
  ├── g/ (g-band photometry)
  │   ├── time (Julian Date)
  │   ├── magnitude (AB magnitudes)
  │   └── mag_error
  ├── r/ (r-band)
  └── i/ (i-band)
```

**Usage Example (Python):**
```python
import h5py
import matplotlib.pyplot as plt
import numpy as np

# Read light curve
with h5py.File('variable_star_lightcurves.h5', 'r') as f:
    star = f['star_001']

    # Get star properties
    print(f"Variable type: {star.attrs['var_type']}")
    print(f"Period: {star.attrs['period']:.3f} days")
    print(f"Amplitude: {star.attrs['amplitude']:.3f} mag")

    # Plot multi-band light curve
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)

    for ax, band in zip(axes, ['g', 'r', 'i']):
        time = star[band]['time'][:]
        mag = star[band]['magnitude'][:]
        err = star[band]['mag_error'][:]

        ax.errorbar(time, mag, yerr=err, fmt='o', markersize=2, alpha=0.6)
        ax.set_ylabel(f'{band}-band (mag)')
        ax.invert_yaxis()  # Magnitudes increase downward
        ax.grid(alpha=0.3)

    axes[-1].set_xlabel('Time (JD)')
    plt.suptitle(f"{star.attrs['var_type']} - Star 001")
    plt.tight_layout()
    plt.show()
```

**Phase-Folded Light Curve:**
```python
with h5py.File('variable_star_lightcurves.h5', 'r') as f:
    star = f['star_001']
    period = star.attrs['period']

    time = star['r']['time'][:]
    mag = star['r']['magnitude'][:]

    # Compute phase
    phase = ((time - time[0]) % period) / period

    # Sort by phase
    sort_idx = np.argsort(phase)

    plt.figure(figsize=(10, 6))
    plt.plot(phase[sort_idx], mag[sort_idx], 'ko', markersize=3)
    plt.xlabel('Phase')
    plt.ylabel('r-band magnitude')
    plt.gca().invert_yaxis()
    plt.title(f'Phase-folded light curve (P = {period:.2f} d)')
    plt.grid(alpha=0.3)
    plt.show()
```

---

## 4. Exoplanet Transit Observations (HDF5)

**File:** `exoplanet_transits.h5` (83 KB)

High-precision photometry of transiting exoplanets.

**Contents:**
- **5 planetary systems**
- **500 data points** per transit
- **Planet radii:** 0.8 - 2.0 Earth radii
- **Orbital periods:** 2 - 20 days
- **Transit durations:** 2 - 6 hours
- **Photometric precision:** 100 ppm

**HDF5 Structure:**
```
/system_01/ through /system_05/
  ├── attributes: planet_radius, star_radius, orbital_period, transit_duration
  ├── time (hours from mid-transit)
  ├── flux (normalized, 1.0 = out-of-transit)
  └── flux_error
```

**Usage Example (Python):**
```python
import h5py
import matplotlib.pyplot as plt

# Read transit data
with h5py.File('exoplanet_transits.h5', 'r') as f:
    system = f['system_01']

    print(f"Planet radius: {system.attrs['planet_radius']:.2f} R_Earth")
    print(f"Orbital period: {system.attrs['orbital_period']:.2f} days")

    time = system['time'][:]
    flux = system['flux'][:]
    error = system['flux_error'][:]

    # Plot transit
    plt.figure(figsize=(12, 6))
    plt.errorbar(time, flux, yerr=error, fmt='k.', markersize=3, alpha=0.5)
    plt.xlabel('Time from mid-transit (hours)')
    plt.ylabel('Normalized Flux')
    plt.title(f"Exoplanet Transit - System 01")
    plt.grid(alpha=0.3)
    plt.axhline(1.0, color='r', linestyle='--', alpha=0.5, label='Out-of-transit')
    plt.legend()
    plt.show()

    # Measure transit depth
    in_transit = (time > -2) & (time < 2)
    out_transit = (np.abs(time) > 3)

    depth = 1.0 - flux[in_transit].mean()
    print(f"Transit depth: {depth*100:.3f}% ({depth*1e6:.0f} ppm)")
```

**Transit Fitting:**
```python
from scipy.optimize import curve_fit

def transit_model(time, depth, duration):
    """Simple box model with linear ingress/egress"""
    flux = np.ones_like(time)
    ingress = 0.5  # hours

    for i, t in enumerate(time):
        if abs(t) < (duration/2 - ingress):
            flux[i] = 1.0 - depth
        elif abs(t) < duration/2:
            frac = (duration/2 - abs(t)) / ingress
            flux[i] = 1.0 - depth * frac

    return flux

# Fit transit
popt, pcov = curve_fit(transit_model, time, flux,
                       p0=[0.01, 4.0], sigma=error)

depth_fit, duration_fit = popt
print(f"Fitted depth: {depth_fit*100:.3f}%")
print(f"Fitted duration: {duration_fit:.2f} hours")

# Plot fit
plt.plot(time, flux, 'k.', markersize=2, alpha=0.5, label='Data')
plt.plot(time, transit_model(time, *popt), 'r-', linewidth=2, label='Model')
plt.legend()
plt.show()
```

---

## File Formats

### FITS (Flexible Image Transport System)
- **Standard:** NASA/IAU FITS standard
- **Structure:** Header + Data Unit(s)
- **Data types:** Images, tables, spectra
- **Compression:** Optional (gzip, Rice)

### HDF5 (Hierarchical Data Format 5)
- **Structure:** Hierarchical groups and datasets
- **Metadata:** Attributes at any level
- **Compression:** gzip (level 6-9)
- **Chunking:** Enabled for efficient partial reads

---

## Tools and Software

### Python
```bash
pip install astropy photutils scipy matplotlib h5py
```

**Key Libraries:**
- `astropy`: FITS I/O, coordinates, units
- `photutils`: Source detection, photometry
- `specutils`: Spectroscopic analysis
- `lightkurve`: Exoplanet and variable star analysis

### Interactive Viewers
- **DS9:** FITS image viewer (http://ds.9.si.edu/)
- **fv:** FITS file viewer from HEASARC
- **Aladin:** Sky atlas and FITS viewer
- **Topcat:** Table viewer for catalogs

### Command Line
```bash
# FITS header
fitsheader galaxy_ngc1234_r_band.fits

# Convert FITS to PNG
fits2png galaxy_ngc1234_r_band.fits

# HDF5 inspection
h5dump -H variable_star_lightcurves.h5
```

---

## Scientific Applications

1. **Galaxy Image:**
   - Photometry and morphology analysis
   - PSF modeling and testing
   - Image processing pipeline validation
   - Star-galaxy separation algorithms

2. **Stellar Spectrum:**
   - Radial velocity measurements
   - Stellar parameter determination
   - Spectral line identification
   - Template matching algorithms

3. **Variable Stars:**
   - Period determination (Lomb-Scargle, Phase Dispersion Minimization)
   - Light curve classification
   - Distance measurements (Cepheids)
   - Time-series analysis methods

4. **Exoplanet Transits:**
   - Transit parameter fitting
   - Planet radius determination
   - Orbital mechanics validation
   - Systematic error modeling

---

## References

### Standards
- **FITS:** https://fits.gsfc.nasa.gov/
- **IAU FITS Working Group:** https://fits.gsfc.nasa.gov/iaufwg/

### Surveys and Archives
- **Sloan Digital Sky Survey:** https://www.sdss.org/
- **Hubble Legacy Archive:** https://hla.stsci.edu/
- **Gaia Archive:** https://gea.esac.esa.int/

### Software
- **Astropy:** https://www.astropy.org/
- **Photutils:** https://photutils.readthedocs.io/
- **Lightkurve:** https://docs.lightkurve.org/

---

## Citation

```
GRC Datasets Repository - Astronomy Collection
Created: October 2024
Formats: FITS, HDF5
Repository: https://github.com/grc-iit/Datasets
```

---

## License

These datasets are provided for educational, research, and benchmarking purposes.
