# Oceanography Datasets

This directory contains oceanographic datasets in NetCDF format following CF (Climate and Forecast) conventions.

## Datasets

### 1. CTD Profiles - North Atlantic (2.4 MB total)

**File:** `ctd_profiles_atlantic_2024.nc` (37 KB)

A realistic CTD (Conductivity-Temperature-Depth) profile dataset from a simulated North Atlantic research cruise.

**Contents:**
- 12 CTD casts along a transect (40°N, 60°W to 20°W)
- 200 depth levels per profile (0-2000 meters)
- Variables: temperature, salinity, pressure
- Temporal spacing: 6 hours between casts
- Features realistic thermocline and halocline structures
- Mediterranean water mass signature in eastern stations

**Dimensions:**
- `profile`: 12 (time-based, one per CTD cast)
- `depth`: 200 levels

**Variables:**
- `temperature` (°C): Sea water temperature with thermocline at ~100m
- `salinity` (PSU): Practical salinity with halocline features
- `pressure` (dbar): Sea water pressure
- `latitude`, `longitude`: Station locations
- `time`: Cast timing in hours since 2024-07-01

**Usage Example (Python):**
```python
import netCDF4 as nc
import matplotlib.pyplot as plt

# Read data
ds = nc.Dataset('ctd_profiles_atlantic_2024.nc')
temp = ds.variables['temperature'][:]
depth = ds.variables['depth'][:]
lat = ds.variables['latitude'][:]

# Plot first profile
plt.plot(temp[0, :], depth, 'b-')
plt.gca().invert_yaxis()
plt.xlabel('Temperature (°C)')
plt.ylabel('Depth (m)')
plt.title(f'CTD Profile at {lat[0]:.2f}°N')
plt.show()
```

**Usage Example (NCO/CDO):**
```bash
# View metadata
ncdump -h ctd_profiles_atlantic_2024.nc

# Extract temperature only
ncks -v temperature ctd_profiles_atlantic_2024.nc temp_only.nc

# Get statistics
ncstat -s '%avg %sdev' -v temperature ctd_profiles_atlantic_2024.nc
```

---

### 2. Ocean Surface Analysis - Global (2.3 MB)

**File:** `ocean_surface_analysis_2024.nc` (2.3 MB)

Global gridded ocean surface dataset with monthly SST and chlorophyll-a concentration.

**Contents:**
- Monthly data for 2024 (12 time steps)
- Global 1° × 1° grid (180 × 360 spatial resolution)
- Sea surface temperature with realistic patterns
- Chlorophyll concentration with seasonal blooms
- Features: Equatorial upwelling, North Atlantic spring bloom, polar productivity

**Dimensions:**
- `time`: 12 (monthly, mid-month values)
- `lat`: 180 (1° resolution, -89.5° to 89.5°)
- `lon`: 360 (1° resolution, -179.5° to 179.5°)

**Variables:**
- `sst` (°C): Sea surface temperature
  - Range: -2°C to 30°C
  - Seasonal variation stronger at mid-latitudes
  - Cold tongue in equatorial Pacific
- `chlorophyll` (mg m⁻³): Mass concentration of chlorophyll
  - Elevated in high latitudes (60°N) and equatorial zones
  - Spring bloom in North Atlantic (March-May)
  - Typical range: 0.05 - 2.0 mg m⁻³

**Usage Example (Python):**
```python
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Read with xarray
ds = xr.open_dataset('ocean_surface_analysis_2024.nc')

# Plot SST for March
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ds.sst.isel(time=2).plot(ax=ax, cmap='RdYlBu_r',
                          vmin=-2, vmax=30,
                          cbar_kwargs={'label': 'SST (°C)'})
ax.coastlines()
ax.set_title('Sea Surface Temperature - March 2024')
plt.show()

# Compute seasonal mean SST
seasonal_sst = ds.sst.groupby('time.season').mean()
print(seasonal_sst)
```

**Usage Example (CDO):**
```bash
# Extract single month (March)
cdo seltimestep,3 ocean_surface_analysis_2024.nc march_2024.nc

# Compute annual mean
cdo timmean ocean_surface_analysis_2024.nc annual_mean.nc

# Extract specific region (North Atlantic)
cdo sellonlatbox,-80,-10,35,65 ocean_surface_analysis_2024.nc north_atlantic.nc

# Compute seasonal climatology
cdo seasmean ocean_surface_analysis_2024.nc seasonal.nc
```

---

## File Formats

All datasets follow **NetCDF-4** (HDF5-based) format with **CF-1.8 conventions**.

### NetCDF Tools

**Python:**
- `netCDF4`: Low-level interface
- `xarray`: High-level, pandas-like interface (recommended)
- `h5netcdf`: Pure Python HDF5-based reader

**Command Line:**
- `ncdump`: View metadata and data
- NCO (NetCDF Operators): `ncks`, `ncwa`, `ncea`, etc.
- CDO (Climate Data Operators): Advanced climate data processing

**Other Languages:**
- MATLAB: `ncread()`, `ncdisp()`
- R: `ncdf4` package
- Julia: `NCDatasets.jl`

---

## Data Characteristics

### CTD Profile Dataset
- **Scientific Domain:** Physical oceanography
- **Application:** Ocean mixing, water mass analysis, CTD sensor testing
- **Format Compliance:** CF-1.8, COARDS
- **Compression:** NetCDF-4 internal compression (level 4)
- **Quality:** Realistic T/S relationships following Atlantic water mass properties

### Ocean Surface Dataset
- **Scientific Domain:** Satellite oceanography, biogeochemistry
- **Application:** Climate analysis, ecosystem modeling, algorithm validation
- **Format Compliance:** CF-1.8
- **Compression:** NetCDF-4 internal compression (level 6)
- **Quality:** Physically realistic SST patterns and biological productivity

---

## References

### Standards
- **CF Conventions:** http://cfconventions.org/
- **NetCDF User Guide:** https://www.unidata.ucar.edu/software/netcdf/docs/

### Related Datasets
- NOAA World Ocean Atlas: https://www.ncei.noaa.gov/products/world-ocean-atlas
- Copernicus Marine Service: https://marine.copernicus.eu/
- NOAA OceanColor: https://oceancolor.gsfc.nasa.gov/

### Tools
- **Xarray:** https://xarray.dev/
- **NCO:** http://nco.sourceforge.net/
- **CDO:** https://code.mpimet.mpg.de/projects/cdo/
- **Panoply:** https://www.giss.nasa.gov/tools/panoply/ (NetCDF viewer)

---

## Citation

If using these datasets for testing or benchmarking:

```
GRC Datasets Repository - Oceanography Collection
Created: October 2024
Format: NetCDF-4 (CF-1.8 compliant)
Repository: https://github.com/grc-iit/Datasets
```

---

## License

These datasets are provided for educational, research, and benchmarking purposes.
