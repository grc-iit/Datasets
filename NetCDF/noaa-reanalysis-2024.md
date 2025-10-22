# NOAA/NCEP Reanalysis Daily Averages - 2024

**Format:** NetCDF-3/4
**Size:** 6.7 MB
**Source:** NOAA Physical Sciences Laboratory
**Access:** [psl.noaa.gov/data/gridded/data.ncep.reanalysis.html](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)

## File

### noaa_air_temp_2024.nc
**Variable:** Air Temperature at Surface (sigma level 0.995)
**Coverage:** Global, 2.5° × 2.5° resolution
**Time Range:** Full year 2024 (daily averages)
**Grid:** 73 latitudes × 144 longitudes × 366 days

## About NCEP/NCAR Reanalysis

The NCEP/NCAR Reanalysis Project uses a state-of-the-art analysis/forecast system to perform data assimilation using past data from 1948 to present. This creates a consistent climate dataset that combines observations with numerical weather prediction models.

**Key Features:**
- **Global coverage:** 2.5° × 2.5° lat/lon grid
- **Time period:** 1948 - present (updated monthly)
- **Variables:** 80+ atmospheric and surface parameters
- **Levels:** Surface and 17 pressure levels
- **Temporal resolution:** 6-hourly and daily averages

## Data Structure

```python
import xarray as xr

ds = xr.open_dataset('noaa_air_temp_2024.nc')
print(ds)
```

Output:
```
Dimensions:  (lat: 73, lon: 144, time: 366)
Coordinates:
  * lat      (lat) float32 90.0 87.5 85.0 ... -87.5 -90.0
  * lon      (lon) float32 0.0 2.5 5.0 ... 355.0 357.5
  * time     (time) datetime64[ns] 2024-01-01 ... 2024-12-31
Data variables:
    air      (time, lat, lon) float32 ...
Attributes:
    Conventions:  CF-1.0
    title:        NCEP/NCAR Reanalysis
    institution:  National Centers for Environmental Prediction
```

## Usage Examples

### Python (xarray - Recommended)
```python
import xarray as xr
import matplotlib.pyplot as plt

# Load data
ds = xr.open_dataset('noaa_air_temp_2024.nc')
air_temp = ds['air']

# Get annual mean temperature
annual_mean = air_temp.mean(dim='time')

# Plot global temperature map
annual_mean.plot(figsize=(12, 6))
plt.title('2024 Annual Mean Surface Air Temperature (K)')
plt.show()

# Extract time series for a location (e.g., Chicago: 41.8°N, 87.6°W)
chicago = air_temp.sel(lat=40, lon=-90, method='nearest')
chicago.plot()
plt.title('Daily Surface Air Temperature - Chicago 2024')
plt.ylabel('Temperature (K)')
plt.show()

# Convert Kelvin to Celsius
air_temp_celsius = air_temp - 273.15

# Calculate monthly climatology
monthly = air_temp.groupby('time.month').mean('time')
monthly.plot(col='month', col_wrap=4, figsize=(15, 10))
```

### Python (netCDF4)
```python
from netCDF4 import Dataset
import numpy as np

# Open file
nc = Dataset('noaa_air_temp_2024.nc', 'r')

# Print file info
print("Variables:", nc.variables.keys())
print("Dimensions:", nc.dimensions.keys())

# Read data
air = nc.variables['air'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
time = nc.variables['time'][:]

# Get variable attributes
units = nc.variables['air'].units
long_name = nc.variables['air'].long_name

print(f"{long_name}: {units}")
print(f"Shape: {air.shape}")

# Calculate global mean for first day
global_mean = np.mean(air[0, :, :])
print(f"Global mean temperature (day 1): {global_mean:.2f} K")

nc.close()
```

### Command Line (NCO)
```bash
# Print file info
ncdump -h noaa_air_temp_2024.nc

# Calculate annual mean
ncwa -a time noaa_air_temp_2024.nc annual_mean.nc

# Extract summer months (Jun-Aug)
ncks -d time,151,243 noaa_air_temp_2024.nc summer_2024.nc

# Extract Northern Hemisphere
ncks -d lat,0.,90. noaa_air_temp_2024.nc northern_hemisphere.nc

# Convert to CSV for specific location
ncks -d lat,40. -d lon,280. -v air noaa_air_temp_2024.nc chicago.nc
```

### Command Line (CDO)
```bash
# Calculate time mean
cdo timmean noaa_air_temp_2024.nc annual_mean.nc

# Calculate seasonal means
cdo seasmean noaa_air_temp_2024.nc seasonal_mean.nc

# Calculate anomalies (need climatology file)
cdo sub noaa_air_temp_2024.nc climatology.nc anomalies.nc

# Regrid to 1° resolution
cdo remapbil,global_1 noaa_air_temp_2024.nc air_temp_1deg.nc
```

## Variable Information

**Air Temperature (air):**
- **Units:** Kelvin (K)
- **Sigma level:** 0.995 (~near-surface, ~50m above ground)
- **Valid range:** ~180-320 K (-93°C to 47°C)
- **Missing value:** -9.96921e+36

**Coordinates:**
- **Latitude:** -90° to 90° (2.5° spacing)
- **Longitude:** 0° to 357.5° (2.5° spacing)
- **Time:** Days since 1800-01-01 00:00:00

## Related Datasets

From the same source, you can download:
- **Precipitation rate** (prate.sfc.gauss.YYYY.nc)
- **Sea level pressure** (slp.YYYY.nc)
- **Geopotential height** (hgt.YYYY.nc)
- **Winds** (uwnd.YYYY.nc, vwnd.YYYY.nc)
- **Relative humidity** (rhum.sig995.YYYY.nc)

## Analysis Ideas

1. **Climate trends:** Compare with previous years
2. **Extreme events:** Identify heat waves and cold snaps
3. **Seasonal patterns:** Analyze temperature seasonality
4. **Regional climate:** Extract and analyze specific regions
5. **Validation:** Compare with observations or other reanalyses

## References

- **NCEP/NCAR Reanalysis:** [psl.noaa.gov/data/gridded/data.ncep.reanalysis.html](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)
- **Documentation:** [psl.noaa.gov/data/gridded/data.ncep.reanalysis.html](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)
- **Kalnay et al. (1996):** The NCEP/NCAR 40-Year Reanalysis Project. Bull. Amer. Meteor. Soc., 77, 437-471
