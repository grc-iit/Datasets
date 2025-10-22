# NetCDF Datasets

NetCDF (Network Common Data Form) is the standard format for climate and atmospheric sciences. NetCDF-4 is built on HDF5 with Climate and Forecast (CF) metadata conventions.

## Available Datasets

*Currently building this collection. Small example datasets will be added here.*

## Shadow Datasets (Large Collections)

For NetCDF datasets too large to host on GitHub, see our shadow dataset documentation:

### Climate & Weather
- **CMIP6** (100s of TB): [Shadow-Datasets/Climate-Earth-Science.md](../Shadow-Datasets/Climate-Earth-Science.md)
- **ERA5** (390 TB): Global reanalysis 1950-present
- **NCAR CESM** (500 TB): Climate model ensemble
- **WRF Regional** (2.14 TB): High-resolution weather simulations

### Ocean & Hydrology
- **HYCOM**: Global ocean model outputs
- **SWOT**: Surface water observations
- **Ocean reanalysis**: Various products in NetCDF

## NetCDF Format Variants

- **NetCDF-3 (Classic):** Simple, limited to 2GB files
- **NetCDF-4:** Built on HDF5, supports compression and parallel I/O
- **NetCDF-4 Classic:** NetCDF-4 with NetCDF-3 data model

## Working with NetCDF

### Python (xarray - Recommended)
```python
import xarray as xr
import matplotlib.pyplot as plt

# Read NetCDF file
ds = xr.open_dataset('climate_data.nc')

# Inspect structure
print(ds)
print(ds.data_vars)

# Select data
temperature = ds['temperature']
subset = temperature.sel(time='2023-01-01', lat=slice(30, 50))

# Compute statistics
mean_temp = temperature.mean(dim=['lat', 'lon'])

# Plot
temperature.isel(time=0).plot()
plt.show()

# Write NetCDF
ds.to_netcdf('output.nc', format='NETCDF4', compression={'zlib': True, 'complevel': 5})
```

### Python (netCDF4 library)
```python
from netCDF4 import Dataset
import numpy as np

# Read
with Dataset('data.nc', 'r') as nc:
    # Get variable
    temp = nc.variables['temperature'][:]

    # Get attributes
    units = nc.variables['temperature'].units

    # Get dimensions
    time = nc.variables['time'][:]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]

# Write
with Dataset('output.nc', 'w', format='NETCDF4') as nc:
    # Create dimensions
    nc.createDimension('time', None)  # unlimited
    nc.createDimension('lat', 180)
    nc.createDimension('lon', 360)

    # Create variables
    times = nc.createVariable('time', 'f8', ('time',))
    temps = nc.createVariable('temperature', 'f4', ('time', 'lat', 'lon'),
                              zlib=True, complevel=5)

    # Add attributes (CF conventions)
    temps.units = 'K'
    temps.long_name = 'Air Temperature'
    temps.standard_name = 'air_temperature'

    # Write data
    temps[0, :, :] = np.random.rand(180, 360)
```

### Command Line Tools (NCO - NetCDF Operators)
```bash
# Get info about file
ncdump -h file.nc

# Extract variable
ncks -v temperature file.nc output.nc

# Average over time dimension
ncwa -a time file.nc average.nc

# Concatenate files
ncrcat file1.nc file2.nc file3.nc combined.nc

# Subset region
ncks -d lat,30.,50. -d lon,-120.,-80. file.nc subset.nc

# Convert to different format
ncks --fl_fmt=netcdf4_classic file.nc output.nc
```

### CDO (Climate Data Operators)
```bash
# Calculate climatology
cdo ymonmean input.nc climatology.nc

# Select years
cdo selyear,2020/2023 input.nc subset.nc

# Compute anomalies
cdo sub input.nc climatology.nc anomalies.nc

# Regrid to different resolution
cdo remapbil,target_grid.txt input.nc regridded.nc

# Merge multiple files
cdo mergetime file*.nc merged.nc
```

## CF Conventions

NetCDF files should follow CF (Climate and Forecast) metadata conventions:

```python
import xarray as xr
import numpy as np

# Create CF-compliant dataset
ds = xr.Dataset(
    {
        'temperature': (
            ['time', 'lat', 'lon'],
            np.random.rand(10, 180, 360),
            {
                'long_name': 'Air Temperature',
                'standard_name': 'air_temperature',
                'units': 'K',
                'coordinates': 'time lat lon'
            }
        )
    },
    coords={
        'time': (
            ['time'],
            np.arange(10),
            {
                'long_name': 'Time',
                'standard_name': 'time',
                'units': 'days since 2023-01-01',
                'calendar': 'gregorian'
            }
        ),
        'lat': (
            ['lat'],
            np.linspace(-89.5, 89.5, 180),
            {
                'long_name': 'Latitude',
                'standard_name': 'latitude',
                'units': 'degrees_north'
            }
        ),
        'lon': (
            ['lon'],
            np.linspace(-179.5, 179.5, 360),
            {
                'long_name': 'Longitude',
                'standard_name': 'longitude',
                'units': 'degrees_east'
            }
        )
    },
    attrs={
        'Conventions': 'CF-1.8',
        'title': 'Example Climate Dataset',
        'institution': 'GRC IIT',
        'source': 'Model simulation',
        'history': '2024-01-01: Created'
    }
)

ds.to_netcdf('cf_compliant.nc')
```

## OPeNDAP Access

Many NetCDF datasets support remote access via OPeNDAP:

```python
import xarray as xr

# Access remote dataset without downloading
url = 'https://data.example.org/thredds/dodsC/dataset.nc'
ds = xr.open_dataset(url)

# Subset is downloaded on-demand
subset = ds['temperature'].sel(time='2023-01-01')
```

## Zarr Alternative

Zarr is emerging as a cloud-optimized alternative to NetCDF:

```python
import xarray as xr

# Read NetCDF
ds = xr.open_dataset('data.nc')

# Convert to Zarr
ds.to_zarr('data.zarr')

# Read Zarr
ds_zarr = xr.open_zarr('data.zarr')
```

## Resources

- **NetCDF Official:** [unidata.ucar.edu/software/netcdf](https://www.unidata.ucar.edu/software/netcdf/)
- **CF Conventions:** [cfconventions.org](http://cfconventions.org/)
- **xarray:** [xarray.pydata.org](https://xarray.pydata.org/)
- **NCO:** [nco.sourceforge.net](http://nco.sourceforge.net/)
- **CDO:** [code.mpimet.mpg.de/projects/cdo](https://code.mpimet.mpg.de/projects/cdo/)
