# NetCDF Files - Climate and Atmospheric Data

NetCDF (Network Common Data Form) is the standard format for climate and atmospheric sciences. NetCDF-4 is built on HDF5 with CF (Climate and Forecast) metadata conventions.

**Documentation:** [xarray.pydata.org](https://xarray.pydata.org/) | **Format Spec:** [unidata.ucar.edu/software/netcdf](https://www.unidata.ucar.edu/software/netcdf/)

## Available Datasets

### noaa_air_temp_2024.nc
- **Size:** 6.7 MB
- **Type:** NCEP/NCAR Reanalysis data
- **Coverage:** Global, 2024 (366 days)
- **Resolution:** 73 lats × 144 lons
- **Source:** NOAA Physical Sciences Laboratory
- **Download:** [psl.noaa.gov/data/gridded/data.ncep.reanalysis.html](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)

### noaa_sea_level_pressure_2023.nc
- **Size:** 6.3 MB
- **Type:** Sea level pressure reanalysis
- **Coverage:** Global, 2023
- **Source:** NOAA/NCEP
- **Download:** [psl.noaa.gov](https://psl.noaa.gov/)

## Where to Find More Data

- **NOAA PSL:** [psl.noaa.gov/data/gridded](https://psl.noaa.gov/data/gridded/) - Reanalysis and gridded climate data
- **CMIP6:** [esgf-node.llnl.gov](https://esgf-node.llnl.gov/) - Climate model intercomparison (100s TB)
- **ERA5:** [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu/) - Global reanalysis (390 TB)
- **NCAR Climate Data Gateway:** [rda.ucar.edu](https://rda.ucar.edu/) - Weather and climate datasets
- **NASA Earthdata:** [earthdata.nasa.gov](https://earthdata.nasa.gov/) - Earth science data in NetCDF/HDF
- **THREDDS Catalogs:** Many institutions host NetCDF via THREDDS/OPeNDAP

## Shadow Datasets

For TB-scale NetCDF climate datasets, see:
- [Shadow-Datasets/Climate-Earth-Science.md](../Shadow-Datasets/Climate-Earth-Science.md)
