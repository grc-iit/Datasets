# Climate and Earth Science Datasets

## CMIP6 - Climate Model Intercomparison Project Phase 6
**Size:** Hundreds of terabytes
**Format:** NetCDF-4 (HDF5-based)

Global climate projections from 23+ general circulation models providing the scientific foundation for IPCC assessment reports.

**Access:**
- Primary: [esgf-node.llnl.gov/projects/cmip6](https://esgf-node.llnl.gov/projects/cmip6)
- NASA NEX-GDDP-CMIP6: [nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6](https://nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6)
- AWS: `registry.opendata.aws/nex-gddp-cmip6`

**Data Includes:**
- Daily and monthly gridded output
- Temperature, precipitation, wind, and 50+ climate variables
- Historical simulations (1850-2014)
- Future scenarios (SSP1-2.6 through SSP5-8.5, 2015-2100)
- Multiple ensemble members per model

**Tools:** Python (xarray, pandas), NCO, CDO

---

## ERA5 - ECMWF Reanalysis
**Size:** 390 TB (complete model-level dataset)
**Format:** NetCDF-4/HDF5 with CF conventions

Hourly atmospheric reanalysis from 1950-present at 31km resolution with 137 vertical levels.

**Access:**
- NCAR Research Data Archive: [rda.ucar.edu/datasets/d633006](https://rda.ucar.edu/datasets/d633006)
- AWS S3: `registry.opendata.aws/ecmwf-era5`
- Copernicus Climate Data Store: [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu)

**Data Includes:**
- Temperature, precipitation, winds at all levels
- Surface variables (2m temperature, 10m wind, etc.)
- Pressure levels and model levels
- Uncertainty estimates

**Tools:** ecmwf-opendata Python package, xarray, Copernicus CDS API

---

## NCAR CESM Large Ensemble
**Size:** 500 TB (70 TB subset on AWS)
**Format:** NetCDF and Zarr (cloud-optimized)

40-member ensemble of Community Earth System Model simulations (1920-2100) under RCP8.5 scenario.

**Access:**
- AWS S3: `registry.opendata.aws/ncar-cesm-lens`
- NCAR: [www.cesm.ucar.edu/community-projects/lens](https://www.cesm.ucar.edu/community-projects/lens)

**Data Includes:**
- 2D and 3D variables across atmosphere, ocean, land, ice
- Monthly, daily, and 6-hourly output
- Complete carbon cycle
- Chemistry and aerosols

---

## NASA Earth Science Data

### GES DISC - Goddard Earth Sciences Data Center
**Size:** Multi-petabyte
**Format:** HDF5, HDF-EOS5, NetCDF

**Access:** [disc.gsfc.nasa.gov](https://disc.gsfc.nasa.gov)

**Key Datasets:**
- **GPM (Global Precipitation Measurement):** Level 1-3 precipitation data in HDF5
- **AIRS:** Atmospheric Infrared Sounder in HDF-EOS5
- **MERRA-2:** Modern-Era Retrospective analysis for Research and Applications
- **OMI:** Ozone Monitoring Instrument atmospheric composition

### PO.DAAC - Physical Oceanography DAAC
**Size:** Multi-petabyte
**Format:** HDF5, NetCDF

**Access:** [podaac.jpl.nasa.gov](https://podaac.jpl.nasa.gov)

**Key Datasets:**
- **Aquarius:** Sea surface salinity (HDF5)
- **SWOT:** Surface Water and Ocean Topography
- **GRACE:** Gravity Recovery and Climate Experiment
- **Jason/Sentinel:** Sea surface height measurements

---

## WRF Regional Climate Model
**Size:** 2.14 TB
**Format:** NetCDF

High-resolution regional climate simulations for North America, Europe, and oceanic regions.

**Access:** [rda.ucar.edu/datasets/ds601.0](https://rda.ucar.edu/datasets/ds601.0)

**Resolution:** 3-hourly and 6-hourly output

---

## Download Instructions

Most datasets provide multiple access methods:

1. **Web portals** with subsetting capabilities
2. **OPeNDAP** for remote data access without full downloads
3. **AWS S3** for cloud-based analysis
4. **Globus** for high-performance bulk transfers
5. **Python APIs** (podaac-data-downloader, cdsapi, etc.)

### Example: ERA5 subset download
```python
import cdsapi

c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'variable': 'temperature',
        'pressure_level': '850',
        'year': '2023',
        'month': '01',
        'time': '00:00',
        'format': 'netcdf'
    },
    'download.nc')
```

### Example: CMIP6 wget download
```bash
wget -i filelist.txt  # Generated from ESGF search interface
```

---

## References

- **CF Conventions:** [cfconventions.org](http://cfconventions.org/)
- **ESGF:** [esgf.llnl.gov](https://esgf.llnl.gov/)
- **NASA Earthdata:** [earthdata.nasa.gov](https://earthdata.nasa.gov/)
