# Shadow Datasets

This directory contains documentation for large-scale scientific datasets that exceed GitHub's file size limits (>50MB). These "shadow datasets" provide comprehensive information, access instructions, and code examples for downloading and using petabyte-scale public datasets from authoritative sources.

## What are Shadow Datasets?

Shadow datasets are references to publicly accessible datasets that are too large to host directly in this repository. Each document provides:

- **Complete access information** with URLs and credentials
- **Download instructions** with code examples
- **Data format specifications** and file structures
- **Python/command-line examples** for working with the data
- **License information** and usage terms

## Available Categories

### [Climate and Earth Science](Climate-Earth-Science.md)
**Total Size:** Hundreds of terabytes to petabytes

- **CMIP6:** Global climate projections (100s of TB, NetCDF-4/HDF5)
- **ERA5:** Hourly atmospheric reanalysis (390 TB, NetCDF-4)
- **NCAR CESM:** Climate model ensemble (500 TB, NetCDF/Zarr)
- **NASA Earth Science:** Multi-PB archives (HDF5, HDF-EOS5)
- **WRF Regional Climate:** High-resolution simulations (2.14 TB, NetCDF)

### [Astrophysics and Cosmology](Astrophysics-Cosmology.md)
**Total Size:** Petabytes

- **IllustrisTNG:** Cosmological simulations (1.1 PB, HDF5)
- **CAMELS:** ML-focused cosmology (200+ TB, HDF5/NumPy)
- **EAGLE:** Galaxy formation simulations (25 TB, HDF5+SQL)
- **CERN Open Data:** LHC collision data (5+ PB, ROOT/HDF5)
- **CosmoFlow:** ML-ready N-body sims (~100 GB, HDF5)

### [Computational Fluid Dynamics](Computational-Fluid-Dynamics.md)
**Total Size:** Petabytes

- **Johns Hopkins Turbulence DB:** 1+ PB DNS/LES data (HDF5/Zarr)
- **DrivAerNet++:** Automotive aerodynamics (39 TB, HDF5/VTK/Parquet)
- **WindsorML:** High-fidelity LES (8 TB, VTU/STL)
- **AirfRANS:** ML-focused airfoil data (Moderate, Parquet)
- **MegaFlow2D:** 2M fluid snapshots (Multi-GB, PyTorch)

### [Molecular Dynamics](Molecular-Dynamics.md)
**Total Size:** Terabytes

- **mdCATH:** All-atom protein simulations (3+ TB, HDF5)
- **GPCRmd:** GPCR receptor database (556.5 μs simulation time)
- **QM9:** Quantum chemistry benchmark (2 GB, XYZ/SDF)
- **ANI-1:** Neural network potentials (100+ GB, HDF5)
- **MD Format Standards:** MDTraj HDF5, H5MD specifications

### [Materials Science](Materials-Science.md)
**Total Size:** Varies (databases with millions of entries)

- **Materials Project:** 154K+ materials (JSON/CIF, AWS S3)
- **OQMD:** 1.3M+ DFT calculations (JSON/SQL, 100s of GB)
- **AFLOW:** 4M+ compounds with phase diagrams (Multi-format)
- **NOMAD:** 19M+ records from 40+ codes (HDF5/JSON)
- **Protein Data Bank:** 230K+ structures (PDB/mmCIF)
- **AlphaFold:** 200M+ AI predictions (PDB/mmCIF)

## Why Shadow Datasets?

**GitHub's 50MB file limit** prevents direct hosting of most scientific datasets, which typically range from gigabytes to petabytes. Shadow datasets provide:

1. **Centralized Documentation:** One-stop reference for accessing major public datasets
2. **Reproducibility:** Exact download instructions and API examples
3. **Discovery:** Learn what data exists before investing in downloads
4. **Integration:** Code examples show how to use data in your workflows

## Format Overview

Most shadow datasets use these formats:

- **HDF5/NetCDF-4:** Dominant in climate, astrophysics, molecular dynamics
- **Parquet:** Growing adoption in ML-focused CFD and analytics
- **JSON/CIF/POSCAR:** Common in materials science databases
- **ROOT:** Particle physics standard (CERN)
- **Zarr:** Emerging cloud-optimized format

## Common Access Patterns

### Cloud-Hosted Data (No Download)
```python
import s3fs
import xarray as xr

# Access ERA5 on AWS S3
s3 = s3fs.S3FileSystem(anon=True)
with s3.open('s3://era5-pds/2023/01/data/air_temperature_at_2_metres.nc') as f:
    ds = xr.open_dataset(f)
```

### API Access
```python
from mp_api.client import MPRester

# Materials Project API
with MPRester("YOUR_API_KEY") as mpr:
    docs = mpr.materials.summary.search(elements=["Li", "Fe", "O"])
```

### Bulk Download
```bash
# Globus transfer for large datasets
globus transfer SOURCE_ENDPOINT_ID DEST_ENDPOINT_ID --recursive

# Or wget with parallel downloads
cat filelist.txt | xargs -P 8 -n 1 wget
```

## Contributing

Found a valuable public dataset? Add it as a shadow dataset:

1. Create a new markdown file in this directory
2. Follow the template format (see existing files)
3. Include: size, format, access methods, code examples
4. Submit a PR

## Resources

- **FAIR Data Principles:** [go-fair.org](https://www.go-fair.org/)
- **HDF5:** [hdfgroup.org](https://www.hdfgroup.org/)
- **NetCDF:** [unidata.ucar.edu/software/netcdf](https://www.unidata.ucar.edu/software/netcdf/)
- **Parquet:** [parquet.apache.org](https://parquet.apache.org/)
- **Zarr:** [zarr.readthedocs.io](https://zarr.readthedocs.io/)
