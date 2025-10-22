# GRC Datasets Repository

Centralized repository for datasets used by the [GRC organization at IIT](https://grc.iit.edu/). Contains scientific simulation datasets and documentation for accessing petabyte-scale public datasets.

## Dataset Statistics

| Category | Datasets | Formats | Total Size | Description |
|----------|----------|---------|------------|-------------|
| **ADIOS** | 23 | BP5 | 755 MB | CFD, MD, Weather simulations |
| **Oceanography** | 2 | NetCDF | 2.4 MB | CTD profiles, surface analysis |
| **Genomics** | 7 | FASTA, HDF5, SAM, VCF, FASTQ | 8.5 MB | Genomes, variants, RNA-seq |
| **Astronomy** | 4 | FITS, HDF5 | 2.7 MB | Images, spectra, light curves |
| **Seismology** | 3 | HDF5 | 16 MB | Earthquake data, noise, RFs |
| **Parquet** | 2 | Parquet | 48 MB | NYC taxi, analytics samples |
| **NetCDF** | 2 | NetCDF | 6.7 MB | NOAA climate data |
| **HDF5** | 5 | HDF5, PDB | 766 KB | OpenPMD, protein structures |
| **ROOT** | 2 | ROOT | 20.5 MB | Higgs analysis, tutorials |
| **FITS** | 2 | FITS | 4.8 MB | Hubble observations |
| **CIF** | 3 | CIF | 85 KB | Crystal structures |
| **Darshan** | Examples | LOG | 36 MB | I/O characterization traces |
| **Shadow** | 50+ | Various | PB-scale | Documentation for public data |

**Total local datasets:** ~840 MB across 50+ datasets
**Total accessible (shadow):** Petabytes of public scientific data

## Index

### Data Formats
- [Adios](Adios/) - ADIOS2 I/O framework datasets
- [HDF5](HDF5/) - Hierarchical Data Format files
- [NetCDF](NetCDF/) - Network Common Data Form files
- [Parquet](Parquet/) - Columnar format files
- [ROOT](ROOT/) - Particle physics data from CERN
- [FITS](FITS/) - Astronomy image and data files
- [CIF](Crystallography/) - Crystallographic Information Files (crystal structures)

### Scientific Domains
- [Oceanography](Oceanography/) - Ocean and marine data (NetCDF)
- [Astronomy](Astronomy/) - Astronomical observations (FITS, HDF5)
- [Seismology](Seismology/) - Earthquake and seismic data (HDF5)
- [Genomics](Genomics/) - Genomics and bioinformatics data (FASTA, HDF5, SAM, VCF, FASTQ)

### Tracking
- [Darshan-Traces](Darshan-Traces/) - HPC I/O characterization traces

### Shadow Datasets
Documentation for petabyte-scale public datasets:
- [Climate & Earth Science](Shadow-Datasets/Climate-Earth-Science.md)
- [Astrophysics & Cosmology](Shadow-Datasets/Astrophysics-Cosmology.md)
- [Computational Fluid Dynamics](Shadow-Datasets/Computational-Fluid-Dynamics.md)
- [Molecular Dynamics](Shadow-Datasets/Molecular-Dynamics.md)
- [Materials Science](Shadow-Datasets/Materials-Science.md)
