# HDF5 Datasets

HDF5 (Hierarchical Data Format version 5) is the dominant format for scientific computing, providing:
- Parallel I/O for HPC applications
- Self-describing hierarchical structure
- Efficient compression
- Random access to large datasets

## Available Datasets

*Currently building this collection. Small example datasets will be added here.*

## Shadow Datasets (Large Collections)

For HDF5 datasets too large to host on GitHub, see our shadow dataset documentation:

### Astrophysics & Cosmology
- **IllustrisTNG** (1.1 PB): [Shadow-Datasets/Astrophysics-Cosmology.md](../Shadow-Datasets/Astrophysics-Cosmology.md)
- **CAMELS** (200+ TB): Cosmology for machine learning
- **CosmoFlow** (100 GB): ML-ready N-body simulations

### Climate Science
- **ERA5** (390 TB as NetCDF-4/HDF5): [Shadow-Datasets/Climate-Earth-Science.md](../Shadow-Datasets/Climate-Earth-Science.md)
- **NASA Earth Science** (Multi-PB in HDF5/HDF-EOS5)
- **GPM Precipitation** (HDF5 format)

### Fluid Dynamics
- **Johns Hopkins Turbulence DB** (1+ PB): [Shadow-Datasets/Computational-Fluid-Dynamics.md](../Shadow-Datasets/Computational-Fluid-Dynamics.md)
- **DrivAerNet++** (39 TB): Automotive CFD with HDF5 fields

### Molecular Dynamics
- **mdCATH** (3+ TB): [Shadow-Datasets/Molecular-Dynamics.md](../Shadow-Datasets/Molecular-Dynamics.md)
- **ANI-1** (100+ GB): Quantum chemistry for neural network potentials

### Materials Science
- **NOMAD** (19M+ records): [Shadow-Datasets/Materials-Science.md](../Shadow-Datasets/Materials-Science.md)

## Working with HDF5

### Python (h5py)
```python
import h5py
import numpy as np

# Read HDF5 file
with h5py.File('dataset.h5', 'r') as f:
    # List all groups
    print("Keys:", list(f.keys()))

    # Access dataset
    data = f['group/dataset'][:]

    # Read attributes
    attrs = dict(f['group'].attrs)

# Write HDF5 file
with h5py.File('output.h5', 'w') as f:
    # Create group
    grp = f.create_group('simulation')

    # Create dataset with compression
    grp.create_dataset('data', data=np.random.rand(1000, 1000),
                       compression='gzip', compression_opts=9)

    # Add metadata
    grp.attrs['description'] = 'Simulation results'
    grp.attrs['timestamp'] = '2024-01-01'
```

### Parallel HDF5 (MPI)
```python
from mpi4py import MPI
import h5py

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Parallel HDF5 write
with h5py.File('parallel.h5', 'w', driver='mpio', comm=comm) as f:
    # Each rank writes to its portion
    dset = f.create_dataset('data', (size, 1000), dtype='f8')
    dset[rank, :] = np.random.rand(1000) * rank
```

### Command Line (h5dump, h5ls)
```bash
# List contents
h5ls -r file.h5

# Dump data as text
h5dump file.h5

# Extract specific dataset
h5dump -d /group/dataset file.h5

# Get file structure
h5stat file.h5
```

## Related Formats

- **NetCDF-4:** Built on HDF5, adds CF conventions (see `../NetCDF/`)
- **HDF-EOS5:** NASA's Earth Science extension of HDF5
- **H5MD:** Molecular Dynamics standardized HDF5 structure
- **MDTraj HDF5:** MD trajectory format

## Resources

- **HDF5 Official:** [hdfgroup.org](https://www.hdfgroup.org/)
- **h5py Documentation:** [docs.h5py.org](https://docs.h5py.org/)
- **HDF5 Tutorial:** [portal.hdfgroup.org/display/HDF5/Learning+HDF5](https://portal.hdfgroup.org/display/HDF5/Learning+HDF5)
- **Parallel HDF5:** [portal.hdfgroup.org/display/HDF5/Parallel+HDF5](https://portal.hdfgroup.org/display/HDF5/Parallel+HDF5)
