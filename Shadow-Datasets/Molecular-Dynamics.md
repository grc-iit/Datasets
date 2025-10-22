# Molecular Dynamics and Protein Simulation Datasets

## mdCATH
**Size:** 3+ TB
**Format:** HDF5 (primary), XTC, PDB

The largest all-atom protein domain MD dataset covering CATH database domains.

**Access:**
- Website: [open.playmolecule.org/mdcath](https://open.playmolecule.org/mdcath)
- HuggingFace: Search "mdCATH"

**Data Includes:**
- 5,398 CATH protein domains
- 5 temperatures (320-450K) × 5 replicates each
- 134,950 total trajectories
- **Unique: Instantaneous forces** alongside coordinates
- Force field: CHARMM22*
- Sampling: Every 1 ns
- Total simulation time: 62+ milliseconds accumulated

**Metrics Provided:**
- Secondary structure (DSSP)
- RMSD (Root Mean Square Deviation)
- RMSF (Root Mean Square Fluctuation)
- Radius of gyration

**HDF5 Structure Example:**
```python
import h5py

with h5py.File('mdcath_domain.h5', 'r') as f:
    coords = f['coordinates'][:]  # Shape: (timesteps, atoms, 3)
    forces = f['forces'][:]       # Shape: (timesteps, atoms, 3)
    time = f['time'][:]
    topology = f['topology'][:]
```

---

## GPCRmd
**Size:** 556.5 microseconds across 371 systems
**Format:** Standard MD formats + web database

Specialized community database for G-protein coupled receptor simulations.

**Access:** [gpcrmd.org](https://gpcrmd.org/)

**Data Includes:**
- 60% of available GPCR structures
- 181 ligand-GPCR complexes
- 190 apo systems
- 3×500 ns per system (standardized protocols)
- 23 international institutions contributing
- Web-based analysis tools
- Interactive visualization

---

## MD File Format Standards

### MDTraj HDF5 Format
**Documentation:** [mdtraj.org](https://mdtraj.org/)

Python library and HDF5 specification widely adopted in the MD community.

**Features:**
- Topology embedding in HDF5
- Flexible compression
- Fast random access
- Integration with major MD packages

**Example:**
```python
import mdtraj as md

# Load trajectory
traj = md.load('trajectory.h5')

# Access data
positions = traj.xyz  # (frames, atoms, 3)
topology = traj.topology

# Compute properties
rg = md.compute_rg(traj)
rmsd = md.rmsd(traj, traj, 0)

# Save subset
traj[::10].save('downsampled.h5')  # Every 10th frame
```

### H5MD - HDF5 for Molecular Dynamics
**Specification:** [nongnu.org/h5md](https://nongnu.org/h5md/)

Standardized HDF5 structure for MD simulations.

**Supported By:**
- LAMMPS
- HOOMD-blue
- ESPResSo
- GROMACS (via plugins)

**H5MD Structure:**
```
/particles/
    /all/
        /position/
            /step  # Time step indices
            /time  # Physical time
            /value # Coordinates (timestep, atom, xyz)
        /velocity/
            /step
            /time
            /value
        /force/
            /step
            /time
            /value
/observables/
    /temperature/
    /pressure/
    /energy/
```

---

## Quantum Chemistry Datasets

### QM9
**Size:** ~2 GB
**Format:** XYZ, SDF, database files

Benchmark dataset for ML molecular property prediction.

**Access:** [quantum-machine.org/datasets](http://quantum-machine.org/datasets/)

**Data Includes:**
- 133,885 equilibrium organic molecules
- Up to 9 heavy atoms (C, N, O, F)
- 13+ quantum properties per molecule
- DFT level: B3LYP/6-31G(2df,p)

**Properties:**
- HOMO-LUMO gaps
- Energies (internal, atomization, etc.)
- Dipole moments
- Polarizabilities
- Harmonic frequencies

### ANI-1
**Size:** 100+ GB
**Format:** HDF5
**DOI:** 10.6084/m9.figshare.c.3846712

20 million off-equilibrium conformations for training neural network potentials.

**Access:**
- Figshare: [doi.org/10.6084/m9.figshare.c.3846712](https://doi.org/10.6084/m9.figshare.c.3846712)
- GitHub: [github.com/isayev/ANI1_dataset](https://github.com/isayev/ANI1_dataset)

**Data Includes:**
- 57,462 small organic molecules
- ~20M conformations (100× more than QM9)
- DFT level: ωB97x/6-31G(d)
- Off-equilibrium structures for transferability

**HDF5 Example:**
```python
import h5py

with h5py.File('ani1_dataset.h5', 'r') as f:
    for molecule in f.keys():
        coordinates = f[molecule]['coordinates'][:]
        energies = f[molecule]['energies'][:]
        species = f[molecule]['species'][:]
```

### MultiXC-QM9
**Size:** Large
**Format:** Database files

QM9 molecules computed with 76 DFT functionals and 3 basis sets.

**Access:** [data.dtu.dk/collections/MultiXC-QM9/6185986](https://data.dtu.dk/collections/MultiXC-QM9/6185986)

**Data Includes:**
- 228 energy values per molecule
- Reaction energies for all monomolecular interconversions
- Enables delta learning across methods

### Hessian QM9
**Size:** Moderate
**Format:** HDF5
**DOI:** 10.1038/s41597-024-04361-2

First database with full molecular Hessian matrices.

**Data Includes:**
- 41,645 molecules
- Full numerical Hessian matrices
- Vibrational frequencies
- Normal modes
- DFT level: ωB97x/6-31G*
- Environments: vacuum + 3 solvents

---

## LAMMPS Datasets

See also: `Adios/Lammps/` in this repository for ADIOS2-format LAMMPS data.

**Common LAMMPS Output Formats:**
- **LAMMPS dump:** Text-based trajectory files
- **DCD:** Binary trajectory format
- **NetCDF:** For efficient I/O
- **H5MD:** HDF5-based molecular dynamics format
- **ADIOS2:** For exascale I/O performance

---

## Download Instructions

### mdCATH
```python
# Via HuggingFace datasets
from datasets import load_dataset

dataset = load_dataset("compsciencelab/mdcath")

# Access specific domain
domain_data = dataset['train'][0]
coords = domain_data['coordinates']
forces = domain_data['forces']
```

### GPCRmd
```bash
# Access via web interface with filtering
# Download specific trajectories after selection
# API access available for bulk downloads
```

### QM9
```bash
wget https://ndownloader.figshare.com/files/3195389 -O qm9.tar.bz2
tar -xjf qm9.tar.bz2
```

### ANI-1
```python
# Download from Figshare (large file)
import requests

url = "https://figshare.com/ndownloader/files/9057631"
r = requests.get(url, stream=True)

with open('ani1_dataset.h5', 'wb') as f:
    for chunk in r.iter_content(chunk_size=8192):
        f.write(chunk)
```

---

## Resources

- **MDTraj:** [mdtraj.org](https://mdtraj.org/)
- **MDAnalysis:** [mdanalysis.org](https://www.mdanalysis.org/)
- **OpenMM:** [openmm.org](http://openmm.org/)
- **LAMMPS:** [lammps.org](https://www.lammps.org/)
- **GROMACS:** [gromacs.org](http://www.gromacs.org/)
