# Computational Fluid Dynamics Datasets

## Johns Hopkins Turbulence Database (JHTDB)
**Size:** 1+ Petabyte
**Format:** HDF5, Zarr (recent datasets)
**DOI:** 10.17616/R3FK7B

The premier open-access turbulence resource with multiple DNS and LES datasets accessible via web service APIs.

**Access:** [turbulence.pha.jhu.edu](http://turbulence.pha.jhu.edu/)

**Available Datasets:**
- **Isotropic Turbulence DNS:** 100 TB
- **Forced Turbulent Channel Flow:** 130 TB
- **Transitional Boundary Layer:** 105 TB
- **Wind Farm LES:** 40 TB
- **32,768³ Mega DNS:** ~0.5 PB

**Access Methods:**
- Web service APIs (Python, MATLAB, Fortran, C)
- Cutout service for HDF5/Zarr subsets
- No supercomputer required - analyze from your laptop

**Python Example:**
```python
import pyJHTDB

# Initialize
lJHTDB = pyJHTDB.libJHTDB()
lJHTDB.initialize()

# Get velocity at specific points
time = 0.002
points = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
result = lJHTDB.getData(time, points, sinterp=4, tinterp=0,
                        data_set='isotropic1024coarse',
                        getFunction='getVelocity')
```

---

## Automotive Aerodynamics

### DrivAerNet++
**Size:** 39 TB
**Format:** HDF5, VTK, Parquet
**License:** CC BY-NC 4.0

8,000 parametric car designs with complete CFD fields (3M CPU-hours on MIT Supercloud).

**Access:**
- GitHub: [github.com/Mohamedelrefaie/DrivAerNet](https://github.com/Mohamedelrefaie/DrivAerNet)
- Harvard Dataverse: Search "DrivAerNet"

**Data Includes (per design):**
- 3D meshes (STL)
- Surface pressure and wall shear stress
- Full volumetric CFD fields
- Point clouds
- Aerodynamic coefficients (Cd, Cl, Cs)
- Semantic segmentation (29 car parts)
- 2D renderings

**Car Types:** Fastback, notchback, estateback configurations

### WindsorML
**Size:** 8 TB
**Format:** VTU, STL, STEP, CSV
**License:** CC BY-SA 4.0

355 Windsor body geometric variants with high-fidelity wall-modeled LES.

**Access:**
- Website: [neilashton.github.io/caemldatasets/windsorml](https://neilashton.github.io/caemldatasets/windsorml)
- HuggingFace: Search "WindsorML"

**Computational Details:**
- ~300M cells per case
- Wall-modeled LES
- Comprehensive flow statistics
- Geometry variations for ML training

---

## Machine Learning-Focused CFD

### AirfRANS
**Size:** Moderate (optimized for ML)
**Format:** Parquet
**License:** ODbL

Reynolds-Averaged Navier-Stokes simulations over airfoils optimized for ML pipelines.

**Access:**
- Data: [data.isir.upmc.fr/extrality/NeurIPS_2022](https://data.isir.upmc.fr/extrality/NeurIPS_2022)
- Code: [github.com/Extrality/AirfRANS](https://github.com/Extrality/AirfRANS)
- Docs: [airfrans.readthedocs.io](https://airfrans.readthedocs.io/)

**Data Includes:**
- Pressure fields
- Velocity fields
- Various angles of attack
- Parquet format for fast columnar access

**Python Example:**
```python
import pandas as pd

# Load airfoil data
df = pd.read_parquet('airfrans_dataset.parquet')

# Access pressure and velocity fields
pressure = df['pressure'].values
velocity = df[['u', 'v']].values
```

### MegaFlow2D
**Size:** Multi-GB
**Format:** PyTorch Geometric

2 million+ snapshots from 3,000 different 2D fluid configurations with multi-fidelity data.

**Data Includes:**
- Various Reynolds numbers
- Different geometries
- Time-series data
- Super-resolution training pairs

**Framework:** FEniCS/Oasis with PyTorch Geometric interface

### Curated Turbulence Modelling Dataset
**Size:** Moderate
**Format:** NumPy arrays, OpenFOAM cases
**DOI:** 10.1038/s41597-021-01034-2

First open-source dataset specifically for ML-augmented turbulence closure models.

**Access:**
- Nature Data: [nature.com/articles/s41597-021-01034-2](https://www.nature.com/articles/s41597-021-01034-2)
- Kaggle: Search "Turbulence Modelling"

**Data Includes:**
- 895,640 data points per model
- Models: k-ε, k-ω, k-ω SST, k-ε-ϕt-f
- 47 tensor basis invariants
- OpenFOAM case files

---

## Download Instructions

### JHTDB Cutout Service
```python
import pyJHTDB

# Get a spatial cutout
lJHTDB = pyJHTDB.libJHTDB()
lJHTDB.initialize()

# Define region of interest
start = [0, 0, 0]
end = [128, 128, 128]
step = [1, 1, 1]

# Download as HDF5
result = lJHTDB.getCutout(start, end, step,
                          data_set='isotropic1024coarse',
                          field='velocity',
                          time=0.002,
                          output_format='hdf5',
                          filename='turbulence_cutout.h5')
```

### DrivAerNet Download
```bash
# Clone repository
git clone https://github.com/Mohamedelrefaie/DrivAerNet.git

# Download specific dataset components
# (Follow repository instructions for Dataverse access)
```

### AirfRANS
```python
from airfrans import dataset

# Load training data
train_data = dataset.load('train')

# Access specific samples
sample = train_data[0]
print(f"Pressure: {sample['pressure'].shape}")
print(f"Velocity: {sample['velocity'].shape}")
```

---

## Resources

- **JHTDB GitHub:** [github.com/idies/pyJHTDB](https://github.com/idies/pyJHTDB)
- **OpenFOAM:** [openfoam.org](https://openfoam.org/)
- **PyTorch Geometric:** [pytorch-geometric.readthedocs.io](https://pytorch-geometric.readthedocs.io/)
