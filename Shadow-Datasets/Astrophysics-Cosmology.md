# Astrophysics and Cosmology Datasets

## IllustrisTNG
**Size:** 1.1 Petabytes
**Format:** HDF5

One of the most comprehensive public cosmological simulation suites with three volumes (TNG50, TNG100, TNG300) at varying resolutions.

**Access:** [tng-project.org/data](https://www.tng-project.org/data)

**Data Includes:**
- 2,000 full volume snapshots
- ~110,000 high time-resolution subbox snapshots
- 100 redshifts (z=0 to z=20)
- Dark matter, gas, stars, and black holes
- Halo/subhalo catalogs and merger trees
- Python, IDL, and MATLAB helper scripts

**Web Interface:** JupyterLab environment for browser-based analysis

**API Example:**
```python
import illustris_python as il

# Get snapshot data
snap = il.snapshot.loadSnapshot(basePath, 99, 'gas')

# Get halo catalog
halos = il.groupcat.loadHalos(basePath, 99, fields=['GroupMass'])
```

---

## CAMELS - Cosmology and Astrophysics with MachineLearning Simulations
**Size:** 200+ TB
**Format:** HDF5

4,233 cosmological simulations systematically varying cosmological and astrophysical parameters for machine learning applications.

**Access:** [camel-simulations.org/data](https://www.camel-simulations.org/data)

**Data Includes:**
- 1,092 IllustrisTNG simulations
- 1,092 SIMBA simulations
- N-body counterparts
- Total: 144,840 snapshots with SUBFIND halo catalogs

### CAMELS Multifield Dataset
**Size:** 70+ TB
**Format:** NumPy arrays (optimized for ML)

**Access:** [camels-multifield-dataset.readthedocs.io](https://camels-multifield-dataset.readthedocs.io/)

**Data Includes:**
- 2D maps (512×512) from 12,000+ simulated universes
- 3D grids (512×512×512)
- Gas, dark matter, and stellar properties
- Multiple redshifts per simulation

**Use Case:** Direct input for convolutional neural networks

---

## EAGLE - Evolution and Assembly of GaLaxies and their Environments
**Size:** 24 large-scale runs (~25 TB total)
**Format:** HDF5 + SQL database

Cosmological hydrodynamics simulations with comprehensive galaxy formation physics.

**Access:** [icc.dur.ac.uk/Eagle/database.php](http://icc.dur.ac.uk/Eagle/database.php)

**Data Includes:**
- Volumes: 25-100 Mpc
- 1 million+ galaxies
- Particle data in HDF5
- Halo catalogs
- SQL database for galaxy properties

**Query Example:**
```sql
SELECT GalaxyID, MassType_Star, StarFormationRate
FROM RefL0100N1504_SubHalo
WHERE MassType_Star > 1e10 AND Redshift < 0.1
```

---

## Millennium Simulation
**Size:** ~25 TB
**Format:** SQL database + particle data

Landmark dark matter simulation with 10 billion particles through a 500 Mpc/h box.

**Access:** [wwwmpa.mpa-garching.mpg.de/galform/millennium](https://wwwmpa.mpa-garching.mpg.de/galform/millennium)

**Data Includes:**
- 9 million semi-analytic galaxies
- Merger trees
- Halo catalogs
- Virtual Observatory SQL interface

---

## CosmoFlow
**Size:** ~100 GB (processed for ML)
**Format:** HDF5

Collection of ~10,000 cosmological N-body simulations designed for deep learning applications.

**Access:** [portal.nersc.gov/project/m3363](https://portal.nersc.gov/project/m3363)

**Data Includes:**
- 512³ particle histograms
- 4 redshifts per simulation
- Varied cosmological parameters (Ωm, σ8)
- Labels for regression/classification tasks

**Generated With:** MUSIC + pyCOLA codes

---

## Additional Dark Matter Datasets

### Dark Matter Flow Analysis
**Size:** Various (GB-scale)
**Format:** HDF5, CSV

**Access:**
- Halo-based statistics: [zenodo.org/records/6757654](https://zenodo.org/records/6757654)
- Correlation-based statistics: [zenodo.org/records/6541231](https://zenodo.org/records/6541231)

---

## Particle Physics: CERN Open Data

### LHC Collision Data
**Size:** 5+ Petabytes
**Format:** ROOT, HDF5

Real collision data from Large Hadron Collider experiments (CMS, ATLAS, ALICE, LHCb).

**Access:** [opendata.cern.ch](http://opendata.cern.ch/)

**Data Includes:**
- Detector data and Monte Carlo simulations
- Analysis tools and virtual machines
- Specialized ML datasets for particle identification
- Event display tools (Phoenix for ATLAS)

**CMS:** Over 2 PB including 100% of research data from certain runs

**ATLAS:** [opendata.atlas.cern](http://opendata.atlas.cern/)

---

## Download Instructions

### IllustrisTNG API Access
```python
import requests

# Get list of available files
baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key": "your_api_key"}

r = requests.get(baseUrl + "TNG100-1/snapshots/99/", headers=headers)
files = r.json()['files']

# Download specific file
file_url = files['particles.0.hdf5']['url']
wget file_url
```

### CAMELS Bulk Download
```bash
# Install camels python package
pip install camels

# Download specific simulation suite
python -c "from camels import download_data; download_data('IllustrisTNG', 'CV_0')"
```

### EAGLE Database Query
```python
import eagle_io

# Connect to database
con = eagle_io.connect('username', 'password')

# Execute query
data = con.execute_query("SELECT * FROM RefL0100N1504_SubHalo WHERE MassType_Star > 1e10 LIMIT 1000")
```

---

## References

- **IllustrisTNG:** [tng-project.org](https://www.tng-project.org/)
- **CAMELS:** [camel-simulations.org](https://www.camel-simulations.org/)
- **EAGLE:** [eagle-dataset.org](http://eagle-dataset.org/)
- **Virgo Consortium:** [virgo.dur.ac.uk](http://virgo.dur.ac.uk/)
