# Materials Science and Electronic Structure Datasets

## The Materials Project
**Size:** 154,000+ materials, 172,000+ molecules
**Format:** JSON, CIF, POSCAR
**License:** Free with registration

Most comprehensive computed materials database using VASP DFT calculations.

**Access:**
- Website: [materialsproject.org](https://materialsproject.org/)
- Next-gen: [next-gen.materialsproject.org](https://next-gen.materialsproject.org/)
- AWS: `registry.opendata.aws/materials-project`

**Data Includes:**
- Formation energies (corrected to experimental accuracy)
- Band structures
- Density of states
- Phase diagrams
- Elastic properties
- Dielectric properties
- DFT functionals: GGA(+U), R²SCAN

**Python API Example:**
```python
from mp_api.client import MPRester

with MPRester("YOUR_API_KEY") as mpr:
    # Search for materials
    docs = mpr.materials.summary.search(
        elements=["Li", "Fe", "O"],
        num_sites=(0, 10)
    )

    # Get structure
    structure = mpr.materials.get_structure_by_material_id("mp-1234")

    # Get band structure
    bs = mpr.materials.get_bandstructure_by_material_id("mp-1234")
```

---

## OQMD - Open Quantum Materials Database
**Size:** 1.3+ million DFT-calculated materials
**Format:** JSON (via API), SQL database
**License:** CC-BY 4.0

Northwestern University's comprehensive materials database.

**Access:** [oqmd.org](https://oqmd.org/)

**Data Includes:**
- ICSD structures + prototype decorations
- Thermodynamic properties
- Structural properties
- DFT level: PBE with GGA+U corrections
- Accuracy: 0.096 eV/atom MAE vs experiments

**Features:**
- Full database download available
- JSON API
- OPTIMADE support

**API Example:**
```python
import requests

# Query OQMD
url = "http://oqmd.org/oqmdapi/formationenergy"
params = {"elements": "Fe-O", "limit": 10}

response = requests.get(url, params=params)
data = response.json()

for entry in data['data']:
    print(f"Formula: {entry['name']}, E_form: {entry['delta_e']} eV/atom")
```

---

## AFLOW
**Size:** 4+ million compounds
**Format:** AFLOW format, POSCAR, CIF, JSON
**License:** Open access

Duke University's largest materials repository.

**Access:**
- Website: [aflowlib.org](http://aflowlib.org/)
- Main: [aflow.org](https://aflow.org/)

**Data Includes:**
- Complete phase diagrams (650+ binary systems)
- Electronic band structures
- Magnetic properties
- 1,783+ crystallographic prototypes
- OPTIMADE support

**AFLOW Prototype Encyclopedia:**
- Systematic categorization of all structures
- Prototype-based material discovery

**REST API Example:**
```bash
# Search for materials
curl "http://aflowlib.org/API/?species(Fe:O),Egap(2:*)"

# Download specific entry
curl "http://aflowlib.org/AFLOWDATA/ICSD_WEB/HEX/Fe1O1_ICSD_15840"
```

---

## NOMAD - Novel Materials Discovery
**Size:** 19+ million materials science records
**Format:** HDF5, JSON
**License:** Open source

Open infrastructure aggregating data from 40+ simulation codes.

**Access:**
- Website: [nomad-lab.eu](https://nomad-lab.eu/)
- Repository: [repository.nomad-coe.eu](https://repository.nomad-coe.eu/)

**Supported Codes:**
- VASP, FHI-aims, Gaussian, Quantum ESPRESSO
- CASTEP, CP2K, GPAW, WIEN2k
- Many more...

**Special Datasets:**

### HSE Hybrid Functional Database
**Size:** 7,024 materials
**Functional:** HSE06 (expensive, rare in databases)

Provides accurate electronic structure beyond standard GGA approximations.

**Features:**
- FAIR principles (Findable, Accessible, Interoperable, Reusable)
- Automatic metadata extraction
- Comprehensive provenance tracking

**Python Example:**
```python
from nomad.client import ArchiveQuery

# Search NOMAD
query = ArchiveQuery(
    required={
        'results.material.elements': ['Li', 'Co', 'O']
    }
)

results = query.execute()
for entry in results:
    print(entry.data)
```

---

## JARVIS - Joint Automated Repository for Various Integrated Simulations
**Size:** Large collection
**Format:** JSON
**License:** Open (NIST)

NIST's comprehensive materials database with DFT, force fields, and ML-ready datasets.

**Access:** [jarvis.nist.gov](https://jarvis.nist.gov/)

**Data Includes:**
- DFT calculations
- Force field parameters
- ML-ready datasets
- Benchmarking tools

**Python Tools:**
```python
from jarvis.db.figshare import data

# Load DFT 3D dataset
dft_3d = data('dft_3d')

# Access materials
for entry in dft_3d:
    print(f"Formula: {entry['formula']}")
    print(f"Formation energy: {entry['formation_energy_peratom']}")
```

---

## Protein Data Bank (PDB)
**Size:** 230,000+ structures
**Format:** PDB (legacy), PDBx/mmCIF (current), PDBML (XML)
**License:** CC0 (public domain)

Global repository of experimentally determined 3D biomolecular structures.

**Access:**
- RCSB: [rcsb.org](https://www.rcsb.org/)
- wwPDB: [wwpdb.org](https://www.wwpdb.org/)

**Data Includes:**
- X-ray crystallography structures
- NMR structures
- Cryo-EM structures
- Proteins, nucleic acids, complex assemblies

**Weekly Updates:** New structures added every Wednesday

**Python Example:**
```python
from Bio.PDB import PDBParser, PDBList

# Download structure
pdbl = PDBList()
pdbl.retrieve_pdb_file('1ABC', pdir='.', file_format='mmCif')

# Parse structure
parser = PDBParser()
structure = parser.get_structure('1ABC', '1abc.cif')

# Access atoms
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom.get_coord())
```

---

## AlphaFold Protein Structure Database
**Size:** 200+ million AI-predicted structures
**Format:** PDB, mmCIF
**License:** CC-BY 4.0

DeepMind's AI-predicted protein structures covering entire proteomes.

**Access:** [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/)

**Coverage:**
- Human proteome
- 48 key organisms
- SwissProt sequences
- Confidence scores included

---

## PubChem
**Size:** 111M+ compounds, 293M+ substances, 1.25M+ bioassays
**Format:** SDF, SMILES, InChI, XML, JSON

World's largest free chemical information resource.

**Access:**
- Website: [pubchem.ncbi.nlm.nih.gov](https://pubchem.ncbi.nlm.nih.gov/)
- FTP: [ftp.ncbi.nlm.nih.gov/pubchem](ftp://ftp.ncbi.nlm.nih.gov/pubchem)

**Data Includes:**
- Chemical structures
- Physical properties
- Biological activities
- Safety data
- Patents and literature
- PubChem3D spatial structures

**REST API Example:**
```bash
# Get compound properties
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/MolecularFormula,MolecularWeight/JSON"

# Download SDF
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF" > aspirin.sdf
```

---

## OPTIMADE - Open Databases Integration for Materials Design

**Standard:** [optimade.org](https://www.optimade.org/)

Unified API for querying across multiple materials databases.

**Participating Databases:**
- Materials Project
- OQMD
- AFLOW
- NOMAD
- COD (Crystallography Open Database)
- Many others

**Example Query:**
```python
import requests

# Query Materials Project via OPTIMADE
url = "https://materialsproject.org/optimade/v1/structures"
params = {
    "filter": 'elements HAS ALL "Si","O" AND nelements=2',
    "response_fields": "chemical_formula_descriptive,nsites"
}

response = requests.get(url, params=params)
data = response.json()
```

---

## Download Instructions

### Bulk Downloads

**Materials Project:**
```bash
# Via AWS S3
aws s3 sync s3://materialsproject-parsed/ ./mp_data/ --no-sign-request
```

**OQMD:**
```bash
# Download full database
wget http://oqmd.org/static/downloads/qmpy_ver3.db.gz
gunzip qmpy_ver3.db.gz
```

**PDB:**
```bash
# Mirror entire PDB
rsync -rlpt -v -z --delete \
    rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./pdb_mirror/
```

---

## Resources

- **pymatgen:** [pymatgen.org](https://pymatgen.org/) - Python materials analysis library
- **ASE:** [wiki.fysik.dtu.dk/ase](https://wiki.fysik.dtu.dk/ase/) - Atomic Simulation Environment
- **OPTIMADE:** [optimade.org](https://www.optimade.org/)
- **FAIR Materials Data:** [fairmat-nfdi.eu](https://fairmat-nfdi.eu/)
