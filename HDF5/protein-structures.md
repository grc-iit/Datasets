# Protein Structure Data (PDB Format)

**Format:** PDB (Protein Data Bank)
**Size:** 49-131 KB
**Source:** RCSB Protein Data Bank [rcsb.org](https://www.rcsb.org/)

## Files

### protein_1CRN.pdb - Crambin
**PDB ID:** 1CRN
**Size:** 49 KB
**Resolution:** 0.945 Å (sub-angstrom!)
**Organism:** *Crambe abyssinica* (Abyssinian crambe)

Crambin is a small plant seed storage protein with 46 amino acid residues. It's one of the first proteins solved at sub-angstrom resolution, making it a valuable reference structure for crystallography methods.

**Highlights:**
- 327 atoms total
- 6 sulfur atoms forming 3 disulfide bonds
- Very high-resolution structure (0.945 Å)
- Classic test case for crystallographic refinement

### lysozyme_2LYZ.pdb - Lysozyme
**PDB ID:** 2LYZ
**Size:** 131 KB
**Resolution:** 1.33 Å
**Organism:** *Gallus gallus* (chicken)

Chicken egg-white lysozyme is one of the most studied enzymes. It catalyzes the hydrolysis of bacterial cell walls and was one of the first enzyme structures determined by X-ray crystallography (1965, Nobel Prize).

**Highlights:**
- 1,001 atoms
- 129 amino acid residues
- Active site with catalytic residues Glu35 and Asp52
- Classic model system for protein folding and dynamics

## Data Structure (PDB Format)

```
HEADER    ...
TITLE     ...
COMPND    MOL_ID: 1; MOLECULE: ...
SOURCE    ORGANISM_SCIENTIFIC: ...
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
...
```

**Columns:**
- Record type (ATOM, HETATM, etc.)
- Atom serial number
- Atom name
- Residue name
- Chain ID
- Residue sequence number
- X, Y, Z coordinates (Å)
- Occupancy
- Temperature factor (B-factor)
- Element symbol

## Usage

### Python (Biopython)
```python
from Bio.PDB import PDBParser

parser = PDBParser()
structure = parser.get_structure('crambin', 'protein_1CRN.pdb')

# Get all atoms
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(f"{atom.name}: {atom.coord}")

# Calculate center of mass
atoms = list(structure.get_atoms())
coords = [atom.coord for atom in atoms]
import numpy as np
center = np.mean(coords, axis=0)
print(f"Center of mass: {center}")
```

### Python (MDAnalysis)
```python
import MDAnalysis as mda

u = mda.Universe('protein_1CRN.pdb')

print(f"Number of atoms: {len(u.atoms)}")
print(f"Number of residues: {len(u.residues)}")

# Select protein atoms
protein = u.select_atoms('protein')

# Calculate radius of gyration
rg = protein.radius_of_gyration()
print(f"Radius of gyration: {rg:.2f} Å")

# Select by residue
active_site = u.select_atoms('resid 35 or resid 52')
```

### Visualization
```bash
# PyMOL
pymol protein_1CRN.pdb

# VMD
vmd protein_1CRN.pdb

# Chimera
chimera protein_1CRN.pdb
```

## Convert to HDF5

```python
from Bio.PDB import PDBParser
import h5py
import numpy as np

# Parse PDB
parser = PDBParser()
structure = parser.get_structure('crambin', 'protein_1CRN.pdb')

# Extract data
atom_names = []
coords = []
bfactors = []

for atom in structure.get_atoms():
    atom_names.append(atom.name)
    coords.append(atom.coord)
    bfactors.append(atom.bfactor)

# Save to HDF5
with h5py.File('crambin.h5', 'w') as f:
    f.create_dataset('coordinates', data=np.array(coords))
    f.create_dataset('atom_names', data=np.array(atom_names, dtype='S4'))
    f.create_dataset('bfactors', data=np.array(bfactors))
    f.attrs['pdb_id'] = '1CRN'
    f.attrs['resolution'] = 0.945
```

## Applications

- **Molecular dynamics:** Starting structures for simulations
- **Drug discovery:** Docking studies and structure-based design
- **Protein folding:** Reference structures for prediction methods
- **Crystallography:** Test cases for refinement algorithms
- **Education:** Teaching protein structure and biochemistry

## References

- **RCSB PDB:** [rcsb.org](https://www.rcsb.org/)
- **Biopython:** [biopython.org](https://biopython.org/)
- **MDAnalysis:** [mdanalysis.org](https://www.mdanalysis.org/)
- **PDB Format Guide:** [wwpdb.org/documentation/file-format](https://www.wwpdb.org/documentation/file-format)
