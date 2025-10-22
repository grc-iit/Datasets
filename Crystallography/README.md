# Crystallography Files (CIF Format)

Crystallographic Information File (CIF) is the standard format for crystallographic data. Developed by the International Union of Crystallography (IUCr), CIF files store crystal structures, atomic coordinates, and experimental metadata for materials science, chemistry, and structural biology.

## Available Datasets

### crambin_1CRN.cif
**Size:** 68 KB
**Source:** Crystallography Open Database (COD)
**Type:** Protein structure
**Resolution:** 0.945 Å (ultra-high resolution)

Crystal structure of crambin, a small plant seed protein. This ultra-high resolution structure showcases the power of X-ray crystallography in determining protein structures at near-atomic detail.

**PDB ID:** 1CRN
**Space Group:** P21
**Unit Cell:** a=40.96 Å, b=18.65 Å, c=22.52 Å, β=90.77°

### quartz_1000000.cif
**Size:** 13 KB
**Source:** Crystallography Open Database
**Mineral:** α-Quartz (SiO₂)
**Structure Type:** Hexagonal

Structure of alpha-quartz, one of the most abundant minerals on Earth. Demonstrates the hexagonal crystal system and provides a standard reference for silicate structures.

**Space Group:** P3221
**Formula:** SiO₂
**Applications:** Reference structure, materials science, geochemistry

### calcite_9008460.cif
**Size:** 4.5 KB
**Source:** Crystallography Open Database
**Mineral:** Calcite (CaCO₃)
**Structure Type:** Rhombohedral

Structure of calcite, the most stable polymorph of calcium carbonate. Important for geology, biomineralization, and carbon sequestration research.

**Space Group:** R-3c
**Formula:** CaCO₃
**Applications:** Biomineralization, carbon cycle, materials engineering

## About CIF Format

CIF files contain:
- **Crystal data:** Unit cell parameters (a, b, c, α, β, γ)
- **Symmetry:** Space group, symmetry operations
- **Atomic positions:** Fractional coordinates, occupancies, thermal parameters
- **Metadata:** Experimental conditions, citations, refinement statistics
- **Structure factors:** (in some files) Observed diffraction intensities

**Key Features:**
- Human-readable text format
- Self-describing with data tags (e.g., `_cell_length_a`)
- Standardized by IUCr
- Used across materials science, chemistry, and biology
- Machine-parseable with well-defined grammar

## Installation

```bash
# Python - pymatgen (comprehensive materials analysis)
pip install pymatgen

# Python - ASE (Atomic Simulation Environment)
pip install ase

# Python - gemmi (fast CIF/PDB parser for crystallography)
pip install gemmi

# Python - PyCifRW (pure CIF reader/writer)
pip install PyCifRW

# Command-line tools
# cif2cell - converts CIF to various simulation formats
pip install cif2cell
```

## Usage Examples

### Python (pymatgen - Materials Science)
```python
from pymatgen.io.cif import CifParser
import matplotlib.pyplot as plt
from pymatgen.analysis.diffraction.xrd import XRDCalculator

# Parse CIF file
parser = CifParser("quartz_1000000.cif")
structure = parser.get_structures()[0]

# Print structure information
print(f"Formula: {structure.composition.reduced_formula}")
print(f"Space group: {structure.get_space_group_info()}")
print(f"Lattice:\n{structure.lattice}")
print(f"Volume: {structure.volume:.2f} Ų")
print(f"Density: {structure.density:.2f} g/cm³")

# Print atomic positions
print("\nAtomic sites:")
for site in structure:
    print(f"{site.species_string:3s} {site.frac_coords}")

# Calculate powder XRD pattern
calculator = XRDCalculator()
pattern = calculator.get_pattern(structure)

# Plot XRD pattern
plt.figure(figsize=(10, 6))
plt.plot(pattern.x, pattern.y)
plt.xlabel("2θ (degrees)")
plt.ylabel("Intensity (a.u.)")
plt.title(f"Calculated XRD Pattern: {structure.composition.reduced_formula}")
plt.grid(True, alpha=0.3)
plt.savefig("xrd_pattern.png", dpi=150)
```

### Python (ASE - Atomic Simulations)
```python
from ase.io import read, write
from ase.visualize import view
import numpy as np

# Read CIF file
atoms = read("calcite_9008460.cif")

# Print structure info
print(f"Chemical formula: {atoms.get_chemical_formula()}")
print(f"Number of atoms: {len(atoms)}")
print(f"Cell volume: {atoms.get_volume():.2f} ų")
print(f"Cell parameters: {atoms.cell.cellpar()}")

# Get atomic positions
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()
print("\nAtomic coordinates (Cartesian, Å):")
for symbol, pos in zip(symbols, positions):
    print(f"{symbol:3s} {pos[0]:8.3f} {pos[1]:8.3f} {pos[2]:8.3f}")

# Create supercell
supercell = atoms.repeat((2, 2, 2))
print(f"\nSupercell atoms: {len(supercell)}")

# Visualize (if GUI available)
# view(atoms)

# Export to other formats
write("calcite.xyz", atoms)          # XYZ format
write("calcite_poscar.vasp", atoms)  # VASP POSCAR
write("calcite.pdb", atoms)          # PDB format
```

### Python (gemmi - Fast CIF Parser)
```python
import gemmi

# Read CIF file
doc = gemmi.cif.read_file("crambin_1CRN.cif")
block = doc.sole_block()

# Access CIF data items
print(f"Structure name: {block.name}")
print(f"Cell a: {block.find_value('_cell_length_a')}")
print(f"Cell b: {block.find_value('_cell_length_b')}")
print(f"Cell c: {block.find_value('_cell_length_c')}")
print(f"Space group: {block.find_value('_symmetry_space_group_name_H-M')}")

# Parse to structure
structure = gemmi.read_structure("crambin_1CRN.cif")
print(f"\nModel: {structure.name}")
print(f"Number of models: {len(structure)}")

# Access atoms
model = structure[0]
for chain in model:
    for residue in chain:
        print(f"Residue: {residue.name}")
        for atom in residue:
            print(f"  {atom.name:4s} {atom.pos.x:7.3f} {atom.pos.y:7.3f} {atom.pos.z:7.3f}")
        break  # Just show first residue
    break
```

### Pure CIF Parsing (PyCifRW)
```python
from CifFile import ReadCif

# Read CIF
cf = ReadCif("quartz_1000000.cif")
block = cf[list(cf.keys())[0]]

# Extract cell parameters
a = float(block['_cell_length_a'])
b = float(block['_cell_length_b'])
c = float(block['_cell_length_c'])
alpha = float(block['_cell_angle_alpha'])
beta = float(block['_cell_angle_beta'])
gamma = float(block['_cell_angle_gamma'])

print(f"Unit cell: a={a:.3f} b={b:.3f} c={c:.3f}")
print(f"Angles: α={alpha:.2f}° β={beta:.2f}° γ={gamma:.2f}°")

# Extract atomic sites
labels = block['_atom_site_label']
x = block['_atom_site_fract_x']
y = block['_atom_site_fract_y']
z = block['_atom_site_fract_z']

print("\nAtomic positions (fractional):")
for label, fx, fy, fz in zip(labels, x, y, z):
    print(f"{label:4s} {float(fx):7.4f} {float(fy):7.4f} {float(fz):7.4f}")
```

## Command-Line Tools

```bash
# View CIF info
gemmi cif2json quartz_1000000.cif | head -20

# Convert CIF to other formats
# To XYZ
cif2cell quartz_1000000.cif -p xyz -o quartz.xyz

# To VASP POSCAR
cif2cell quartz_1000000.cif -p vasp -o POSCAR

# Validate CIF file
gemmi validate calcite_9008460.cif

# Extract specific data items
grep "_cell_length" quartz_1000000.cif
grep "_atom_site" quartz_1000000.cif | head
```

## Structure Visualization

### Python (matplotlib + ASE)
```python
from ase.io import read
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read structure
atoms = read("quartz_1000000.cif")

# Create 3D plot
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot atoms
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()

colors = {'Si': 'gold', 'O': 'red', 'C': 'gray', 'Ca': 'green'}
for symbol in set(symbols):
    mask = [s == symbol for s in symbols]
    pos = positions[mask]
    ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2],
               c=colors.get(symbol, 'blue'),
               label=symbol, s=100)

ax.set_xlabel('X (Å)')
ax.set_ylabel('Y (Å)')
ax.set_zlabel('Z (Å)')
ax.legend()
plt.title(atoms.get_chemical_formula())
plt.savefig("structure_3d.png", dpi=150)
```

### Using nglview (Jupyter)
```python
import nglview as nv
from ase.io import read

# Read and visualize
atoms = read("crambin_1CRN.cif")
view = nv.show_ase(atoms)
view
```

## Structure Analysis

### Calculate Bond Lengths and Angles
```python
from pymatgen.io.cif import CifParser
from pymatgen.analysis.local_env import CrystalNN

parser = CifParser("calcite_9008460.cif")
structure = parser.get_structures()[0]

# Nearest neighbor analysis
nn = CrystalNN()
for i, site in enumerate(structure):
    neighbors = nn.get_nn_info(structure, i)
    print(f"\nSite {i} ({site.species_string}):")
    for neighbor in neighbors:
        distance = neighbor['weight']
        site_idx = neighbor['site_index']
        neighbor_species = structure[site_idx].species_string
        print(f"  -> {neighbor_species} at {distance:.3f} Å")
```

### Generate Supercell and Export
```python
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar

# Read structure
parser = CifParser("quartz_1000000.cif")
structure = parser.get_structures()[0]

# Create supercell
supercell = structure.copy()
supercell.make_supercell([2, 2, 2])

print(f"Original: {len(structure)} atoms")
print(f"Supercell: {len(supercell)} atoms")

# Export to VASP format
poscar = Poscar(supercell)
poscar.write_file("POSCAR_supercell")
```

## Applications

CIF files are used in:
- **Materials science:** Crystal structure databases, structure prediction
- **Chemistry:** Small molecule crystallography, molecular structures
- **Structural biology:** Protein and nucleic acid structures
- **Mineralogy:** Mineral identification and classification
- **Solid-state physics:** Electronic structure calculations, phonons
- **Drug discovery:** Protein-ligand complexes, polymorphism studies

## Major Databases

- **Crystallography Open Database (COD):** 500,000+ structures - [crystallography.net](http://www.crystallography.net/)
- **Protein Data Bank (PDB):** 200,000+ biomolecular structures - [rcsb.org](https://www.rcsb.org/)
- **Cambridge Structural Database (CSD):** 1,200,000+ organic/organometallic
- **Inorganic Crystal Structure Database (ICSD):** 250,000+ inorganic structures
- **Materials Project:** 150,000+ computed structures - [materialsproject.org](https://materialsproject.org/)

## Shadow Datasets

For large-scale crystallography databases, see:
- [Shadow-Datasets/Materials-Science.md](../Shadow-Datasets/Materials-Science.md)

## Resources

- **CIF Standard:** [iucr.org/resources/cif](https://www.iucr.org/resources/cif)
- **pymatgen:** [pymatgen.org](https://pymatgen.org/)
- **ASE:** [wiki.fysik.dtu.dk/ase](https://wiki.fysik.dtu.dk/ase/)
- **gemmi:** [gemmi.readthedocs.io](https://gemmi.readthedocs.io/)
- **COD:** [crystallography.net](http://www.crystallography.net/)
- **PDB:** [rcsb.org](https://www.rcsb.org/)
