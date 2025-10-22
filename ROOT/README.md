# ROOT Files - Particle Physics Data

ROOT is a data analysis framework developed at CERN for particle physics. It's the standard format for storing collision data from the Large Hadron Collider (LHC) and other particle physics experiments.

## Available Datasets

### higgs_sample.root
**Size:** 20 MB
**Source:** ROOT Project / CERN
**Analysis:** Higgs to Tau-Tau decay channel
**Dataset:** Reduced analysis sample from H→ττ search

Real physics analysis data from searches for the Higgs boson decaying to tau lepton pairs. This is a reduced dataset from the GluGlu→H→TauTau analysis, containing TTrees with reconstructed physics objects and analysis variables.

**Access:** [root.cern/files](https://root.cern/files/)

### hsimple_tutorial.root
**Size:** 506 KB
**Source:** ROOT Project Official Tutorials

Classic ROOT tutorial file containing example histograms, trees, and other ROOT objects. Used for learning ROOT basics and testing ROOT installations.

## About ROOT Format

ROOT files store:
- **TTree:** Columnar data structures (like tables with branches)
- **TH1/TH2/TH3:** 1D/2D/3D histograms
- **TGraph:** Graphs and plots
- **Custom classes:** User-defined C++ objects
- **Metadata:** Analysis metadata and configuration

**Key Features:**
- Columnar storage for efficient queries
- Built-in compression (typically ZLIB or LZMA)
- Self-describing format
- C++ integration with code generation
- Multi-TB file support

## Installation

```bash
# Install ROOT via conda (recommended)
conda install -c conda-forge root

# Or use Docker
docker pull rootproject/root

# Or build from source
# See: https://root.cern/install/
```

## Usage Examples

### Python (PyROOT)
```python
import ROOT

# Open ROOT file
f = ROOT.TFile.Open("hsimple_tutorial.root")

# List contents
f.ls()

# Get a histogram
h = f.Get("hpx")
print(f"Histogram entries: {h.GetEntries()}")
print(f"Mean: {h.GetMean():.2f}")
print(f"RMS: {h.GetRMS():.2f}")

# Draw histogram
c = ROOT.TCanvas("c", "Canvas", 800, 600)
h.Draw()
c.SaveAs("histogram.png")

# Access TTree
tree = f.Get("ntuple")
if tree:
    print(f"Tree entries: {tree.GetEntries()}")

    # Print branch names
    branches = tree.GetListOfBranches()
    for b in branches:
        print(f"Branch: {b.GetName()}")

    # Loop over events
    for i, event in enumerate(tree):
        if i >= 10:  # First 10 events
            break
        print(f"Event {i}: px={event.px:.2f}, py={event.py:.2f}")

f.Close()
```

### Python (uproot - Pure Python, No ROOT Required)
```python
import uproot
import matplotlib.pyplot as plt
import numpy as np

# Open ROOT file (no ROOT installation needed!)
file = uproot.open("hsimple_tutorial.root")

# List contents
print("Keys:", file.keys())

# Read histogram
h = file["hpx"]
counts, edges = h.to_numpy()

# Plot with matplotlib
plt.figure(figsize=(10, 6))
plt.stairs(counts, edges, fill=True, alpha=0.7)
plt.xlabel("Value")
plt.ylabel("Counts")
plt.title("Histogram from ROOT file")
plt.grid(True, alpha=0.3)
plt.savefig("uproot_histogram.png")

# Read TTree as arrays
tree = file["ntuple"]
arrays = tree.arrays(["px", "py", "pz"], library="np")
print(f"px: {arrays['px'][:5]}")
print(f"py: {arrays['py'][:5]}")

# Or as pandas DataFrame
df = tree.arrays(library="pd")
print(df.head())
print(df.describe())

# Scatter plot
plt.figure(figsize=(8, 8))
plt.scatter(df['px'], df['py'], alpha=0.5)
plt.xlabel("px")
plt.ylabel("py")
plt.title("Momentum Components")
plt.grid(True, alpha=0.3)
plt.savefig("momentum_scatter.png")
```

### Python (awkward-array for Complex Structures)
```python
import uproot
import awkward as ak

# Read Higgs analysis data
file = uproot.open("higgs_sample.root")

# List available trees
print("Trees:", file.keys())

# Access analysis tree
tree = file["GluGluToHToTauTau"]
branches = tree.arrays(library="ak")

# Explore event structure
print(f"Total events: {len(branches)}")
print(f"Available branches: {tree.keys()}")
```

### C++ (Native ROOT)
```cpp
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <iostream>

void analyze_root() {
    // Open file
    TFile *f = TFile::Open("hsimple_tutorial.root");

    // Get histogram
    TH1F *h = (TH1F*)f->Get("hpx");
    std::cout << "Mean: " << h->GetMean() << std::endl;

    // Get tree
    TTree *tree = (TTree*)f->Get("ntuple");
    Float_t px, py, pz;
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);

    // Loop over entries
    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        std::cout << "Event " << i << ": px=" << px << std::endl;
    }

    f->Close();
}
```

### Command Line (rootls, rootprint)
```bash
# List contents
rootls hsimple_tutorial.root

# Print tree structure
rootls -t hsimple_tutorial.root

# Print first entries
root -l hsimple_tutorial.root -e "ntuple->Show(0)"

# Convert to CSV
root -l -b -q 'hsimple_tutorial.root' <<EOF
.! mkdir -p output
ntuple->Scan("*", "", "", 100, 0)
.q
EOF
```

## Higgs Analysis Example

For the Higgs physics analysis data:

```python
import uproot
import matplotlib.pyplot as plt
import numpy as np

# Open Higgs data
f = uproot.open("higgs_sample.root")

# Get the analysis tree
tree = f["GluGluToHToTauTau"]
print("Branches:", tree.keys())

# Read analysis variables as pandas DataFrame
df = tree.arrays(library="pd")

print(f"Total events: {len(df)}")
print(f"Columns: {df.columns.tolist()}")
print(df.head())

# Example: Plot distributions
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot various kinematic distributions
for ax, col in zip(axes.flat, df.columns[:4]):
    ax.hist(df[col], bins=50, alpha=0.7)
    ax.set_xlabel(col)
    ax.set_ylabel("Events")
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("higgs_distributions.png")
```

## Convert ROOT to Other Formats

```python
import uproot
import h5py
import pandas as pd

# ROOT to HDF5
with uproot.open("hsimple_tutorial.root") as f:
    tree = f["ntuple"]
    arrays = tree.arrays(library="np")

    with h5py.File("converted.h5", "w") as hf:
        for name, arr in arrays.items():
            hf.create_dataset(name, data=arr)

# ROOT to Parquet
df = tree.arrays(library="pd")
df.to_parquet("converted.parquet")

# ROOT to CSV
df.to_csv("converted.csv", index=False)
```

## Applications

ROOT files are used for:
- **Particle physics:** LHC experiments (ATLAS, CMS, LHCb, ALICE)
- **Nuclear physics:** Experiments worldwide
- **Astrophysics:** Gamma-ray and cosmic ray experiments
- **Detector simulations:** Geant4 output
- **Data analysis:** Statistical analysis in physics

## Shadow Dataset: CERN Open Data

For multi-TB particle physics datasets, see:
- [Shadow-Datasets/Astrophysics-Cosmology.md](../Shadow-Datasets/Astrophysics-Cosmology.md#particle-physics-cern-open-data)

Full datasets available:
- **CMS:** 2+ PB of collision data
- **ATLAS:** Petabyte-scale datasets
- Includes Monte Carlo simulations and detector data

## Resources

- **ROOT Documentation:** [root.cern](https://root.cern/)
- **ROOT Tutorials:** [root.cern/doc/master/group__Tutorials.html](https://root.cern/doc/master/group__Tutorials.html)
- **uproot:** [uproot.readthedocs.io](https://uproot.readthedocs.io/)
- **CMS Open Data:** [opendata.cern.ch](http://opendata.cern.ch/)
- **ATLAS Open Data:** [opendata.atlas.cern](http://opendata.atlas.cern/)
