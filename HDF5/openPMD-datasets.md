# openPMD Example Datasets

**Format:** HDF5
**Size:** ~286 KB each
**Source:** [github.com/openPMD/openPMD-example-datasets](https://github.com/openPMD/openPMD-example-datasets)

## About openPMD

The Open Standard for Particle-Mesh Data (openPMD) is a format specification for storing particle and mesh-based simulation data. It's used extensively in plasma physics, particle accelerators, and electromagnetic simulations.

## Files

### openPMD_2d_sample.h5
2D particle-mesh simulation data following the openPMD standard.

### openPMD_3d_sample.h5
3D particle-mesh simulation data following the openPMD standard.

## Data Structure

These files contain:
- **Particle data:** Positions, momenta, charges
- **Mesh/field data:** Electric and magnetic fields
- **Metadata:** Simulation parameters, units, iteration info

## Usage

### Python (openPMD-api)
```python
import openpmd_api as io

series = io.Series("openPMD_3d_sample.h5", io.Access.read_only)

for iteration in series.iterations:
    print(f"Iteration {iteration}")
    i = series.iterations[iteration]

    # Access mesh data
    if "E" in i.meshes:
        E = i.meshes["E"]
        print(f"Electric field components: {list(E.keys())}")

    # Access particle data
    if "electrons" in i.particles:
        electrons = i.particles["electrons"]
        print(f"Electron properties: {list(electrons.keys())}")
```

### Python (h5py - raw access)
```python
import h5py

with h5py.File('openPMD_3d_sample.h5', 'r') as f:
    print("Groups:", list(f.keys()))

    # Explore structure
    def print_structure(name, obj):
        print(name)

    f.visititems(print_structure)
```

## Applications

openPMD format is used by:
- **PIConGPU:** Particle-in-Cell GPU simulation
- **WarpX:** Advanced electromagnetic plasma simulator
- **FBPIC:** Fourier-Bessel Particle-In-Cell code
- **OSIRIS:** Plasma simulation framework

## References

- **openPMD Standard:** [github.com/openPMD/openPMD-standard](https://github.com/openPMD/openPMD-standard)
- **openPMD-api:** [openpmd-api.readthedocs.io](https://openpmd-api.readthedocs.io/)
- **Installation:** `pip install openpmd-api`
