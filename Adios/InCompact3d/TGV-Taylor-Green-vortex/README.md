# Taylor-Green Vortex (TGV)

This dataset contains Direct Numerical Simulation (DNS) results for the Taylor-Green Vortex case using the Incompact3D solver.

## About Taylor-Green Vortex

The Taylor–Green Vortex (TGV) is a well-known benchmark in computational fluid dynamics that models the transition from an initially laminar state to fully turbulent flow. It's attractive as a test case due to its simple setup with periodic and free-slip boundary conditions.

## Simulation Parameters

- **Reynolds number (Re)**: 1,600
- **Simulation type**: Direct Numerical Simulation (DNS)

## Dataset Contents

- `data.bp5/` - ADIOS2 BP5 output files containing simulation results
- `adios2_config.xml` - ADIOS2 configuration file
- `input.i3d` - Incompact3D input parameters

## References

- [Reference data (Re=1,250 to 20,000)](https://zenodo.org/records/2577239#.YsV6GozMI5k)
- [Incompact3D TGV Documentation](https://xcompact3d.readthedocs.io/en/latest/pages/cases/TGV.html)
- [Incompact3D GitHub](https://github.com/xcompact3d/Incompact3d)