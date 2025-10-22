# Pipe Flow

This dataset contains simulation results for turbulent pipe flow using the Incompact3D solver.

## About Pipe Flow

When a fluid flows through a pipe, the internal roughness of the pipe wall can create local eddy currents within the fluid, adding resistance to flow. The velocity profile in a pipe shows that fluid at the center of the stream moves more quickly than fluid towards the edge, creating friction between layers within the fluid.

Pipes with smooth walls (glass, copper, brass, polyethylene) have only a small effect on frictional resistance, while pipes with rough walls (concrete, cast iron, steel) create larger eddy currents with more significant effects. Fluids with high viscosity flow more slowly and generally don't support eddy currents, resulting in laminar flow conditions.

## Dataset Contents

- `data.bp5/` - ADIOS2 BP5 output files containing simulation results
- `adios2_config.xml` - ADIOS2 configuration file
- `input.i3d` - Incompact3D input parameters

## References

- [Incompact3D Pipe Flow Example](https://github.com/xcompact3d/Incompact3d/tree/master/examples/Pipe-Flow)
- [Incompact3D Documentation](https://xcompact3d.readthedocs.io/) 