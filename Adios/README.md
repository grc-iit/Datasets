# ADIOS2 Datasets

ADIOS2 (Adaptable I/O System) is an HPC I/O framework designed for efficient parallel data movement in scientific simulations.

**Documentation:** [adios2.readthedocs.io](https://adios2.readthedocs.io/)

## Index

### Scientific Applications
- [Gray-Scott](Gray-Scott/) - Reaction-diffusion model simulations
- [InCompact3d](InCompact3d/) - Turbulent flow simulations (DNS/LES)
- [Lammps](Lammps/) - Molecular dynamics simulations
- [WRF](WRF/) - Weather Research and Forecasting model

## Gray-Scott

Reaction-diffusion model simulating chemical species interactions producing complex spatial patterns (spots, stripes, spiral waves).

**Datasets:**
- [L64-noise-variations](Gray-Scott/L64-noise-variations/) - Multiple 64×64 grid simulations with varying noise (0.01-0.16)
- [L48-noise-0.16](Gray-Scott/L48-noise-0.16/) - Single 48×48 grid simulation with noise=0.16

**Resources:** [Gray-Scott Model](https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/)

## InCompact3d

High-order finite-difference flow solvers for turbulent flows using Direct and Large Eddy Simulations (DNS/LES) with Immersed Boundary Method.

**Datasets:**
- [ABL-Atmospheric-Boundary-Layer](InCompact3d/ABL-Atmospheric-Boundary-Layer/) - Turbulent atmospheric boundary layer
- [Cavity](InCompact3d/Cavity/) - Lid-driven cavity flow (classic CFD benchmark)
- [Pipe-Flow](InCompact3d/Pipe-Flow/) - Turbulent pipe flow
- [Periodic-hill](InCompact3d/Periodic-hill/) - Flow over periodic hills
- [TGV-Taylor-Green-vortex](InCompact3d/TGV-Taylor-Green-vortex/) - Taylor-Green Vortex (Re=1,600)

**Resources:** [github.com/xcompact3d/Incompact3d](https://github.com/xcompact3d/Incompact3d)

## Lammps

Large-scale Atomic/Molecular Massively Parallel Simulator for classical molecular dynamics.

**Datasets:**
- [2D-lennard-jones-fluid](Lammps/2D-lennard-jones-fluid/) - 2D Lennard-Jones fluid
- [3D-lennard-jones-fluid](Lammps/3D-lennard-jones-fluid/) - 3D Lennard-Jones fluid
- [GO-nanoparticle-water](Lammps/GO-nanoparticle-water/) - Graphene oxide in water
- [carbon-nanotube-water](Lammps/carbon-nanotube-water/) - Carbon nanotube in water
- [ethanol-water](Lammps/ethanol-water/) - Ethanol-water mixture
- [gold-nanoparticle-water](Lammps/gold-nanoparticle-water/) - Gold nanoparticle in water
- [kalj](Lammps/kalj/) - Kob-Andersen binary mixture
- [nacl-water](Lammps/nacl-water/) - NaCl in water

**Resources:** [lammps.org](https://www.lammps.org/)

## WRF

Weather Research and Forecasting model for atmospheric simulations from meters to thousands of kilometers.

**Datasets:**
- [v4.4_bench_conus12km](WRF/v4.4_bench_conus12km/) - WRF v4.4 benchmark, CONUS 12km domain, 1-hour simulation

**Resources:** [github.com/wrf-model/WRF](https://github.com/wrf-model/WRF) | [NERSC WRF Benchmark](https://docs.nersc.gov/applications/wrf/wrf_benchmark/)
