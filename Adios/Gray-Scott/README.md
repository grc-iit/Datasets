# Gray-Scott Reaction-Diffusion Datasets

This directory contains datasets generated using the Gray-Scott reaction-diffusion model, a classic system in computational science used for benchmarking I/O performance and testing visualization tools.

## About Gray-Scott

The Gray-Scott model simulates the interaction between two chemical species through reaction and diffusion processes. The model produces complex spatial patterns including spots, stripes, and spiral waves depending on the parameters used.

The governing equations are:

```
∂u/∂t = Du∇²u - uv² + F(1-u)
∂v/∂t = Dv∇²v + uv² - (F+k)v
```

Where:
- **Du, Dv**: Diffusion coefficients for species u and v
- **F**: Feed rate
- **k**: Kill rate

## Available Datasets

### L64-noise-variations
Multiple simulations with a 64x64 grid and varying noise levels (0.01, 0.02, 0.04, 0.08, 0.16) to study how initial perturbations affect pattern formation.

### L48-noise-0.16
Single simulation with a 48x48 grid and high noise level (0.16).

## Data Format

All datasets are stored using the ADIOS2 BP5 engine, which provides:
- High-performance parallel I/O
- Self-describing data format
- Built-in compression support
- Random access capabilities

## References

- [Gray-Scott Model Information](https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/)
- [ADIOS2 Documentation](https://adios2.readthedocs.io/)
