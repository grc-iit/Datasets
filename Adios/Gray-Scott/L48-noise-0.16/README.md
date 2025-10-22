# Gray-Scott Reaction-Diffusion Model - L48 High Noise

This dataset contains a Gray-Scott reaction-diffusion simulation output with high noise level and a 48x48 grid size.

## Configuration

- **Diffusion coefficient Du**: 0.2
- **Diffusion coefficient Dv**: 0.1
- **Feed rate F**: 0.01
- **Kill rate k**: 0.05
- **Grid size L**: 48
- **Noise**: 0.16 (high)
- **Time step dt**: 2.0
- **Steps**: 40
- **Processes**: 16
- **Engine**: ADIOS2 BP5

## Dataset Contents

- `data.bp5/` - ADIOS2 BP5 output files
- `config.txt` - Full configuration parameters used for this simulation

## About Gray-Scott

The Gray-Scott model is a reaction-diffusion system that models the interaction between two chemical species. It's commonly used as a benchmark for I/O systems and visualization tools in computational science.

## References

- [Gray-Scott Model Information](https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/)
