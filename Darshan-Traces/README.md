# Darshan I/O Traces

Darshan is a lightweight I/O characterization tool for HPC applications. It captures detailed information about I/O behavior including file access patterns, bandwidth achieved, and operation counts.

## Available Datasets

### Example Logs (36 MB)
**Location:** [`examples/`](examples/)

Curated collection of example Darshan logs demonstrating various I/O patterns and edge cases. These logs are maintained by the Darshan HPC project and include comprehensive README descriptions for each log.

**Source:** [github.com/darshan-hpc/darshan-logs](https://github.com/darshan-hpc/darshan-logs)

**Contents:**
- IOR POSIX and DAOS benchmark logs
- MPI-IO test logs across multiple Darshan versions
- Imbalanced I/O patterns
- HDF5 and NetCDF application traces
- Parallel compression examples
- STDIO usage patterns

**Tools:** Install Darshan utilities to parse logs: [darshan.readthedocs.io](https://darshan.readthedocs.io/)

---

## Production Trace Collections (Shadow Datasets)

These massive trace collections require significant storage and are provided as references with download instructions.

### ALCF Polaris Darshan Logs (Continuously Updated)
**Size:** 1.2M+ logs (growing ~3,000/day)

The most comprehensive production HPC I/O trace collection, continuously updated from all jobs on Argonne's 560-node Polaris system.

**Access:** [zenodo.org/records/15353810](https://zenodo.org/records/15353810)

**Use Cases:**
- Understanding modern GPU-accelerated HPC I/O patterns
- Storage system performance analysis
- Workload characterization for exascale systems

### Blue Waters Darshan Dataset
**Size:** Multi-TB

Historical traces from NCSA's decommissioned Cray XE6/XK7 system including system monitoring metrics and Lustre performance data.

**Access:** [bluewaters.ncsa.illinois.edu/data-sets](https://bluewaters.ncsa.illinois.edu/data-sets) (via Globus)

### Mira I/O Data Repository
**Size:** Large-scale traces

I/O traces from Argonne's IBM Blue Gene/Q system.

**Access:** [ftp.mcs.anl.gov/pub/darshan/data](ftp://ftp.mcs.anl.gov/pub/darshan/data)

---

## Resources

- **Darshan Documentation:** [darshan.readthedocs.io](https://darshan.readthedocs.io/)
- **Darshan GitHub:** [github.com/darshan-hpc/darshan](https://github.com/darshan-hpc/darshan)
- **PyDarshan (Python interface):** [github.com/darshan-hpc/darshan/tree/main/darshan-util/pydarshan](https://github.com/darshan-hpc/darshan/tree/main/darshan-util/pydarshan)
