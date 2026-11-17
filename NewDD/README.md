# NewDD: RAMSES Z-Directional Slab Decomposition Package

A tool for reading RAMSES simulation snapshot data and decomposing the domain into slabs along the Z-axis direction. Designed for efficient post-processing of large-scale cosmological simulation data.

---

## Related Documentation

| Document | Description |
|----------|-------------|
| [QUICKSTART.md](QUICKSTART.md) | Get started in 5 minutes |
| [CONFIGURATION.md](CONFIGURATION.md) | Advanced configuration options |
| [PIPELINE.md](PIPELINE.md) | PGalF pipeline integration guide |
| [README4GADGET.md](README4GADGET.md) | GADGET HDF5 adaptation guide |
| [SWIFT_GADGET.md](SWIFT_GADGET.md) | Implemented SWIFT/GADGET HDF5 usage guide |
| [NOTE_SWIFT_h_units.md](NOTE_SWIFT_h_units.md) | SWIFT has no h in its units: how NewDD gets h, and the pre-fix h=1 bug |

---

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Directory Structure](#directory-structure)
4. [Dependencies](#dependencies)
5. [Build Instructions](#build-instructions)
6. [Usage](#usage)
7. [Input File Formats](#input-file-formats)
8. [Output File Formats](#output-file-formats)
9. [Core Data Structures](#core-data-structures)
10. [Algorithm Details](#algorithm-details)
11. [Compilation Options](#compilation-options)
12. [Source File Descriptions](#source-file-descriptions)
13. [Performance and Memory](#performance-and-memory)
14. [Examples](#examples)
15. [Physical Constants](#physical-constants)
16. [Troubleshooting](#troubleshooting)

---

## Overview

**NewDD** (New Domain Decomposition) is a tool for processing RAMSES (Rasterized N-body Adaptive Refinement Tree) simulation output data. It decomposes the simulation domain into vertical slabs along the Z-axis, reorganizing data into a format optimized for parallel analysis and Friends-of-Friends (FoF) galaxy finding algorithms.

### Key Features

- **Z-directional slab decomposition**: Uniformly partitions the simulation box along the Z-axis
- **MPI parallelization**: Efficient parallel I/O for large-scale datasets
- **Multi-component processing**: Simultaneous handling of dark matter, stars, gas cells, and sinks (black holes)
- **Unit conversion**: Automatic conversion from simulation units to physical units (cMpc/h, km/s, Msun/h)
- **AMR support**: Extraction of leaf (finest resolution) cells from Adaptive Mesh Refinement structures

---

## Features

### Data Reading
- RAMSES info files (simulation metadata)
- AMR mesh structures (Adaptive Mesh Refinement)
- Particle data (dark matter + stars)
- Hydrodynamics data (gas density, velocity, temperature)
- Sink particles (black holes/point masses)

### Data Processing
- Leaf gas cell extraction
- Sorting particles/cells by Z-coordinate
- Domain decomposition into user-specified number of slabs
- Unit conversion (code units → physical units)

### Data Output
- Separate binary files for each slab
- Metadata info file for each slab

---

## Directory Structure

```
NewDD/
├── Makefile            # Build configuration
├── README.md           # This document
├── ramses.h            # Core data structure definitions
├── Memory.h            # Memory manager header
├── newdd.c             # Main program (entry point)
├── rd_info.c           # Metadata reader
├── rd_amr.c            # AMR mesh reader
├── rd_part.c           # Particle data reader
├── rd_hydro.c          # Hydrodynamics data reader
├── rd_sink.c           # Sink particle reader
├── rd_grav.c           # Gravity/potential reader (optional)
├── find_leaf_gas.c     # Leaf gas cell extractor
├── utils.c             # Utility functions & slab decomposition
├── header.c            # Header file I/O
└── Memory2.c           # Custom memory allocator
```

---

## Dependencies

### Required
- **C Compiler**: GCC or compatible compiler
- **MPI Library**: OpenMPI, MPICH, or compatible implementation (for parallel execution)

### Standard Libraries
- `stdio.h`, `stdlib.h`, `string.h`
- `math.h`
- `unistd.h`, `sys/stat.h`

### Input Data
- RAMSES simulation output snapshots (`output_XXXXX/` directories)

---

## Build Instructions

### Using configure (Recommended)

The recommended way to build NewDD is through the top-level `configure` script, which generates the Makefile from `Makefile.in`:

```bash
# From the top-level GalaxyFinder directory:
./configure                     # Default settings
./configure --debug             # Debug build
./configure --newdd-nmeg=30000  # Custom NMEG

# Build NewDD only
make newdd

# Or build everything
make all
```

#### Key configure options for NewDD

| Option | Default | Description |
|--------|---------|-------------|
| `--cc=CC` | `mpiicx` | C compiler |
| `--opt=FLAGS` | `-O3` | Optimization flags |
| `--debug` | - | Use `-g` debug flags |
| `--newdd-nmeg=N` | `20000` | Memory pool size (MB) |
| `--newdd-nchem=N` | `9` | Number of chemical species |
| `--newdd-ndust=N` | `4` | Number of dust species |
| `--newdd-wgroupsize=N` | `5` | I/O sync group size |
| `--nener=N` | `0` | Number of radiation energy groups |
| `--npre=N` | `8` | Precision indicator |

Run `./configure --help` for the full list of options.

### Direct Build (Alternative)

You can also build directly within the NewDD directory if the Makefile already exists:

```bash
cd NewDD
make clean
make all
```

### Build Targets

| Target | Description |
|--------|-------------|
| `make` or `make default` | Build libmyram.a library only |
| `make all` | Build library + newdd.exe executable |
| `make newdd` | Rebuild newdd.exe only |
| `make clean` | Remove all objects, archives, and executables |

---

## Usage

### Command Format

```bash
mpirun -np <num_processes> ./newdd.exe <ISTEP> <NSPLIT>
```

### Parameters

| Parameter | Description | Examples |
|-----------|-------------|----------|
| `<num_processes>` | Number of MPI processes (must be a divisor of ncpu) | 8, 16, 32 |
| `<ISTEP>` | Snapshot number (XXXXX from output_XXXXX) | 5, 10, 100 |
| `<NSPLIT>` | Number of slab divisions | 32, 64, 128 |

### Execution Examples

```bash
# Single processor execution
./newdd.exe 5 64

# Parallel execution with 8 MPI processes
mpirun -np 8 ./newdd.exe 5 64

# Execution on SLURM cluster
srun -n 32 ./newdd.exe 100 128
```

### Input Path Configuration

By default, the program looks for `output_XXXXX/` directories in the parent directory.
To change the path, modify the path settings in `rd_info.c`.

---

## Input File Formats

### Required RAMSES Output Files

NewDD reads the following RAMSES output files:

```
output_XXXXX/
├── info_XXXXX.txt          # Simulation metadata (ASCII)
├── amr_XXXXX.out00001      # AMR mesh structure (Fortran binary)
├── amr_XXXXX.out00002
├── ...
├── amr_XXXXX.outNCPU
├── part_XXXXX.out00001     # Particle data (Fortran binary)
├── part_XXXXX.out00002
├── ...
├── hydro_XXXXX.out00001    # Hydrodynamics data (Fortran binary)
├── hydro_XXXXX.out00002
├── ...
└── sink_XXXXX.out00001     # Sink data (Fortran binary, single file)
```

### Info File Format (ASCII)

```
ncpu                = 512
ndim                = 3
levelmin            = 4
nlevelmax           = 13
ngridmax            = 262144
nstep_coarse        = 10000
boxlen              = 100.0
time                = 1.0
aexp                = 1.0
H0                  = 68.0
omega_m             = 0.3
omega_l             = 0.7
omega_k             = 0.0
omega_b             = 0.05
unit_l              = 3.085678e+24
unit_d              = 1.0e-30
unit_t              = 1.0e+17
ordering type       = quadhilbert
```

### Fortran Binary Format

All binary files use Fortran unformatted sequential I/O format:
- Each record is preceded and followed by 4-byte size information
- Read using the `F77read()` macro

### Diagnostic Macros (ramses.h)

| Macro | Output | Condition | Description |
|-------|--------|-----------|-------------|
| `DEBUGPRINT(fmt, ...)` | stderr | `DEBUG` is true | Debug messages with file/line info |
| `LOGPRINT(fmt, ...)` | stdout | `LOG` is true | Log messages with file/line info |
| `ERRORPRINT(fmt, ...)` | stderr | Always | Error messages with file/line info |
| `ERRORPRINT0(fmt)` | stderr | Always | Error messages (no format args) |
| `YES` / `NO` | - | - | Boolean constants (1/0) |

---

## Output File Formats

### Output Directory Structure

```
FoF_Data/NewDD.XXXXX/
├── SN.XXXXX.00000.info         # Slab 0 metadata
├── SN.XXXXX.DM.00000.dat       # Slab 0 dark matter
├── SN.XXXXX.STAR.00000.dat     # Slab 0 stars
├── SN.XXXXX.GAS.00000.dat      # Slab 0 gas
├── SN.XXXXX.SINK.00000.dat     # Slab 0 sinks
├── SN.XXXXX.00001.info         # Slab 1 metadata
├── SN.XXXXX.DM.00001.dat       # Slab 1 dark matter
├── ...
└── SN.XXXXX.SINK.NNNNN.dat     # Last slab sinks
```

### Metadata Files (.info)

ASCII files containing information for each slab:
- Number of particles/cells in the slab
- Slab boundaries (xmin, xmax)
- Cosmological parameters

### Binary Data Files (.dat)

Each file stores an array of the corresponding data structure directly:

```c
// Reading example
DmType *dm = malloc(sizeof(DmType) * ndm);
fread(dm, sizeof(DmType), ndm, fp);
```

---

## Core Data Structures

### RamsesType (Main Container)

```c
typedef struct RamsesType {
    // Particle/cell counts
    int npart;              // Total particle count (DM + stars)
    int ndm, nstar;         // Dark matter, star counts
    int nsink;              // Sink count
    int ngas;               // Gas cell count

    // Simulation parameters
    int ncpu, ndim;         // Number of CPU domains, dimensions
    int nlevelmax;          // Maximum AMR level

    // Cosmological parameters
    dptype omega_m, omega_l, omega_b, omega_k;
    dptype H0, aexp, boxlen_ini;

    // Unit conversion factors
    dptype scale_l, scale_d, scale_t, scale_v;
    dptype kmscale_v;       // km/s conversion
    dptype mpcscale_l;      // cMpc/h conversion

    // Hydro variable indices
    int nvar, nboundary, ncoarse, nener;
    int imetal, ichem, idust, iion, nelt, inener;

    // Data arrays
    PmType *particle;       // DM + stars
    GasType *gas;           // Gas cells
    SinkType *sink;         // Sinks
    MeshType mesh;          // AMR structure
    HydroType hydro;        // Hydrodynamics data

    // Output range
    int xmin, xmax;         // Slab boundaries
} RamsesType;
```

### PmType (Particles: Dark Matter/Stars)

```c
typedef struct PmType {
    dptype x, y, z;         // Position [cMpc/h]
    dptype vx, vy, vz;      // Velocity [km/s]
    dptype mass;            // Mass [Msun/h]
    idtype id;              // Unique ID
    int levelp;             // AMR level
    familytype family, tag; // Type: 1=DM, 2=star; tag
#ifdef OUTPUT_PARTICLE_POTENTIAL
    dptype potent;          // Gravitational potential
#endif
#ifndef NBODY
    dptype tp;              // Birth time (stars)
    dptype zp;              // Metallicity (stars)
#ifdef NCHEM
    dptype chem[NCHEM];     // Chemical abundances
#endif
    dptype mass0;           // Initial mass [Msun/h]
    dptype birth_d;         // Dust birth parameter
    int partp;              // Particle pointer/index
#endif
} PmType;

typedef PmType DmType;      // Dark matter (family=1)
typedef PmType StarType;    // Stars (family=2)
```

### GasType (Gas Cells)

```c
typedef struct GasType {
    dptype x, y, z;         // Position [cMpc/h]
    dptype cellsize;        // Cell size [cMpc/h]
    float vx, vy, vz;       // Velocity [km/s]
    dptype den;             // Density [simulation units]
    float temp;             // Temperature [K/mu]
    float metallicity;      // Metallicity fraction
#ifdef NCHEM
    float chem[NCHEM];      // Chemical composition
#endif
#ifdef NDUST
    float dust[NDUST];      // Dust species mass fractions
#endif
    float mass;             // Cell mass [Msun/h]
    dptype potent;          // Gravitational potential
    dptype fx, fy, fz;      // Gravitational force [optional]
} GasType;
```

### SinkType (Sinks/Black Holes)

```c
typedef struct SinkType {
    dptype x, y, z;         // Position [cMpc/h]
    dptype vx, vy, vz;      // Velocity [km/s]
    dptype mass;            // Mass [Msun/h]
    dptype tbirth;          // Birth time
    dptype Jx, Jy, Jz;      // Angular momentum
    dptype Sx, Sy, Sz;      // Spin vector
    dptype dMsmbh;          // SMBH accretion rate
    dptype dMBH_coarse;     // Coarse accretion rate
    dptype dMEd_coarse;     // Eddington accretion rate
    dptype Esave, Smag, eps;// AGN parameters
    int id;                 // Unique ID
} SinkType;
```

---

## Algorithm Details

### Overall Workflow

```
Input: RAMSES snapshot (output_XXXXX/)
  │
  ├─ 1. Read metadata (rd_info)
  │      └─ Cosmological parameters, compute unit conversion factors
  │
  ├─ 2. Loop over each CPU domain [1, ncpu]:
  │      ├─ a) Read AMR mesh (rd_amr)
  │      │      └─ Build spatial index, cell hierarchy
  │      │
  │      ├─ b) Read particles (rd_part)
  │      │      └─ Filter DM (family=1), stars (family=2)
  │      │      └─ Apply unit conversions
  │      │
  │      ├─ c) Read hydrodynamics (rd_hydro)
  │      │      └─ Extract density, velocity, temperature for each cell
  │      │
  │      ├─ d) Read sinks (rd_sink, CPU#1 only)
  │      │
  │      ├─ e) Extract leaf gas cells (find_leaf_gas)
  │      │      └─ Extract only finest resolution gas cells
  │      │      └─ Compute cell masses and properties
  │      │
  │      ├─ f) Sort by Z-coordinate
  │      │      └─ dmsortx(), starsortx(), gassortx(), sinksortx()
  │      │
  │      └─ g) Slab decomposition and output (SplitDump)
  │
  └─ Output: FoF_Data/NewDD.XXXXX/
             └─ Per-slab DM, STAR, GAS, SINK files
```

### Slab Decomposition Algorithm (SplitDump)

```c
1. Divide Z-axis into nsplit bins:
   step = (xmax - xmin) / nsplit
   xpos[i] = xmin + i * step

2. For each object (particle/gas cell):
   ibin = (x - xmin) / step
   if xmin <= x < xmax: add to bin ibin

3. MPI synchronization (if MPI enabled):
   - Gather counts from all ranks
   - Compute write offsets
   - Pre-allocate file space
   - Synchronize rank groups (WGROUPSIZE)

4. Parallel write:
   for rank in group:
       for ibin in range(nsplit):
           if rank has data in bin:
               seek to offset[rank][ibin]
               write objects

5. Output: SN.XXXXX.[TYPE].[IBIN].dat
```

### Unit Conversion

```
Code units → Physical units:

Position:  x_phys = x_code * boxlen * scale_l / Mpc * h
           Result: cMpc/h (comoving Mpc per h)

Velocity:  v_phys = v_code * scale_l / scale_t / 1e5
           Result: km/s

Mass:      M_phys = M_code * scale_d * scale_l^3 / Msun * h
           Result: Msun/h (Solar masses per h)

Density:   rho_phys = rho_code * scale_d
           Result: g/cm^3
```

### Leaf Cell Extraction (find_leaf_gas)

```
For each AMR level:
    For each grid in level:
        For each of 2^3 child cells:
            if son[cell] == 0 (no children):
                Extract cell properties
                mass = density * volume * scale_d
                Store in gas[] array
```

---

## Compilation Options

### Makefile Flags

These flags are automatically set by `./configure`. To change them, use the corresponding configure options (e.g., `--newdd-nmeg=30000`):

```makefile
OPT = -O3 -DNENER=0 -DNPRE=8 -DNMEG=20000 -DWGROUPSIZE=5 \
      -DUSE_MPI -DQUADHILBERT -DREAD_SINK -DNCHEM=9 -DNDUST=4
```

### Option Descriptions

| Option | Default | Description |
|--------|---------|-------------|
| `-O3` | Enabled | Optimization level |
| `-DUSE_MPI` | Defined | Enable MPI parallelization |
| `-DQUADHILBERT` | Defined | Use long double precision for Hilbert curves |
| `-DREAD_SINK` | Defined | Enable sink particle reading |
| `-DNENER=0` | 0 | Number of radiation energy groups |
| `-DNMEG=20000` | 20000 | Memory pool size (MB) |
| `-DWGROUPSIZE=5` | 5 | I/O synchronization group size |
| `-DNCHEM=9` | 9 | Number of chemical species |
| `-DNDUST=4` | 4 | Number of dust species |
| `-DNPRE=8` | 8 | Precision indicator |

### Custom Builds

Use configure options for specific simulation configurations:

```bash
# Simulation without chemical species
./configure --newdd-nchem=0 --newdd-ndust=0

# Smaller memory pool
./configure --newdd-nmeg=10000

# Custom WGROUPSIZE for I/O tuning
./configure --newdd-wgroupsize=10
```

---

## Source File Descriptions

### Core Files

| File | Lines | Description |
|------|-------|-------------|
| `newdd.c` | ~280 | Main program. MPI initialization, workflow coordination |
| `ramses.h` | ~330 | Core data structures and macro definitions |
| `Memory.h` | ~60 | Memory manager interface |

### Data Reading Modules

| File | Lines | Description |
|------|-------|-------------|
| `rd_info.c` | ~220 | Info file parsing, unit conversion factor computation |
| `rd_amr.c` | ~330 | AMR mesh structure reading, cell hierarchy construction |
| `rd_part.c` | ~110 | Particle data reading, DM/star filtering |
| `rd_hydro.c` | ~180 | Hydrodynamics data reading |
| `rd_sink.c` | ~85 | Sink (black hole) data reading |
| `rd_grav.c` | ~120 | Gravity/potential data reading (optional) |

### Processing Modules

| File | Lines | Description |
|------|-------|-------------|
| `find_leaf_gas.c` | ~250 | Leaf gas cell extraction algorithm |
| `utils.c` | ~450 | Sorting functions, SplitDump (slab decomposition) |
| `header.c` | ~105 | Header file I/O |
| `Memory2.c` | ~80 | Custom memory allocator (stack-based) |

---

## Performance and Memory

### Memory Requirements

| Scale | NMEG Setting | Description |
|-------|--------------|-------------|
| Small | 5000 | ~5 GB |
| Medium | 10000 | ~10 GB |
| Large | 20000 | ~20 GB (default) |

Approximate memory usage per CPU domain:
- AMR mesh: ~40 MB
- Hydrodynamics: ~800 MB
- Particles: ~200 MB

### Performance Characteristics

- **I/O Pattern**: Sequential processing of each CPU domain, parallel writing per slab
- **Sorting Complexity**: O(n log n) (per data type)
- **Slab Decomposition**: O(n) binning
- **Scalability**: Recommended to use MPI process counts that evenly divide ncpu

### Optimization Tips

1. **MPI Process Count**: Use a number that evenly divides ncpu
2. **WGROUPSIZE**: Adjust based on I/O bottlenecks (default: 5)
3. **NMEG**: Adjust according to system memory
4. **nsplit**: Set to match parallelism of downstream analysis tools

---

## Examples

### Basic Execution

```bash
# Decompose snapshot 5 into 64 slabs
./newdd.exe 5 64
```

### MPI Parallel Execution

```bash
# Decompose snapshot 100 into 128 slabs using 8 processes
mpirun -np 8 ./newdd.exe 100 128
```

### SLURM Batch Script

```bash
#!/bin/bash
#SBATCH -J newdd
#SBATCH -N 2
#SBATCH -n 32
#SBATCH -t 02:00:00
#SBATCH -p normal

module load mpi/openmpi

srun ./newdd.exe 50 256
```

### Reading Output Data (C)

```c
#include "ramses.h"

int main() {
    FILE *fp;
    int ndm;
    DmType *dm;

    // Read particle count from metadata
    // (requires parsing info file)
    ndm = 1000000;  // example

    // Read binary data
    fp = fopen("FoF_Data/NewDD.00005/SN.00005.DM.00000.dat", "rb");
    dm = malloc(sizeof(DmType) * ndm);
    fread(dm, sizeof(DmType), ndm, fp);
    fclose(fp);

    // Use data
    for (int i = 0; i < ndm; i++) {
        printf("DM %d: pos=(%.3f, %.3f, %.3f) mass=%.3e\n",
               i, dm[i].x, dm[i].y, dm[i].z, dm[i].mass);
    }

    free(dm);
    return 0;
}
```

### Reading Output Data (Python)

```python
import numpy as np

# Define data type (must match ramses.h)
# Adjust NCHEM, NDUST values to match your compilation flags
NCHEM = 9
dm_dtype = np.dtype([
    ('x', 'f8'), ('y', 'f8'), ('z', 'f8'),
    ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'),
    ('mass', 'f8'),
    ('id', 'i8'),
    ('levelp', 'i4'),
    ('family', 'i1'), ('tag', 'i1'),
    # Fields below are present when compiled without -DNBODY
    ('tp', 'f8'), ('zp', 'f8'),
    # If NCHEM=9
    ('chem', 'f8', NCHEM),
    ('mass0', 'f8'),
    ('birth_d', 'f8'),
    ('partp', 'i4'),
])

# Read binary file
dm = np.fromfile('FoF_Data/NewDD.00005/SN.00005.DM.00000.dat',
                 dtype=dm_dtype)

print(f"Number of dark matter particles: {len(dm)}")
print(f"Mass range: {dm['mass'].min():.3e} - {dm['mass'].max():.3e} Msun/h")
```

---

## Physical Constants

Physical constants used in the code (ramses.h):

```c
#define Mpc        3.085678e24    // 1 Mpc [cm]
#define Msun       1.98892e33     // Solar mass [g]
#define rhoc       1.8791e-29     // Critical density [g/cm^3] (h=1)
#define mH         1.6600000e-24  // Hydrogen atom mass [g]
#define kB         1.3806200e-16  // Boltzmann constant [erg/K]
#define X          0.76           // Hydrogen mass fraction
#define Y          0.24           // Helium mass fraction
```

---

## Troubleshooting

### Common Errors

**1. File not found**
```
Error: Cannot open info file
```
→ Check input path. Verify that `output_XXXXX/` directory exists in the correct location.

**2. Out of memory**
```
Error: Memory allocation failed
```
→ Reduce `NMEG` value in Makefile or check system memory.

**3. MPI error**
```
MPI_Bcast: Invalid communicator
```
→ Verify that the number of MPI processes is a divisor of ncpu.

**4. Segmentation fault**
```
Segmentation fault (core dumped)
```
→ Check data type sizes. Verify that `NCHEM`, `NDUST` settings match the simulation.

### Debugging Tips

1. Test with single processor first: `./newdd.exe 1 4`
2. Start with small nsplit values: 4, 8, 16
3. Check output directory permissions
4. Analyze core dumps with `gdb`

---

## Related Projects

This package is part of the PGalF (Parallel Galaxy Finder) pipeline:

- **NewDD**: Domain decomposition (this package)
- **opFoF**: OpenMP parallel Friends-of-Friends halo finder
- **NewGalFinder**: Galaxy identification within halos
- **GalCenter**: Precise galaxy center finding

For detailed pipeline information, see [PIPELINE.md](PIPELINE.md).

---

## Additional Documentation

- **[QUICKSTART.md](QUICKSTART.md)**: Step-by-step quick start guide for new users
- **[CONFIGURATION.md](CONFIGURATION.md)**: Advanced configuration and performance tuning
- **[PIPELINE.md](PIPELINE.md)**: Complete PGalF pipeline integration guide
- **[README4GADGET.md](README4GADGET.md)**: Guide for adapting NewDD to read GADGET HDF5 format

---

## License

This code was developed for research purposes.

---

## Contact

For questions, please contact the code maintainer.
