# PGalF Suite: Galaxy Finding Pipeline for RAMSES Simulations

A comprehensive suite for identifying galaxies and halos from RAMSES cosmological simulation outputs. The pipeline processes raw simulation data through domain decomposition, halo finding, and galaxy identification.

## Pipeline Overview

```
RAMSES Output → NewDD → opFoF → NewGalFinder → Galaxy Catalogs
               (Slab)   (Halo)   (Galaxy)
```

## Components

| Directory | Description |
|-----------|-------------|
| **NewDD** | RAMSES domain decomposition into Z-directional slabs for parallel processing |
| **opFoF** | MPI-parallel Friends-of-Friends (FoF) halo finder |
| **NewGalFinder** | Parallel galaxy/subhalo finder using density peaks and water-shedding |
| **GalCenter** | Galaxy center identification utilities |

---

## Building & Running on the cluster (GCC 14 / OpenMPI 5, SWIFT HDF5)

This is the **actual working setup** for this fork on `grammar` (OpenHPC,
GNU 14.2 + OpenMPI 5). The committed Makefiles are already configured for this
toolchain and this cluster's library paths — on the same cluster you usually do
**not** need to run `./configure`; just load the modules and `make`.

### 1. Environment (modules + runtime libs)

```bash
module purge
module load gnu14/14.2.0 openmpi5/5.0.7 hwloc/2.12.0 libfabric/1.18.0 gsl/2.8 hdf5/1.14.6
module unload ucx
module load ucx/1.18.0-mt
export OMP_NUM_THREADS=1                                  # REQUIRED at runtime (see notes)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/fftw3/lib
```

### 2. Build (order matters; compile on the login node)

```bash
cd NewGalFinder && make clean && make all && cd ..       # -> gfind.exe
cd GalCenter    && make new                   && cd ..   # -> galcenter.exe
cd NewDD        && make all                   && cd ..   # -> newdd.exe (+ NewDD/libmyram.a)
cp NewDD/libmyram.a ./                                   # opFoF links ../libmyram.a
cd opFoF        && make clean && make this    && cd ..   # -> opfof.exe
```

> **Do NOT run `./configure` for NewDD** — the HDF5 support (`-DGADGET_HDF5`,
> `rd_gadget.o`, HDF5 include/lib paths) lives only in the committed
> `NewDD/Makefile`; re-running configure regenerates it from `Makefile.in` and
> wipes HDF5. If you re-run configure for the others, re-add the GCC flags below.

### 3. Toolchain porting notes (Intel → GCC/OpenMPI)

Compilers were switched from Intel (`mpiicx`/`mpiifx`) to `mpicc`/`mpif90`, and
GCC 14 needs these "relaxing" flags (already in the committed Makefiles):

| Language | Flags | Why |
|----------|-------|-----|
| C (`mpicc`) | `-fcommon -fpermissive` | GCC 14 promotes implicit-decl / int-conversion to errors |
| Fortran (`mpif90`) | `-std=legacy -fallow-argument-mismatch` | legacy arg-mismatch is an error in gfortran 10+ |
| OpenMP | `-fopenmp` (not Intel `-qopenmp`) | |

A GCC-incompatible `aa = (void*)star = malloc(...)` (cast-as-lvalue) in
`opFoF/ramses2read.c` was changed to `aa = star = malloc(...)`.

### 4. Key settings — where they live (edit for a different cluster/run)

| Setting | Value (this run) | File(s) |
|---------|------------------|---------|
| HDF5 include/lib | `/opt/ohpc/pub/libs/gnu14/hdf5/1.14.6/{include,lib}` | `NewDD/Makefile` (`HDF5_INC`/`HDF5_LIB`) |
| FFTW (float+omp) | `/home/dbi224/local/omp5/fftw3` | `NewGalFinder/Makefile`, `configure` |
| Cosmology (Ωm, Ωb, H0, ΩΛ) | 0.287845, 0.047143, 69.33, 0.712155 | `REQUIRED_*` in `opFoF/opfof.c`, `GalCenter/galcenter.c`, **`NewGalFinder/gfind.c`** (all three!) |
| `NCHEM` (particle struct size, must match across stages) | 9 | `NewDD/Makefile`, `opFoF/Rules.make`, `GalCenter/Makefile`, `NewGalFinder/Makefile` |
| `NMEG` (per-rank memory pool, MB) | NewDD 20000, opFoF 17000, gfind 12000 | respective Makefiles / `opFoF/Rules.make` |

> Cosmology is **compile-time** (`#define REQUIRED_OMEGA_M ...`). All three C
> programs override the header with these — change in **all three** and rebuild.

### 5. Running the pipeline (`run_pgalf.sh`)

SWIFT snapshots (e.g. `flamingo_0042.hdf5`) are symlinked to the name NewDD
expects (`snap_0042.hdf5`) automatically by the driver script:

```bash
./run_pgalf.sh -i <sim_dir> -o <out_dir> -s <snap|a-b|all> \
               -n <nsplit> -N <gfind_ranks> -z <zmax> -m "srun -n" [-r] [-x prefix]
```

- `-n` = NewDD slabs **and** opFoF ranks (opFoF requires `np == nsplit`).
- `-N` = gfind ranks (independent). **Keep small (≤8)** — gfind master/slave
  **deadlocks at ~32 ranks** under OpenMPI 5; 4–8 works.
- `-z` skips snapshots with `z > zmax` (FoF percolates / is meaningless at very
  high z). `-r` resumes (skips snapshots whose final catalog already exists).
- `-m "srun -n"` on SLURM; `-m "mpirun -np"` on a login node.
- `export OMP_NUM_THREADS=1` first, always — the binaries are `-fopenmp` +
  `fftw3f_omp`; unset, each rank spawns a full thread pool and exhausts the node
  (OOM / `libgomp: Thread creation failed`).

Output per snapshot `NNNNN` → `<out_dir>/FoF_Data/FoF.NNNNN/`:
`FoF_halo_cat.NNNNN` (halos), `GALCATALOG.LIST.NNNNN` (galaxies/subhalos),
`GALFIND.CENTER.NNNNN` (centers). See `read_pgalf_catalog.py` for the formats.

### 6. SLURM batch

Request `--ntasks=<nsplit>` (NOT `--cpus-per-task=<n>`): opFoF is MPI, not
OpenMP. Example header: `-N 1 --ntasks=32 --cpus-per-task=1 --exclusive --mem=0`,
then `export OMP_NUM_THREADS=1` and call `run_pgalf.sh ... -m "srun -n"`. For
many snapshots use a job array (`#SBATCH --array=a-b%N`, one snapshot per task
via `-s $SLURM_ARRAY_TASK_ID`).

---

## Quick Start

### 0. Configure Build System
```bash
# Default configuration (Intel OneAPI mpiicx/mpiifx compilers, optimized build)
./configure

# Debug build
./configure --debug

# Custom compilers and FFTW path
./configure --cc=mpicc --fc=mpifort --fftw=/path/to/fftw

# See all options
./configure --help
```

### 1. Build All Components
```bash
make all          # Build everything
# Or build individually:
make galcenter    # Build GalCenter only
make galfinder    # Build NewGalFinder only
make newdd        # Build NewDD only
```

### 2. Domain Decomposition
```bash
mpirun -np 8 ./NewDD/newdd.exe <snapshot> <nsplit>
```

### 3. Halo Finding
```bash
cd opFoF
make this
mpirun -np 8 ./opfof.exe <snapshot> <nfiles>
```

### 4. Galaxy Finding
```bash
mpirun -np 64 ./NewGalFinder/gfind.exe <snapshot>
```

## Build System

The project uses a `configure` script to generate Makefiles for all subdirectories from `Makefile.in` templates.

| File | Description |
|------|-------------|
| `configure` | Build configuration script |
| `Makefile.in` | Top-level Makefile template |
| `GalCenter/Makefile.in` | GalCenter Makefile template |
| `NewGalFinder/Makefile.in` | NewGalFinder Makefile template |
| `NewDD/Makefile.in` | NewDD Makefile template |

### Configure Options

| Option | Default | Description |
|--------|---------|-------------|
| `--cc=CC` | `mpiicx` | C compiler for GalCenter/NewDD |
| `--fc=FC` | `mpiifx` | Fortran compiler for GalCenter/NewDD |
| `--galfinder-cc=CC` | `mpiicx` | C compiler for NewGalFinder |
| `--galfinder-fc=FC` | `mpiifx` | Fortran compiler for NewGalFinder |
| `--fftw=PATH` | `/home/kjhan/local` | FFTW installation path |
| `--opt=FLAGS` | `-O3` | Optimization flags |
| `--debug` | - | Use `-g` debug flags |
| `--openmp=FLAGS` | `-qopenmp` | OpenMP flags |
| `--no-openmp` | - | Disable OpenMP |

Per-component options (e.g., `--galcenter-nmeg=N`, `--galfinder-nchem=N`, `--newdd-ndust=N`) are also available. Run `./configure --help` for the full list.

## Shared Files

- `ramses.h` - RAMSES data structure definitions
- `params.h` - Algorithm parameters
- `libmyram.a` - Custom memory management library

## Documentation

Each subdirectory contains detailed documentation:
- [NewDD/README.md](NewDD/README.md) - Domain decomposition details
- [opFoF/README.md](opFoF/README.md) - FoF algorithm and usage
- [NewGalFinder/README.md](NewGalFinder/README.md) - Galaxy finder algorithm

## Requirements

- MPI (OpenMPI, MPICH, or Intel MPI)
- Intel or GCC compilers with OpenMP support
- FFTW3 library (for NewGalFinder)

## Output Units

| Quantity | Unit |
|----------|------|
| Position | comoving Mpc/h |
| Velocity | km/s |
| Mass | Msun/h |

## License

For scientific research purposes.
