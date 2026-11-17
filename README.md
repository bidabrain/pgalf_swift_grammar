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
| `NCHEM`, `NDUST`, `NPRE` (particle struct size — **must match across ALL stages**, incl. NewGalFinder) | 9, 4, 8 | `NewDD/Makefile`, `opFoF/Rules.make`, `GalCenter/Makefile`, `NewGalFinder/Makefile` |
| `NMEG` (per-rank memory pool, MB) | NewDD 20000, opFoF 17000, gfind 12000 | respective Makefiles / `opFoF/Rules.make` |

> Cosmology is **compile-time** (`#define REQUIRED_OMEGA_M ...`). All three C
> programs override the header with these — change in **all three** and rebuild.

> **SWIFT and little h:** SWIFT snapshots store lengths/masses *without* h and
> keep h in `/Cosmology/h` (no `/Header/HubbleParam`). NewDD converts to
> cMpc/h and Msun/h using that h; before the fix in
> [`NewDD/NOTE_SWIFT_h_units.md`](NewDD/NOTE_SWIFT_h_units.md) it fell back to
> h=1 on SWIFT input, so older outputs are really in cMpc / Msun (e.g. FoF
> `box_size` = 25 instead of 17.33). The `REQUIRED_H0` override above does
> **not** affect this conversion. Check the NewDD log for
> `rd_gadget_info: h=... -> boxlen=... cMpc/h`.

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
`GALFIND.DATA.NNNNN` (per-subhalo member particles), `GALFIND.CENTER.NNNNN`
(centers). See `read_pgalf_catalog.py` for the formats.

Size cuts with the default build: opFoF keeps FoF groups with ≥20 particles
(all types, `MinNumMem` in `opFoF/fof.h`), and gfind only processes hosts with
≥30 (`gfind.c`), so every `GALCATALOG.LIST` host has ≥30 particles. There is
**no** minimum particle count for a subhalo/galaxy: the `MINCORENMEM` cut in
`NewGalFinder/params.h` sits under `#ifdef OLD` and is not compiled. Apply your
own cut (e.g. on `npstar`) when analysing low-mass galaxies.

gfind also switches peak finder per FoF halo (`subhaloden.mod6.c`, `subhalo_den`):
halos with ≤32 star particles (`NUMNEIGHBOR`) or M*<`MINSTELLARMASS` use an
unweighted all-particle SPH number density (DM-dominated; one peak ⇒ the whole
FoF halo becomes one galaxy, no unbinding; no peak above `PEAKTHRESHOLD` ⇒ the
halo is dropped), the rest use a mass-weighted stellar TSC density with
Gaussian smoothing (`-DADV`). The same `PEAKTHRESHOLD` is compared against
different quantities in the two branches. Expect a feature in the GSMF near
M* ≈ 32 × (star particle mass).

### 6. SLURM batch

Request `--ntasks=<nsplit>` (NOT `--cpus-per-task=<n>`): opFoF is MPI, not
OpenMP. Example header: `-N 1 --ntasks=32 --cpus-per-task=1 --exclusive --mem=0`,
then `export OMP_NUM_THREADS=1` and call `run_pgalf.sh ... -m "srun -n"`. For
many snapshots use a job array (`#SBATCH --array=a-b%N`, one snapshot per task
via `-s $SLURM_ARRAY_TASK_ID`).

### 7. Analysis plots

Run from the repo root (they import `read_pgalf_catalog.py`). `--root` is the
`run_pgalf.sh` `<out_dir>`.

> Alternative finder: [HBT-HERONS-GRAMMAR](https://github.com/bidabrain/HBT-HERONS-GRAMMAR)
> (`grammar/run_hbt.sh`, `grammar/plot_gsmf_hbt.py`, ...) produces the same
> plots in the same units from HBT-HERONS catalogues, reusing the reference
> curves and binning from the scripts below; keep this repo next to it or set
> `PGALF_REPO`.

```bash
# Galaxy stellar mass function vs. a Schechter calibration
python plot_gsmf.py --root <out_dir> --snapshot 274 --dlog 0.2 \
    --calib driver2022 --aperture-kpc 0 30 50 100 --output gsmf.png

# Stellar-to-halo mass relation (Moster+13; Behroozi+19 with --um-dir)
python plot_stellar_halo_relation.py --root <out_dir> --snapshot 274 \
    --population central --aperture-kpc 0 30 50 100 --output shmr.png
```

- `--aperture-kpc` (default `50`): M* inside a 3D proper-kpc aperture around
  the subhalo centre of mass, computed from `GALFIND.DATA`; `0` (or <0) = the
  catalogue's total bound M*. Several values are overplotted: GSMF curves
  per aperture; for the SHMR, light scatter plus a binned median per aperture
  (`--dlog`, default 0.2 dex). `GALFIND.DATA` is read once for all apertures.
- `--center com|galcenter` (both scripts, default `com`): aperture centre.
  `com` is the subhalo centre of mass over all particle types (GALCATALOG
  `x,y,z`), which can sit far from the stellar body for single-peak FoF halos
  (whole FoF kept, no unbinding), mergers or stripped satellites. `galcenter`
  uses the GalCenter stellar centre from `GALFIND.CENTER` (shrinking sphere on
  the star particles for >30 stars, star centre of mass for 1–30), matched by
  `(halo_id, subhalo_id)`; subhalos without an entry or stars fall back to the
  COM. Needs the galcenter stage. Periodic boxes are fine: opFoF (run with
  `<snap> <nsplit>` as in `run_pgalf.sh`, i.e. periodic mode) unwraps every
  FoF halo crossing the box edge before writing its particles, so both
  centres see contiguous positions (possibly <0). GalCenter's own x-wrap in
  `find_centers.c` (which has a `MAX(xmin, …)` typo and assumes boxes >700
  cMpc/h) only matters for non-periodic opFoF runs.
- `--min-npstar N` (both scripts, default 0 = off): keep only galaxies with
  at least N bound star particles (`npstar` in `GALCATALOG.LIST`), for every
  aperture; in the SHMR the cut is applied before the central is chosen. The
  script prints the typical star-particle mass (median `mstar/npstar`) and the
  M* this corresponds to. The cut is on the *total* bound count, so aperture
  curves are not sharply truncated at N × m★. Below ~100 star particles the
  GSMF is resolution-dominated (and gfind switches peak finders at 32 star
  particles per FoF halo, see "Size cuts" above), so compare runs with e.g.
  `--min-npstar 0`, `32`, `100`, `200` before trusting the low-mass end.
- `--calib` (GSMF): `driver2022`, `baldry2012` (z≈0), `davidzon2017`
  (0.2<z<5.5, picks the snapshot's bin), `none`; `--obs-file` overlays points.
- Units are sim-native: M* in Msun/h, φ in (Mpc/h)⁻³ dex⁻¹; published fits are
  converted from h=0.7.
- `GALFIND.DATA` record sizes depend on the gfind build flags (`NCHEM`,
  `NDUST`, …). The reader auto-detects them from the file size and
  `GALCATALOG.LIST`, so no edits are needed after a rebuild. A
  `Truncated subhalo header` error means an old (pre-`3d0b545`) reader with
  hard-coded sizes: `git pull`.

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

These hold for SWIFT input only when NewDD found h (`/Cosmology/h`); outputs
made before that fix are in comoving Mpc and Msun. See
[`NewDD/NOTE_SWIFT_h_units.md`](NewDD/NOTE_SWIFT_h_units.md).

## License

For scientific research purposes.
