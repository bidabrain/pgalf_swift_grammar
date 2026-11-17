# SWIFT/GADGET HDF5 Mode for NewDD

This document describes how to run `NewDD` on SWIFT/GADGET-style HDF5 snapshots for halo finding.

## What Was Added

- Optional build mode with `-DGADGET_HDF5`
- HDF5 snapshot readers in `rd_gadget.c` for:
  - DM + stars (`PartType1/2/3/4`)
  - gas (`PartType0`)
  - black holes/sinks (`PartType5`)
- Runtime support in `newdd.c` for GADGET/SWIFT snapshot loops

## Units and little h

SWIFT stores lengths and masses without h and keeps h in `/Cosmology/h`
(there is no `/Header/HubbleParam`). `rd_gadget_info()` reads h from
`/Header/HubbleParam`, then `/Cosmology/h`, and converts positions/box to
cMpc/h and masses to Msun/h; if neither exists it warns and uses h=1. Earlier
versions skipped the `/Cosmology` fallback, so SWIFT outputs made with them are
in cMpc / Msun. Details and consequences: [NOTE_SWIFT_h_units.md](NOTE_SWIFT_h_units.md).

## Build

From `GalaxyFinder/NewDD`:

```bash
make clean
make all FORMAT_OPT=-DGADGET_HDF5 HDF5_DIR=/path/to/hdf5
```

Notes:
- `HDF5_DIR` must contain:
  - headers at `${HDF5_DIR}/include/hdf5.h`
  - library at `${HDF5_DIR}/lib/libhdf5.*`
- If your environment module already sets include/lib paths, you can keep `HDF5_DIR` default and pass extra flags through `OPT`/`LDFLAGS` if needed.

## Input Snapshot Paths

For snapshot step `N`, `newdd.exe` tries:

1. `./snapdir_NNN/snapshot_NNN(.file).hdf5`
2. fallback: `./snapshot_NNN(.file).hdf5`

where `NNN` is zero-padded to 3 digits.

Examples:

- single-file: `snapshot_010.hdf5`
- multi-file: `snapshot_010.0.hdf5`, `snapshot_010.1.hdf5`, ...

## Run

```bash
./newdd.exe <ISTEP> <NSPLIT>
```

Example:

```bash
./newdd.exe 10 64
```

- `ISTEP=10` looks for `snapshot_010*`
- `NSPLIT=64` writes 64 x-slabs

## Output

Outputs stay compatible with the existing NewDD/FoF pipeline and are written under:

`FoF_Data/NewDD.XXXXX/`

with `.DM`, `.STAR`, `.GAS`, and `.SINK` slab files.

## Mapping Notes

- Particles:
  - `PartType1/2/3` -> `family=1` (DM)
  - `PartType4` -> `family=2` (STAR)
- Gas:
  - `SmoothingLength` -> `cellsize`
  - `Density` -> `den`
  - `InternalEnergy` -> temperature-like field via `scale_T2`
- BH:
  - `BH_Mass` (or `Masses`) -> sink mass
  - `BH_Mdot` -> `dMsmbh`

## Current Assumptions

- Snapshot header has standard GADGET attributes (`Header` group).
- Velocity conversion follows GADGET/SWIFT convention using `sqrt(a)`.
- Chemical/dust extended fields are not fully populated in GADGET path yet (base fields are).

## Quick Validation Checklist

1. Confirm HDF5 link:
   - `ldd newdd.exe | grep hdf5` (Linux)
   - `otool -L newdd.exe | grep hdf5` (macOS)
2. Run one small snapshot and check:
   - `FoF_Data/NewDD.XXXXX/` exists
   - DM/STAR/GAS files are non-empty when expected
3. Run downstream FoF and compare counts against a reference.
