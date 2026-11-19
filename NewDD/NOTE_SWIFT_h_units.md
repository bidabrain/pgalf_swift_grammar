# Note: SWIFT snapshots and little-h units

## TL;DR

- **SWIFT snapshots carry no h.** Coordinates and `BoxSize` are comoving Mpc,
  masses are 10^10 Msun (see `/InternalCodeUnits`). h is only stored as a
  cosmological parameter in **`/Cosmology/h`**; there is **no
  `/Header/HubbleParam`** (that is a GADGET attribute).
- **PGalF works in h units.** NewDD multiplies lengths and masses by h when it
  reads the snapshot (`mpcscale_l`, `scale_m` in `rd_gadget.c`), so every PGalF
  output is in **cMpc/h** and **Msun/h** (velocities stay km/s). The opFoF and
  gfind parameters (linking length, `PEAKTHRESHOLD`, `MINSTELLARMASS`, ...) and
  all Python analysis scripts assume these units.
- **Bug (fixed in this commit):** NewDD only looked for `/Header/HubbleParam`.
  On SWIFT files it was missing, so h silently stayed at its default **1.0** and
  *no h conversion was applied*. Everything downstream was really in cMpc and
  Msun while being treated as cMpc/h and Msun/h.
- **Fix:** `rd_gadget_info()` now reads h (and Ωm, ΩΛ) in this order:
  1. `/Header/HubbleParam`, `Omega0`, `OmegaLambda` (GADGET)
  2. `/Cosmology/h`, `Omega_m`, `Omega_lambda` (SWIFT)
  3. otherwise h = 1 **with an error message** that lengths/masses will not be
     in h units.

  It prints, for the first file of each snapshot:
  `rd_gadget_info: h=0.6933 Omega_m=... BoxSize=25 -> boxlen=17.3324 cMpc/h`.

## How to tell whether an output is affected

```bash
h5dump -a /Header/HubbleParam snap.hdf5   # "unable to open attribute" on SWIFT
h5dump -a /Cosmology/h        snap.hdf5   # e.g. 0.6933
```
```python
from read_pgalf_catalog import read_fof_catalog
hdr, _ = read_fof_catalog(".../FoF.NNNNN/FoF_halo_cat.NNNNN")
print(hdr["box_size"])   # affected: = BoxSize in Mpc (e.g. 25.0)
                         # fixed:    = BoxSize * h     (e.g. 17.33)
```

The `hubble` field of the FoF header is **not** a good indicator: opFoF and
gfind overwrite the cosmology with the compile-time `REQUIRED_*` values
(H0 = 69.33), but that override does not touch box size, positions or masses.

## Consequences for outputs made before the fix (SWIFT input)

Actual units were cMpc and Msun (h = 0.6933 in the example run):

| Where | Effect |
|---|---|
| GSMF volume (`box_size`³) | over-estimated by h⁻³ ≈ 3.0 → φ **too low by ≈0.48 dex** |
| M* read as Msun/h | too high by 1/h → GSMF shifted **+0.16 dex in M*** |
| Aperture radius (pkpc → cMpc/h uses ×h, coords in cMpc) | "50 pkpc" was really **≈35 pkpc** (30 → 21, 100 → 69) |
| opFoF linking length `0.2·(m/ρ̄)^{1/3}` | effective **b ≈ 0.2·h^{2/3} ≈ 0.157** instead of 0.2 |
| gfind stellar density vs `PEAKTHRESHOLD` (h² Msun/ckpc³) | density lower by h² → threshold effectively **≈2.1× higher** |
| TSC cell / Gaussian smoothing (`TSC_CELL_SIZE` 0.004, `Gaussian_Smoothing_Length` 0.008, code length) | 4 / 8 ckpc instead of 4 / 8 ckpc/h → physically **smaller by h** (finer grid, less smoothing) |
| `MINSTELLARMASS` 2e6 (TYPE_ALL/TYPE_STAR switch) | applied to Msun instead of Msun/h → switch at ≈1.4e6 Msun/h (**lower by h**); the 32-star-particle criterion is unaffected |
| Unbinding (G m / r) | m and r off by the same h → roughly unaffected |

### What can and cannot be corrected afterwards

- **Unit mislabelling (first three rows)** is a fixed power of h for every
  object: box/positions ×h, masses ×h, φ ×h⁻³ with the corrected volume, and
  apertures can be recomputed from the particles in `GALFIND.DATA` with the
  right radius. These could in principle be fixed in post-processing.
- **The finder parameters (FoF linking length, `PEAKTHRESHOLD`, TSC grid and
  smoothing, `MINSTELLARMASS`)** changed *which particles form a galaxy*. FoF
  groups are smaller and more fragmented (host masses and total M* biased
  low), and fewer density peaks pass the threshold (satellites missed or
  merged into centrals, biasing the low-mass end). **This cannot be undone by
  rescaling**; only a re-run fixes it.
- Unit-corrected pre-fix catalogues are roughly usable for qualitative looks
  at the high-mass end (M* ≳ 10^10 Msun/h), where stellar peaks are strong. The
  low-mass end and satellite counts carry a systematic bias of unknown size
  (part of the low-mass GSMF shape, e.g. the step near 10^8, may come from
  it). **Never use them for emulator calibration**;
  `emulator_grammar/measure_gsmf.py` refuses them (`check_little_h`).

### Recommendation

**Re-run from NewDD** (all four stages) after rebuilding NewDD
(`cd NewDD && make clean && make all`). Check that the NewDD log prints
`h=0.6933 ... -> boxlen=17.33 cMpc/h` and that the FoF header `box_size` is
17.33, not 25. Write to a new
output directory (or remove the old one): `run_pgalf.sh -r` would treat the old
outputs as finished and skip them. Calibration runs processed with the old
NewDD must be re-run too.
