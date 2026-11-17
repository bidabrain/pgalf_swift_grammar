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
| gfind stellar density vs `PEAKTHRESHOLD` (h² Msun/ckpc³) | density lower by h² → threshold effectively **≈2.1× higher**; TSC cell / smoothing physically 1/h larger |
| Unbinding (G m / r) | m and r off by the same h → roughly unaffected |

**Re-run from NewDD** (all four stages) after rebuilding NewDD. Write to a new
output directory (or remove the old one): `run_pgalf.sh -r` would treat the old
outputs as finished and skip them. Calibration runs processed with the old
NewDD must be re-run too.
