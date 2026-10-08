#!/usr/bin/env python3
"""Galaxy Stellar Mass Function (GSMF) from PGalF catalogs, with calibration.

GSMF = galaxies per comoving volume per dex of stellar mass. Plotted in
sim-native units: M* in Msun/h, phi in (Mpc/h)^-3 dex^-1.

Calibration (choose with --calib):
  baldry2012   : Baldry+2012 GAMA double-Schechter, z~0
  davidzon2017 : Davidzon+2017 COSMOS2015, z-dependent (0.2<z<5.5), auto-picks
                 the redshift bin of the snapshot
  none
Plus optional --obs-file (overlaid as points):
    columns  logM[Msun]  logPhi[Mpc^-3 dex^-1]  [err_dex]

Published fits are in h=0.7 (Msun, Mpc^-3); converted to the sim h internally.
Run from the directory with read_pgalf_catalog.py.

    python plot_gsmf.py --root /gpfs/.../pgalf/eagle_spin_jet --snapshot 42 \
        --calib davidzon2017 --output gsmf_00042.png
"""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_pgalf_catalog import read_fof_catalog, read_galaxy_catalog


# Davidzon+2017 COSMOS2015, TOTAL sample (Table 1, h=0.7).
# (zlo, zhi, logMstar, a1, phi1[1e-3], a2, phi2[1e-3]); phi2=None -> single Schechter.
DAVIDZON2017 = [
    (0.2, 0.5, 10.78, -1.38, 1.187, -0.43, 1.92),
    (0.5, 0.8, 10.77, -1.36, 1.070,  0.03, 1.68),
    (0.8, 1.1, 10.56, -1.31, 1.428,  0.51, 2.19),
    (1.1, 1.5, 10.62, -1.28, 1.069,  0.29, 1.21),
    (1.5, 2.0, 10.51, -1.28, 0.969,  0.82, 0.64),
    (2.0, 2.5, 10.60, -1.57, 0.295,  0.07, 0.45),
    (2.5, 3.0, 10.59, -1.67, 0.228, -0.08, 0.21),
    (3.0, 3.5, 10.83, -1.76, 0.090,  None, None),
    (3.5, 4.5, 11.10, -1.98, 0.016,  None, None),
    (4.5, 5.5, 11.30, -2.11, 0.003,  None, None),
]


def _h(hubble):
    return hubble / 100.0 if (hubble and hubble > 1.0) else hubble


def _schechter(logM_msun, logMstar, phi1, a1, phi2=None, a2=None):
    """(Double-)Schechter number density [Mpc^-3 dex^-1] at logM (Msun)."""
    x = 10.0 ** (np.asarray(logM_msun, float) - logMstar)
    term = phi1 * x ** (a1 + 1)
    if phi2 is not None:
        term = term + phi2 * x ** (a2 + 1)
    return np.log(10.0) * np.exp(-x) * term


def baldry2012(logM_msun):
    """Baldry+2012 GAMA z<0.06 GSMF [Mpc^-3 dex^-1], masses Msun (h=0.7)."""
    return _schechter(logM_msun, 10.66, 3.96e-3, -0.35, 0.79e-3, -1.47)


def davidzon2017(logM_msun, redshift):
    """Davidzon+2017 total GSMF [Mpc^-3 dex^-1] at the bin containing `redshift`.
    Returns (phi, label) or (None, None) if z is outside 0.2-5.5."""
    for zlo, zhi, lMs, a1, p1, a2, p2 in DAVIDZON2017:
        if zlo <= redshift < zhi:
            p2v = p2 * 1e-3 if p2 is not None else None
            phi = _schechter(logM_msun, lMs, p1 * 1e-3, a1, p2v, a2)
            return phi, f"Davidzon+2017 ({zlo}<z<{zhi})"
    return None, None


def compute_gsmf(root, snapshot, mass_lo=7.0, mass_hi=12.0, dlog=0.25, min_count=1):
    """Sim GSMF (all galaxies with mstar>0). Returns
    (logM_centers[Msun/h], logPhi[(Mpc/h)^-3 dex^-1], err_dex, header)."""
    step = f"{snapshot:05d}"
    fof_dir = Path(root) / "FoF_Data" / f"FoF.{step}"
    header, _ = read_fof_catalog(fof_dir / f"FoF_halo_cat.{step}")
    _, subs, _ = read_galaxy_catalog(fof_dir / f"GALCATALOG.LIST.{step}")

    mstar = subs["mstar"][subs["mstar"] > 0]          # Msun/h
    volume = header["box_size"] ** 3                  # (Mpc/h)^3
    edges = np.arange(mass_lo, mass_hi + 1e-9, dlog)
    counts, _ = np.histogram(np.log10(mstar), bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    good = counts >= min_count
    phi = counts / (volume * dlog)                    # (Mpc/h)^-3 dex^-1
    err_dex = 0.434 / np.sqrt(np.where(counts > 0, counts, np.nan))
    return centers[good], np.log10(phi[good]), err_dex[good], header


def overlay_calib(ax, xrange, h, z, which):
    """Overplot a Schechter calibration, converted to Msun/h & (Mpc/h)^-3."""
    if which == "none":
        return
    grid_h = np.linspace(xrange[0], xrange[1], 200)
    logM_msun = grid_h - np.log10(h)                  # Msun/h -> Msun
    if which == "baldry2012":
        phi_mpc = baldry2012(logM_msun)
        label = "Baldry+2012 (GAMA, z~0)"
        if abs(z) > 0.5:
            label += f"  [!! snapshot z={z:.2g}]"
    elif which == "davidzon2017":
        phi_mpc, label = davidzon2017(logM_msun, z)
        if phi_mpc is None:
            print(f"davidzon2017: z={z:.3g} outside 0.2-5.5, no curve"); return
    else:
        raise ValueError(which)
    logphi_h = np.log10(phi_mpc) - 3.0 * np.log10(h)  # Mpc^-3 -> (Mpc/h)^-3
    ax.plot(grid_h, logphi_h, color="crimson", lw=2, zorder=3, label=label)


def overlay_obs_file(ax, path, h, label=None):
    """Obs GSMF file: logM[Msun] logPhi[Mpc^-3/dex] [err_dex]."""
    d = np.loadtxt(path)
    x = d[:, 0] + np.log10(h)                         # Msun -> Msun/h
    y = d[:, 1] - 3.0 * np.log10(h)                   # Mpc^-3 -> (Mpc/h)^-3
    err = d[:, 2] if d.shape[1] > 2 else None
    ax.errorbar(x, y, yerr=err, fmt="s", color="navy", ms=4, capsize=2,
                zorder=4, label=label or f"obs: {Path(path).name}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--snapshot", type=int, required=True)
    p.add_argument("--dlog", type=float, default=0.25, help="mass bin width [dex]")
    p.add_argument("--calib", choices=("baldry2012", "davidzon2017", "none"),
                   default="davidzon2017")
    p.add_argument("--obs-file", type=Path, default=None,
                   help="logM[Msun] logPhi[Mpc^-3/dex] [err_dex]")
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    centers, logphi, err, header = compute_gsmf(args.root, args.snapshot, dlog=args.dlog)
    if len(centers) == 0:
        raise SystemExit("No galaxies with mstar>0 found")
    z = header["redshift"]; h = _h(header["hubble"])

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.errorbar(centers, logphi, yerr=err, fmt="o", color="C0", ms=5,
                capsize=2, label=f"PGalF {Path(args.root).name}")
    overlay_calib(ax, (centers.min(), centers.max()), h, z, args.calib)
    if args.obs_file:
        overlay_obs_file(ax, args.obs_file, h)

    ax.set_xlabel(r"$\log_{10}(M_\star/[M_\odot/h])$")
    ax.set_ylabel(r"$\log_{10}(\phi/[(h/{\rm Mpc})^3\,{\rm dex}^{-1}])$")
    ax.set_title(f"GSMF  snapshot {args.snapshot:04d}, z={z:.3g}")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    fig.tight_layout()
    out = args.output or Path(f"gsmf_{args.snapshot:05d}.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160); plt.close(fig)
    print(f"z={z:.4g}; calib={args.calib}; wrote {out.resolve()}")


if __name__ == "__main__":
    main()
