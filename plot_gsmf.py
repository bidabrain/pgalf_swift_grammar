#!/usr/bin/env python3
"""Galaxy Stellar Mass Function (GSMF) from PGalF catalogs, with calibration.

GSMF = number of galaxies per comoving volume per dex of stellar mass.
Masses are in Msun/h and number densities in (Mpc/h)^-3 dex^-1 (sim-native),
consistent with the other plotting scripts.

Calibration overlays:
  - built-in Baldry+2012 (GAMA) double-Schechter at z~0 (shown with a warning
    if the snapshot is far from z=0 -- use a z-matched obs set instead);
  - optional --obs-file with columns:  logM[Msun]  logPhi[Mpc^-3 dex^-1]  [err_dex]
    (converted to Msun/h and (Mpc/h)^-3 internally using the sim h).

Run from the directory containing read_pgalf_catalog.py.

    python plot_gsmf.py --root /gpfs/.../pgalf/eagle_spin_jet --snapshot 42 \
        --output gsmf_00042.png [--obs-file davidzon2017_z3.txt]
"""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_pgalf_catalog import read_fof_catalog, read_galaxy_catalog


def _h(hubble):
    return hubble / 100.0 if (hubble and hubble > 1.0) else hubble


def schechter_double(logM_msun, logMstar, phi1, a1, phi2, a2):
    """Double-Schechter number density [Mpc^-3 dex^-1] at logM (Msun)."""
    x = 10.0 ** (np.asarray(logM_msun, float) - logMstar)
    return np.log(10.0) * np.exp(-x) * (phi1 * x ** (a1 + 1) + phi2 * x ** (a2 + 1))


def baldry2012(logM_msun):
    """Baldry+2012 GAMA z<0.06 GSMF [Mpc^-3 dex^-1], masses in Msun (h=0.7)."""
    return schechter_double(logM_msun, 10.66, 3.96e-3, -0.35, 0.79e-3, -1.47)


def compute_gsmf(root, snapshot, mass_lo=7.0, mass_hi=12.0, dlog=0.25, min_count=1):
    """Sim GSMF from the PGalF galaxy catalog.

    Returns (logM_centers[Msun/h], logPhi[(Mpc/h)^-3 dex^-1], err_dex, header).
    Counts ALL galaxies/subhalos with mstar>0 (centrals + satellites).
    """
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
    return (centers[good], np.log10(phi[good]), err_dex[good], header)


def overlay_baldry(ax, logM_hunits_range, h, redshift):
    """Overplot Baldry+2012 converted to Msun/h and (Mpc/h)^-3 dex^-1."""
    grid_h = np.linspace(*logM_hunits_range, 200)
    logM_msun = grid_h - np.log10(h)                  # Msun/h -> Msun
    phi_mpc = baldry2012(logM_msun)
    logphi_h = np.log10(phi_mpc) - 3.0 * np.log10(h)  # Mpc^-3 -> (Mpc/h)^-3
    lab = "Baldry+2012 (GAMA, z~0)"
    if abs(redshift) > 0.5:
        lab += f"  [!! snapshot z={redshift:.2g}, use z-matched obs]"
    ax.plot(grid_h, logphi_h, color="crimson", lw=2, zorder=3, label=lab)


def overlay_obs_file(ax, path, h, label=None):
    """Overplot an obs GSMF file: logM[Msun] logPhi[Mpc^-3/dex] [err_dex]."""
    d = np.loadtxt(path)
    logM_msun, logphi_mpc = d[:, 0], d[:, 1]
    err = d[:, 2] if d.shape[1] > 2 else None
    x = logM_msun + np.log10(h)                       # Msun -> Msun/h
    y = logphi_mpc - 3.0 * np.log10(h)                # Mpc^-3 -> (Mpc/h)^-3
    ax.errorbar(x, y, yerr=err, fmt="s", color="navy", ms=4, capsize=2,
                zorder=4, label=label or f"obs: {Path(path).name}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--snapshot", type=int, required=True)
    p.add_argument("--dlog", type=float, default=0.25, help="mass bin width [dex]")
    p.add_argument("--obs-file", type=Path, default=None,
                   help="obs GSMF: logM[Msun] logPhi[Mpc^-3/dex] [err_dex]")
    p.add_argument("--no-baldry", action="store_true")
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    centers, logphi, err, header = compute_gsmf(args.root, args.snapshot, dlog=args.dlog)
    if len(centers) == 0:
        raise SystemExit("No galaxies with mstar>0 found")
    z = header["redshift"]; h = _h(header["hubble"])

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.errorbar(centers, logphi, yerr=err, fmt="o", color="C0", ms=5,
                capsize=2, label=f"PGalF {Path(args.root).name}")
    if not args.no_baldry:
        overlay_baldry(ax, (centers.min(), centers.max()), h, z)
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
    print(f"galaxies binned; z={z:.4g}; wrote {out.resolve()}")


if __name__ == "__main__":
    main()
