#!/usr/bin/env python3
"""Overplot two (or more) PGalF runs on the SAME stellar-halo and BH-stellar
figures, for one snapshot. Reuses the loaders/calibrations of the single-model
scripts, so run it from the same directory as those .py files.

Usage:
    python plot_compare_models.py \
        --root eagle=/gpfs/.../pgalf/eagle_spin_jet \
        --root flamingo=/gpfs/.../pgalf/flamingo_spin_jet_th_kin25 \
        --snapshot 42 \
        --um-dir /gpfs/.../pgalf/umachine-dr1 \
        --outdir /gpfs/.../pgalf

Produces shmr_compare_<snap>.png and mbh_mstar_compare_<snap>.png.
"""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_stellar_halo_relation import (
    load_relation, moster13, behroozi19_um, _um_param_file,
)
from plot_bh_stellar_relation import load_bh_stellar, reines_volonteri15
from plot_gsmf import compute_gsmf, overlay_baldry


def parse_roots(root_args):
    """Each --root is 'label=path' (or just 'path', label = dir name)."""
    out = []
    for r in root_args:
        if "=" in r:
            label, path = r.split("=", 1)
        else:
            label, path = Path(r).name, r
        out.append((label, Path(path)))
    return out


def _h(hubble):
    return hubble / 100.0 if (hubble and hubble > 1.0) else hubble


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", action="append", required=True,
                   help="label=path (repeatable)")
    p.add_argument("--snapshot", type=int, required=True)
    p.add_argument("--population", choices=("central", "subhalo"), default="central")
    p.add_argument("--um-dir", type=Path, default=None)
    p.add_argument("--outdir", type=Path, default=Path("."))
    args = p.parse_args()

    models = parse_roots(args.root)
    args.outdir.mkdir(parents=True, exist_ok=True)
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    # ------------------- Stellar-halo -------------------
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    hub = z = None
    allx = []
    for i, (label, root) in enumerate(models):
        try:
            header, mh, ms = load_relation(root, args.snapshot, args.population)
        except (FileNotFoundError, ValueError) as e:
            print(f"{label}: skip ({e})"); continue
        if len(mh) == 0:
            print(f"{label}: no galaxies"); continue
        ax.scatter(np.log10(mh), np.log10(ms), s=10, alpha=0.5, linewidths=0,
                   color=colors[i], label=f"{label} (N={len(mh)})")
        hub, z = header["hubble"], header["redshift"]
        allx.append(np.log10(mh))
    h = _h(hub)
    if allx:
        grid = np.linspace(np.concatenate(allx).min(), np.concatenate(allx).max(), 200)
        ax.plot(grid, moster13(grid - np.log10(h), z) + np.log10(h),
                color="crimson", lw=2, zorder=3, label="Moster+2013")
        if args.um_dir:
            pf = _um_param_file(args.um_dir, args.population)
            if pf:
                m_um, sm_um = behroozi19_um(pf, z)
                if m_um is not None:
                    ax.plot(m_um + np.log10(h), sm_um + np.log10(h),
                            color="navy", lw=2, ls="--", zorder=3,
                            label="Behroozi+2019 (UM)")
    ax.set_xlabel(r"$\log_{10}(M_{\rm halo}/[M_\odot/h])$")
    ax.set_ylabel(r"$\log_{10}(M_\star/[M_\odot/h])$")
    ax.set_title(f"Stellar-to-halo, snapshot {args.snapshot:04d}, "
                 f"z={z:.3g} ({args.population})")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    out1 = args.outdir / f"shmr_compare_{args.snapshot:05d}.png"
    fig.tight_layout(); fig.savefig(out1, dpi=160); plt.close(fig)

    # ------------------- BH-stellar -------------------
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    hub = z = None
    allx = []
    for i, (label, root) in enumerate(models):
        try:
            header, ms, mbh = load_bh_stellar(root, args.snapshot, args.population)
        except (FileNotFoundError, ValueError) as e:
            print(f"{label}: skip ({e})"); continue
        if len(ms) == 0:
            print(f"{label}: no BH galaxies"); continue
        ax.scatter(np.log10(ms), np.log10(mbh), s=10, alpha=0.5, linewidths=0,
                   color=colors[i], label=f"{label} (N={len(ms)})")
        hub, z = header["hubble"], header["redshift"]
        allx.append(np.log10(ms))
    h = _h(hub)
    if allx:
        grid = np.linspace(np.concatenate(allx).min(), np.concatenate(allx).max(), 100)
        xmsun = grid - np.log10(h)
        ax.plot(grid, reines_volonteri15(xmsun, "agn") + np.log10(h),
                color="crimson", lw=2, label="Reines & Volonteri 2015, AGN (local)")
        ax.plot(grid, reines_volonteri15(xmsun, "ell") + np.log10(h),
                color="darkorange", lw=2, ls="--",
                label="Reines & Volonteri 2015, ellipticals (local)")
    ax.set_xlabel(r"$\log_{10}(M_\star/[M_\odot/h])$")
    ax.set_ylabel(r"$\log_{10}(M_{\rm BH}/[M_\odot/h])$")
    ax.set_title(f"M_BH-M_star, snapshot {args.snapshot:04d}, "
                 f"z={z:.3g} ({args.population})")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    out2 = args.outdir / f"mbh_mstar_compare_{args.snapshot:05d}.png"
    fig.tight_layout(); fig.savefig(out2, dpi=160); plt.close(fig)

    # ------------------- GSMF -------------------
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    hub = z = None
    xlo, xhi = [], []
    for i, (label, root) in enumerate(models):
        try:
            centers, logphi, err, header = compute_gsmf(root, args.snapshot)
        except (FileNotFoundError, ValueError) as e:
            print(f"{label}: skip gsmf ({e})"); continue
        if len(centers) == 0:
            print(f"{label}: no galaxies for gsmf"); continue
        ax.errorbar(centers, logphi, yerr=err, fmt="o", ms=4, capsize=2,
                    color=colors[i], label=label)
        hub, z = header["hubble"], header["redshift"]
        xlo.append(centers.min()); xhi.append(centers.max())
    h = _h(hub)
    if xlo:
        overlay_baldry(ax, (min(xlo), max(xhi)), h, z)
    ax.set_xlabel(r"$\log_{10}(M_\star/[M_\odot/h])$")
    ax.set_ylabel(r"$\log_{10}(\phi/[(h/{\rm Mpc})^3\,{\rm dex}^{-1}])$")
    ax.set_title(f"GSMF, snapshot {args.snapshot:04d}, z={z:.3g}")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    out3 = args.outdir / f"gsmf_compare_{args.snapshot:05d}.png"
    fig.tight_layout(); fig.savefig(out3, dpi=160); plt.close(fig)

    print(f"wrote {out1.resolve()}")
    print(f"wrote {out2.resolve()}")
    print(f"wrote {out3.resolve()}")


if __name__ == "__main__":
    main()
