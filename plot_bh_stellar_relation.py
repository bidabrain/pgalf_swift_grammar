#!/usr/bin/env python3
"""Plot the black-hole mass -- stellar mass relation from PGALF catalogs.

M_BH is taken from the subhalo 'msink' field (BHs are sink particles in
SWIFT/EAGLE); M_star from 'mstar'. One point per galaxy.

Calibration overlay: Reines & Volonteri (2015), ApJ 813, 82.
  AGN (total stellar mass):   log(M_BH) = 7.45 + 1.05 log(M*/1e11)
  Ellipticals/classical bulges: log(M_BH) = 8.95 + 1.40 log(M*/1e11)
  (masses in Msun; converted to Msun/h here to match the catalog units.)
  NOTE: coefficients are the canonical published values -- double-check
        against the paper before using in a publication.

Single snapshot:
    python plot_bh_stellar_relation.py --root <OUT> --snapshot 35
Several snapshots overlaid, coloured by redshift:
    python plot_bh_stellar_relation.py --root <OUT> --snapshots 20-35
"""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_pgalf_catalog import read_fof_catalog, read_galaxy_catalog


def reines_volonteri15(log_mstar_msun, sample="agn"):
    """R&V15 mean relation. Input/output in log10 of Msun (no little-h)."""
    if sample == "agn":
        alpha, beta = 7.45, 1.05
    else:  # ellipticals + classical bulges
        alpha, beta = 8.95, 1.40
    return alpha + beta * (np.asarray(log_mstar_msun, float) - 11.0)


def load_bh_stellar(root, snapshot, population):
    """Return (header, M_star[], M_BH[]) in Msun/h for one snapshot."""
    step = f"{snapshot:05d}"
    fof_dir = Path(root) / "FoF_Data" / f"FoF.{step}"
    header, _ = read_fof_catalog(fof_dir / f"FoF_halo_cat.{step}")
    hosts, subs, host_ids = read_galaxy_catalog(fof_dir / f"GALCATALOG.LIST.{step}")

    mstar, mbh = [], []
    for host_id in range(len(hosts)):
        members = subs[host_ids == host_id]
        members = members[(members["mstar"] > 0) & (members["msink"] > 0)]
        if len(members) == 0:
            continue
        if population == "central":          # most massive stellar subhalo = central
            members = members[[int(np.argmax(members["mstar"]))]]
        mstar.extend(members["mstar"])
        mbh.extend(members["msink"])
    return header, np.asarray(mstar), np.asarray(mbh)


def expand_snapshots(spec):
    """'20-35' -> [20..35]; '20,25,30' -> [20,25,30]; '35' -> [35]."""
    out = []
    for part in str(spec).split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def add_rv15_lines(ax, hubble):
    """Overplot R&V15 AGN + elliptical relations, converted to Msun/h axes."""
    h = hubble / 100.0 if hubble > 1.0 else hubble      # H0 (km/s/Mpc) -> h
    lo, hi = ax.get_xlim()
    xh = np.linspace(lo, hi, 100)                        # log10(M*/[Msun/h])
    xs = xh - np.log10(h)                                # -> log10(M*/[Msun])
    for sample, style, lab in [
        ("agn", dict(color="crimson", lw=2),
         "Reines & Volonteri (2015), AGN"),
        ("ell", dict(color="darkorange", lw=2, ls="--"),
         "Reines & Volonteri (2015), ellipticals"),
    ]:
        yh = reines_volonteri15(xs, sample) + np.log10(h)   # -> Msun/h
        ax.plot(xh, yh, zorder=3, label=lab, **style)
    ax.set_xlim(lo, hi)
    ax.legend(frameon=False, loc="upper left", fontsize=9)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True,
                   help="PGALF run dir containing FoF_Data")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--snapshot", type=int, help="single snapshot number")
    g.add_argument("--snapshots", type=str,
                   help="several snapshots: 'a-b' or 'a,b,c' (overlaid, coloured by z)")
    p.add_argument("--population", choices=("central", "subhalo"), default="central",
                   help="most massive stellar subhalo per host (central) or all subhalos")
    p.add_argument("--no-calib", action="store_true", help="do not draw R&V15 lines")
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    snaps = [args.snapshot] if args.snapshot is not None else expand_snapshots(args.snapshots)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    npts_total = 0
    hubble = None

    if len(snaps) == 1:
        header, mstar, mbh = load_bh_stellar(args.root, snaps[0], args.population)
        hubble = header["hubble"]
        if len(mstar):
            ax.scatter(np.log10(mstar), np.log10(mbh), s=10, alpha=0.55, linewidths=0)
            npts_total = len(mstar)
        title = (f"PGALF M_BH-M_star  snapshot {snaps[0]:04d}, "
                 f"z={header['redshift']:.3g} ({args.population})")
    else:
        data = []
        for s in snaps:
            try:
                header, mstar, mbh = load_bh_stellar(args.root, s, args.population)
            except (FileNotFoundError, ValueError):
                print(f"skip snapshot {s}: catalog missing/unreadable")
                continue
            if len(mstar):
                data.append((header["redshift"], mstar, mbh))
                npts_total += len(mstar)
                hubble = header["hubble"]
        if not data:
            raise SystemExit("No usable snapshots with BH+stellar galaxies")
        zs = [d[0] for d in data]
        norm = plt.Normalize(min(zs), max(zs))
        cmap = plt.cm.viridis
        for z, mstar, mbh in data:
            ax.scatter(np.log10(mstar), np.log10(mbh), s=10, alpha=0.55,
                       linewidths=0, color=cmap(norm(z)))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
        fig.colorbar(sm, ax=ax, label="redshift z")
        title = (f"PGALF M_BH-M_star  snapshots {snaps[0]:04d}-{snaps[-1]:04d} "
                 f"({args.population})")

    if npts_total == 0:
        raise SystemExit("No galaxies with both M_BH>0 and M_star>0 found")

    ax.set_xlabel(r"$\log_{10}(M_\star/[M_\odot/h])$")
    ax.set_ylabel(r"$\log_{10}(M_{\rm BH}/[M_\odot/h])$")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    if not args.no_calib and hubble is not None:
        add_rv15_lines(ax, hubble)
    fig.tight_layout()

    if args.output is None:
        tag = f"{snaps[0]:05d}" if len(snaps) == 1 else f"{snaps[0]:05d}-{snaps[-1]:05d}"
        args.output = Path(f"mbh_mstar_{args.population}_{tag}.png")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    print(f"points={npts_total}  output={args.output.resolve()}")


if __name__ == "__main__":
    main()
