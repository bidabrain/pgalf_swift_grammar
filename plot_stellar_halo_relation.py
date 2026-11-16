#!/usr/bin/env python3
"""Plot the z=0 stellar-to-halo mass relation from PGALF catalogs.

By default, one point is plotted per PGALF host: the most massive stellar
subhalo is used as the central galaxy.  Use --population subhalo to plot all
stellar-mass subhalos instead.

Several --aperture-kpc values (e.g. 0 30 50 100) are overplotted, each as a
light scatter plus its binned median; GALFIND.DATA is read only once.
"""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_pgalf_catalog import read_galaxy_catalog
from plot_gsmf import load_stellar_masses


def _select(hosts, host_ids, mstar_all, population):
    """Pair subhalo M* with its host's total mass. Keeps mstar>0 subhalos of
    hosts with totm>0; 'central' keeps only the most massive one per host."""
    host_m = hosts["totm"][host_ids]
    keep = (mstar_all > 0) & (host_m > 0)
    if population == "central":
        idx = np.flatnonzero(keep)
        # host_ids is non-decreasing (GALCATALOG order): sort by (host, -M*)
        order = idx[np.lexsort((-mstar_all[idx], host_ids[idx]))]
        first = np.ones(len(order), dtype=bool)
        first[1:] = host_ids[order][1:] != host_ids[order][:-1]
        sel = order[first]
    else:
        sel = np.flatnonzero(keep)
    return host_m[sel], mstar_all[sel]


def load_relations(root, snapshot, population, apertures_kpc):
    """Return (header, [(halo_mass, stellar_mass), ...]) in Msun/h, one pair
    per entry of apertures_kpc.

    Halo mass is the host total mass (no aperture). Stellar mass uses a 3D
    proper-kpc aperture (from GALFIND.DATA) when the entry is >0, else the
    catalogue's total bound M*."""
    step = f"{snapshot:05d}"
    galaxy_file = Path(root) / "FoF_Data" / f"FoF.{step}" / f"GALCATALOG.LIST.{step}"
    hosts, _, host_ids = read_galaxy_catalog(galaxy_file)
    masses, header = load_stellar_masses(root, snapshot, apertures_kpc)
    return header, [_select(hosts, host_ids, m, population) for m in masses]


def load_relation(root, snapshot, population, aperture_kpc=None):
    """Single-aperture load_relations: (header, halo_mass, stellar_mass)."""
    header, [(halo_mass, stellar_mass)] = load_relations(
        root, snapshot, population, [aperture_kpc])
    return header, halo_mass, stellar_mass


def binned_median(log_x, log_y, dlog=0.2, min_count=5):
    """Median of log_y in log_x bins with >= min_count objects."""
    edges = np.arange(np.floor(log_x.min() / dlog) * dlog, log_x.max() + dlog, dlog)
    idx = np.digitize(log_x, edges) - 1
    xc, ym = [], []
    for i in range(len(edges) - 1):
        sel = idx == i
        if sel.sum() >= min_count:
            xc.append(0.5 * (edges[i] + edges[i + 1]))
            ym.append(np.median(log_y[sel]))
    return np.asarray(xc), np.asarray(ym)


def _aperture_label(a):
    return f"{a:g} pkpc" if a and a > 0 else "total $M_\\star$"


def moster13(log_mhalo, redshift=0.0):
    """Moster, Naab & White (2013) mean stellar-to-halo mass relation.

    Parameters and units follow the paper: log_mhalo is log10(M_200c/Msun)
    and the return value is log10(M_star/Msun), i.e. no little-h factors.
    """
    zz = redshift / (1.0 + redshift)
    log_m1 = 11.590 + 1.195 * zz
    norm = 0.0351 - 0.0247 * zz
    beta = 1.376 - 0.826 * zz
    gamma = 0.608 + 0.329 * zz
    x = 10.0 ** (np.asarray(log_mhalo, dtype=float) - log_m1)
    return log_mhalo + np.log10(2.0 * norm / (x ** (-beta) + x ** gamma))


# ---- Behroozi+2019 (UniverseMachine DR1) central SHMR --------------------
_UM_NAMES = ("EFF_0 EFF_0_A EFF_0_A2 EFF_0_Z M_1 M_1_A M_1_A2 M_1_Z "
             "ALPHA ALPHA_A ALPHA_A2 ALPHA_Z BETA BETA_A BETA_Z "
             "DELTA GAMMA GAMMA_A GAMMA_Z").split()


def _um_param_file(um_dir, population):
    """Locate the UniverseMachine SMHM parameter file under um_dir."""
    name = "smhm_med_cen_params.txt" if population == "central" else "smhm_med_params.txt"
    p = Path(um_dir)
    if p.is_file():
        return p
    for cand in (p / name, p / "params" / name,
                 p / "data" / "smhm" / "params" / name, p / "smhm" / "params" / name):
        if cand.is_file():
            return cand
    return None


def behroozi19_um(param_file, redshift):
    """UniverseMachine DR1 median SMHM at a given z (masses in Msun, no h).

    Returns (log_mpeak, log_mstar) arrays, or (None, None) if z is beyond the
    tabulated range (the fit upper mass limit 14.5-0.35z drops below 10.5).
    Formula reproduced from data/smhm/params/gen_smhm.py.
    """
    vals = []
    with open(param_file) as fh:
        for line in fh:
            tok = line.split()
            if not tok or tok[0].startswith("#"):
                continue
            vals.append(float(tok[1]))
            if len(vals) == len(_UM_NAMES):
                break
    p = dict(zip(_UM_NAMES, vals))
    z = float(redshift)
    a1 = 1.0 / (1.0 + z) - 1.0
    lna = np.log(1.0 / (1.0 + z))
    m_1 = p["M_1"] + a1 * p["M_1_A"] - lna * p["M_1_A2"] + z * p["M_1_Z"]
    sm_0 = m_1 + p["EFF_0"] + a1 * p["EFF_0_A"] - lna * p["EFF_0_A2"] + z * p["EFF_0_Z"]
    alpha = p["ALPHA"] + a1 * p["ALPHA_A"] - lna * p["ALPHA_A2"] + z * p["ALPHA_Z"]
    beta = p["BETA"] + a1 * p["BETA_A"] + z * p["BETA_Z"]
    delta = p["DELTA"]
    gamma = 10.0 ** (p["GAMMA"] + a1 * p["GAMMA_A"] + z * p["GAMMA_Z"])
    m_hi = 14.5 - 0.35 * z
    if m_hi <= 10.5:
        return None, None
    m = np.arange(10.5, m_hi + 1e-9, 0.05)          # log10(Mpeak/Msun)
    dm = m - m_1
    sm = sm_0 - np.log10(10 ** (-alpha * dm) + 10 ** (-beta * dm)) \
        + gamma * np.exp(-0.5 * (dm / delta) ** 2)  # log10(M*/Msun)
    return m, sm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/gpfs/bkoh/eagle_pgalf/eagle"),
        help="PGALF run directory containing FoF_Data",
    )
    parser.add_argument("--snapshot", type=int, default=274)
    parser.add_argument(
        "--population",
        choices=("central", "subhalo"),
        default="central",
        help="Plot the most massive stellar subhalo per host or all subhalos",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="PNG output path (default: mstar_mhalo_<population>_<snapshot>.png)",
    )
    parser.add_argument(
        "--um-dir",
        type=Path,
        default=None,
        help="UniverseMachine DR1 dir (or params dir/file) to overlay "
             "Behroozi+2019 central SHMR",
    )
    parser.add_argument(
        "--aperture-kpc",
        type=float,
        nargs="+",
        default=[50.0],
        help="3D stellar aperture(s) [proper kpc] from GALFIND.DATA (COLIBRE/obs "
             "use 50); several values are overplotted with binned medians, "
             "e.g. --aperture-kpc 0 30 50 100; <=0 means total bound M*. "
             "Halo mass is unaffected.",
    )
    parser.add_argument("--dlog", type=float, default=0.2,
                        help="halo-mass bin width [dex] for the medians "
                             "(multi-aperture only)")
    args = parser.parse_args()

    header, relations = load_relations(
        args.root, args.snapshot, args.population, args.aperture_kpc
    )
    if all(len(hm) == 0 for hm, _ in relations):
        raise SystemExit("No positive-stellar-mass PGALF galaxies found")
    multi = len(args.aperture_kpc) > 1

    output = args.output or Path(
        f"mstar_mhalo_{args.population}_{args.snapshot:05d}.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    redshift = header["redshift"]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for k, (a, (halo_mass, stellar_mass)) in enumerate(
            zip(args.aperture_kpc, relations)):
        if len(halo_mass) == 0:
            print(f"aperture {_aperture_label(a)}: no positive-stellar-mass galaxies")
            continue
        lx, ly = np.log10(halo_mass), np.log10(stellar_mass)
        if not multi:
            ax.scatter(lx, ly, s=10, alpha=0.55, linewidths=0)
            continue
        ax.scatter(lx, ly, s=4, alpha=0.2, linewidths=0, color=f"C{k}")
        xc, ym = binned_median(lx, ly, args.dlog)
        ax.plot(xc, ym, "o-", color=f"C{k}", ms=4, lw=2, zorder=4,
                label=f"{_aperture_label(a)} (median)")
    # Some PGALF headers store H0 in km/s/Mpc rather than the dimensionless h.
    hubble = header["hubble"]
    if hubble > 1.0:
        hubble /= 100.0

    all_lx = np.log10(np.concatenate([hm for hm, _ in relations]))
    log_mhalo_grid = np.linspace(all_lx.min(), all_lx.max(), 200)
    ax.plot(
        log_mhalo_grid,
        moster13(log_mhalo_grid - np.log10(hubble), redshift) + np.log10(hubble),
        color="crimson",
        lw=2,
        zorder=3,
        label="Moster, Naab & White (2013), MNRAS 428, 3121",
    )

    # Behroozi+2019 (UniverseMachine) central SHMR, converted Msun -> Msun/h
    if args.um_dir is not None:
        pf = _um_param_file(args.um_dir, args.population)
        if pf is None:
            print(f"warning: UniverseMachine params not found under {args.um_dir}")
        else:
            m_um, sm_um = behroozi19_um(pf, redshift)
            if m_um is None:
                print(f"UniverseMachine: z={redshift:.3g} beyond fit range (z>~11), "
                      "skipping overlay")
            else:
                ax.plot(
                    m_um + np.log10(hubble),
                    sm_um + np.log10(hubble),
                    color="navy", lw=2, ls="--", zorder=3,
                    label="Behroozi+2019 (UniverseMachine DR1)",
                )

    ax.legend(frameon=False, loc="upper left", fontsize=9)

    ax.set_xlabel(r"$\log_{10}(M_{\rm halo}/[M_\odot/h])$")
    ax.set_ylabel(r"$\log_{10}(M_\star/[M_\odot/h])$")
    if multi:
        ap = "multi-aperture"
    else:
        a = args.aperture_kpc[0]
        ap = f"{a:g} pkpc aperture" if a and a > 0 else "total $M_\\star$"
    ax.set_title(
        f"PGALF stellar-to-halo mass relation "
        f"snapshot {args.snapshot:04d}, z={redshift:.3g}\n({args.population}, {ap})"
    )
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)

    print(f"snapshot={args.snapshot:04d} z={redshift:.6g}")
    for a, (halo_mass, _) in zip(args.aperture_kpc, relations):
        print(f"population={args.population} aperture={_aperture_label(a)} "
              f"points={len(halo_mass)}")
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()
