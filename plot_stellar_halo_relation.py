#!/usr/bin/env python3
"""Plot the z=0 stellar-to-halo mass relation from PGALF catalogs.

By default, one point is plotted per PGALF host: the most massive stellar
subhalo is used as the central galaxy.  Use --population subhalo to plot all
stellar-mass subhalos instead.
"""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from read_pgalf_catalog import read_fof_catalog, read_galaxy_catalog


def load_relation(root, snapshot, population):
    """Return host-halo and stellar masses in Msun/h."""
    step = f"{snapshot:05d}"
    fof_dir = Path(root) / "FoF_Data" / f"FoF.{step}"
    fof_file = fof_dir / f"FoF_halo_cat.{step}"
    galaxy_file = fof_dir / f"GALCATALOG.LIST.{step}"

    header, _ = read_fof_catalog(fof_file)
    hosts, subhalos, host_ids = read_galaxy_catalog(galaxy_file)

    halo_masses = []
    stellar_masses = []
    for host_id, host in enumerate(hosts):
        members = subhalos[host_ids == host_id]
        members = members[members["mstar"] > 0]
        if len(members) == 0 or host["totm"] <= 0:
            continue

        if population == "central":
            members = members[[int(np.argmax(members["mstar"]))]]

        halo_masses.extend([host["totm"]] * len(members))
        stellar_masses.extend(members["mstar"])

    return header, np.asarray(halo_masses), np.asarray(stellar_masses)


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
    args = parser.parse_args()

    header, halo_mass, stellar_mass = load_relation(
        args.root, args.snapshot, args.population
    )
    if len(halo_mass) == 0:
        raise SystemExit("No positive-stellar-mass PGALF galaxies found")

    output = args.output or Path(
        f"mstar_mhalo_{args.population}_{args.snapshot:05d}.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    redshift = header["redshift"]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(
        np.log10(halo_mass),
        np.log10(stellar_mass),
        s=10,
        alpha=0.55,
        linewidths=0,
    )
    # Some PGALF headers store H0 in km/s/Mpc rather than the dimensionless h.
    hubble = header["hubble"]
    if hubble > 1.0:
        hubble /= 100.0

    log_mhalo_grid = np.linspace(
        np.log10(halo_mass).min(), np.log10(halo_mass).max(), 200
    )
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
    ax.set_title(
        f"PGALF stellar-to-halo mass relation "
        f"snapshot {args.snapshot:04d}, z={redshift:.3g} ({args.population})"
    )
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)

    print(f"snapshot={args.snapshot:04d} z={redshift:.6g}")
    print(f"population={args.population} points={len(halo_mass)}")
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()
