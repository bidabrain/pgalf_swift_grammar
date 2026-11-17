#!/usr/bin/env python3
"""Galaxy Stellar Mass Function (GSMF) from PGalF catalogs, with calibration.

GSMF = galaxies per comoving volume per dex of stellar mass. Plotted in
sim-native units: M* in Msun/h, phi in (Mpc/h)^-3 dex^-1. M* defaults to a 3D
50 proper-kpc aperture (COLIBRE/obs convention, from GALFIND.DATA); use
--aperture-kpc to change it, or --aperture-kpc 0 for the total bound M*.
Several values (e.g. --aperture-kpc 0 30 50 100) are overplotted; GALFIND.DATA
is read only once for all of them.

Calibration (choose with --calib):
  baldry2012   : Baldry+2012 GAMA double-Schechter, z~0
  driver2022   : Driver+2022 GAMA DR4 double-Schechter, z<0.1 (total)
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

from read_pgalf_catalog import (read_fof_catalog, read_galaxy_catalog,
                                 aperture_stellar_mass, galcenter_positions)


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


def driver2022(logM_msun):
    """Driver+2022 GAMA DR4 z<0.1 total GSMF [Mpc^-3 dex^-1], Msun (h=0.7).
    Double-Schechter, Table 7 'All': logM*=10.745, logphi1*=-2.437, a1=-0.466,
    logphi2*=-3.201, a2=-1.530."""
    return _schechter(logM_msun, 10.745, 10 ** -2.437, -0.466, 10 ** -3.201, -1.530)


def davidzon2017(logM_msun, redshift):
    """Davidzon+2017 total GSMF [Mpc^-3 dex^-1] at the bin containing `redshift`.
    Returns (phi, label) or (None, None) if z is outside 0.2-5.5."""
    for zlo, zhi, lMs, a1, p1, a2, p2 in DAVIDZON2017:
        if zlo <= redshift < zhi:
            p2v = p2 * 1e-3 if p2 is not None else None
            phi = _schechter(logM_msun, lMs, p1 * 1e-3, a1, p2v, a2)
            return phi, f"Davidzon+2017 ({zlo}<z<{zhi})"
    return None, None


def aperture_centers(fof_dir, step, center, n_subs):
    """Aperture centres for `center`: None for 'com' (each subhalo's
    centre of mass from GALFIND.DATA), else the GalCenter stellar centres
    aligned to GALCATALOG order (NaN rows fall back to the COM)."""
    if center == "com":
        return None
    if center != "galcenter":
        raise ValueError(center)
    center_file = fof_dir / f"GALFIND.CENTER.{step}"
    if not center_file.exists():
        raise SystemExit(f"--center galcenter but {center_file} is missing "
                         "(run the galcenter stage)")
    cen = galcenter_positions(center_file, fof_dir / f"GALCATALOG.LIST.{step}")
    if len(cen) != n_subs:
        raise SystemExit(f"GALFIND.CENTER aligned to {len(cen)} subhalos != "
                         f"GALCATALOG ({n_subs})")
    nmiss = int((~np.isfinite(cen).all(axis=1)).sum())
    print(f"--center galcenter: {n_subs - nmiss}/{n_subs} subhalos use the GalCenter "
          f"stellar centre; {nmiss} (no entry / no stars) fall back to the COM")
    return cen


def load_stellar_masses(root, snapshot, apertures_kpc, min_npstar=0, center="com"):
    """Per-subhalo M* [Msun/h] for each entry of `apertures_kpc`. An entry >0
    is a 3D proper-kpc aperture (from GALFIND.DATA, COLIBRE/obs convention);
    <=0 or None means the catalogue's total bound M*. GALFIND.DATA is read
    once for all apertures. Subhalos with fewer than `min_npstar` bound star
    particles get M*=0 (i.e. are dropped from the GSMF) in every array.
    center='com' centres apertures on the subhalo centre of mass (all particle
    types); 'galcenter' uses the GalCenter stellar centre (GALFIND.CENTER).
    Returns (list of arrays, header)."""
    step = f"{snapshot:05d}"
    fof_dir = Path(root) / "FoF_Data" / f"FoF.{step}"
    header, _ = read_fof_catalog(fof_dir / f"FoF_halo_cat.{step}")
    _, subs, _ = read_galaxy_catalog(fof_dir / f"GALCATALOG.LIST.{step}")

    ap_kpc = [a for a in apertures_kpc if a and a > 0]
    by_ap = {}
    if ap_kpc:
        h = _h(header["hubble"])
        radii = [(a / 1000.0) / header["anow"] * h for a in ap_kpc]  # pkpc -> cMpc/h
        data_file = fof_dir / f"GALFIND.DATA.{step}"
        if not data_file.exists():
            raise SystemExit(f"--aperture-kpc set but {data_file} is missing")
        centers = aperture_centers(fof_dir, step, center, len(subs))
        mstar_ap = aperture_stellar_mass(data_file, radii, header["box_size"],
                                         centers=centers)
        if len(mstar_ap) != len(subs):
            raise SystemExit(f"GALFIND.DATA subhalos ({len(mstar_ap)}) != "
                             f"GALCATALOG ({len(subs)})")
        by_ap = {a: mstar_ap[:, j] for j, a in enumerate(ap_kpc)}
    masses = [by_ap[a] if (a and a > 0) else subs["mstar"] for a in apertures_kpc]
    if min_npstar > 0:
        keep = subs["npstar"] >= min_npstar
        masses = [np.where(keep, m, 0.0) for m in masses]
    return masses, header


def star_particle_mass(root, snapshot):
    """Median M*/npstar over galaxies with stars [Msun/h]: a typical star
    particle mass, to translate a --min-npstar cut into an M* scale."""
    step = f"{snapshot:05d}"
    _, subs, _ = read_galaxy_catalog(
        Path(root) / "FoF_Data" / f"FoF.{step}" / f"GALCATALOG.LIST.{step}")
    s = subs[subs["npstar"] > 0]
    return float(np.median(s["mstar"] / s["npstar"])) if len(s) else np.nan


def gsmf_from_masses(mstar, box_size, mass_lo=7.0, mass_hi=12.0, dlog=0.25,
                     min_count=1):
    """GSMF of galaxies with mstar>0. Returns
    (logM_centers[Msun/h], logPhi[(Mpc/h)^-3 dex^-1], err_dex)."""
    mstar = mstar[mstar > 0]                          # Msun/h
    volume = box_size ** 3                            # (Mpc/h)^3
    edges = np.arange(mass_lo, mass_hi + 1e-9, dlog)
    counts, _ = np.histogram(np.log10(mstar), bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    good = counts >= min_count
    phi = counts / (volume * dlog)                    # (Mpc/h)^-3 dex^-1
    err_dex = 0.434 / np.sqrt(np.where(counts > 0, counts, np.nan))
    return centers[good], np.log10(phi[good]), err_dex[good]


def compute_gsmf(root, snapshot, mass_lo=7.0, mass_hi=12.0, dlog=0.25, min_count=1,
                 aperture_kpc=None):
    """Sim GSMF for one aperture (see load_stellar_masses). Returns
    (logM_centers[Msun/h], logPhi[(Mpc/h)^-3 dex^-1], err_dex, header)."""
    (mstar,), header = load_stellar_masses(root, snapshot, [aperture_kpc])
    return (*gsmf_from_masses(mstar, header["box_size"], mass_lo, mass_hi, dlog,
                              min_count), header)


def _aperture_label(a):
    return f"{a:g} pkpc" if a and a > 0 else "total $M_\\star$"


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
    elif which == "driver2022":
        phi_mpc = driver2022(logM_msun)
        label = "Driver+2022 (GAMA DR4, z<0.1)"
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
    p.add_argument("--calib",
                   choices=("baldry2012", "driver2022", "davidzon2017", "none"),
                   default="davidzon2017")
    p.add_argument("--obs-file", type=Path, default=None,
                   help="logM[Msun] logPhi[Mpc^-3/dex] [err_dex]")
    p.add_argument("--aperture-kpc", type=float, nargs="+", default=[50.0],
                   help="3D stellar aperture(s) [proper kpc] from GALFIND.DATA "
                        "(COLIBRE/obs use 50); several values are overplotted, "
                        "e.g. --aperture-kpc 0 30 50 100; <=0 means total bound M*")
    p.add_argument("--min-npstar", type=int, default=0,
                   help="only count galaxies with at least this many bound star "
                        "particles (GALCATALOG npstar); applied to every aperture")
    p.add_argument("--center", choices=("com", "galcenter"), default="com",
                   help="aperture centre: 'com' = subhalo centre of mass of all "
                        "particle types; 'galcenter' = GalCenter stellar centre "
                        "(GALFIND.CENTER), falling back to the COM where missing")
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    masses, header = load_stellar_masses(args.root, args.snapshot, args.aperture_kpc,
                                         min_npstar=args.min_npstar,
                                         center=args.center)
    if args.min_npstar > 0:
        m_sp = star_particle_mass(args.root, args.snapshot)
        print(f"--min-npstar {args.min_npstar}: star particle mass ~{m_sp:.3g} Msun/h, "
              f"i.e. total M* >~ {args.min_npstar * m_sp:.3g} Msun/h "
              f"(log {np.log10(args.min_npstar * m_sp):.2f})")
    z = header["redshift"]; h = _h(header["hubble"])
    multi = len(args.aperture_kpc) > 1

    fig, ax = plt.subplots(figsize=(7, 5.5))
    xmin, xmax = np.inf, -np.inf
    for k, (a, mstar) in enumerate(zip(args.aperture_kpc, masses)):
        centers, logphi, err = gsmf_from_masses(mstar, header["box_size"],
                                                dlog=args.dlog)
        if len(centers) == 0:
            print(f"aperture {_aperture_label(a)}: no galaxies with mstar>0")
            continue
        xmin, xmax = min(xmin, centers.min()), max(xmax, centers.max())
        label = f"PGalF {Path(args.root).name}"
        if multi:
            label = f"{_aperture_label(a)}"
        ax.errorbar(centers, logphi, yerr=err, fmt="o-" if multi else "o",
                    color=f"C{k}", ms=4 if multi else 5, lw=1.2, capsize=2,
                    label=label)
    if not np.isfinite(xmin):
        raise SystemExit("No galaxies with mstar>0 found")
    overlay_calib(ax, (xmin, xmax), h, z, args.calib)
    if args.obs_file:
        overlay_obs_file(ax, args.obs_file, h)

    ax.set_xlabel(r"$\log_{10}(M_\star/[M_\odot/h])$")
    ax.set_ylabel(r"$\log_{10}(\phi/[(h/{\rm Mpc})^3\,{\rm dex}^{-1}])$")
    if multi:
        ap = f"PGalF {Path(args.root).name}"
    else:
        a = args.aperture_kpc[0]
        ap = f"{a:g} pkpc aperture" if a and a > 0 else "total $M_\\star$"
    if args.min_npstar > 0:
        ap += f", $N_\\star\\geq{args.min_npstar}$"
    if args.center != "com" and any(a and a > 0 for a in args.aperture_kpc):
        ap += f", centre={args.center}"
    ax.set_title(f"GSMF  snapshot {args.snapshot:04d}, z={z:.3g}  ({ap})")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    fig.tight_layout()
    out = args.output or Path(f"gsmf_{args.snapshot:05d}.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160); plt.close(fig)
    print(f"z={z:.4g}; calib={args.calib}; wrote {out.resolve()}")


if __name__ == "__main__":
    main()
