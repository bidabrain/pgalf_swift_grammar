#!/usr/bin/env python3
"""Read PGALF FoF halo and galaxy/subhalo catalogs."""

from pathlib import Path
import argparse

import numpy as np


# These dtypes match the little-endian Linux binaries used for this run.
# The explicit padding fields are required to match the C struct sizes.
FOF_HEADER_DTYPE = np.dtype(
    [
        ("box_size", "<f4"),
        ("hubble", "<f4"),
        ("omega_m", "<f4"),
        ("omega_b", "<f4"),
        ("omega_l", "<f4"),
        ("amax", "<f4"),
        ("anow", "<f4"),
    ]
)

FOF_HALO_DTYPE = np.dtype(
    [
        ("np", "<u8"),
        ("npstar", "<u8"),
        ("npgas", "<u8"),
        ("npdm", "<u8"),
        ("npsink", "<u8"),
        ("x", "<f8"),
        ("y", "<f8"),
        ("z", "<f8"),
        ("mass", "<f8"),
        ("mstar", "<f8"),
        ("mgas", "<f8"),
        ("mdm", "<f8"),
        ("msink", "<f8"),
        ("vx", "<f4"),
        ("vy", "<f4"),
        ("vz", "<f4"),
        ("_pad", "V4"),
    ]
)

PGALF_HALO_DTYPE = np.dtype(
    [
        ("nsub", "<i4"),
        ("ndm", "<i4"),
        ("nstar", "<i4"),
        ("nsink", "<i4"),
        ("ngas", "<i4"),
        ("npall", "<i4"),
        ("totm", "<f8"),
        ("mdm", "<f8"),
        ("mgas", "<f8"),
        ("msink", "<f8"),
        ("mstar", "<f8"),
        ("x", "<f8"),
        ("y", "<f8"),
        ("z", "<f8"),
        ("vx", "<f8"),
        ("vy", "<f8"),
        ("vz", "<f8"),
    ]
)

PGALF_SUBHALO_DTYPE = np.dtype(
    [
        ("npdm", "<i4"),
        ("npgas", "<i4"),
        ("npsink", "<i4"),
        ("npstar", "<i4"),
        ("npall", "<i4"),
        ("_pad", "V4"),
        ("totm", "<f8"),
        ("mdm", "<f8"),
        ("mgas", "<f8"),
        ("msink", "<f8"),
        ("mstar", "<f8"),
        ("x", "<f8"),
        ("y", "<f8"),
        ("z", "<f8"),
        ("vx", "<f8"),
        ("vy", "<f8"),
        ("vz", "<f8"),
    ]
)

GALCENTER_DTYPE = np.dtype(
    [
        ("halo_id", "<i4"),
        ("subhalo_id", "<i4"),
        ("gal", "<f8", (4,)),
        ("dmhalo", "<f8", (4,)),
        ("gas", "<f8", (4,)),
    ]
)


# GALFIND.DATA layout. Halo and subhalo headers are both 112 bytes
# (6*int32 + 11*float64). Each subhalo is followed by its DM, gas, sink then
# star particle records (gfind.c writes sizeof(DmType)/GasType/SinkType/StarType
# blocks in that order). Those struct sizes depend on the gfind build flags
# (NCHEM, NDUST, ...), so the DM/star stride (DmType==StarType==PmType) and the
# gas stride are AUTO-DETECTED from the file size + GALCATALOG.LIST counts and
# validated against the catalogue's mstar -- never hard-coded. SinkType has no
# conditional fields (fixed 168 bytes). In every build the first PmType fields
# are x,y,z,vx,vy,vz,mass, so a star's position is at byte 0 and its mass at
# byte 48 of its record regardless of flags.
_GALFIND_HDR_BYTES = 112
_SINK_BYTES = 168
_PM_POS_OFFSET = 0
_PM_MASS_OFFSET = 48
_COUNT_CAP = 200_000_000        # reject absurd counts from a misaligned read


def _galcatalog_for(data_file):
    p = Path(data_file)
    return p.with_name(p.name.replace("GALFIND.DATA", "GALCATALOG.LIST"))


def _iter_galfind(filename, s_pm, s_gas, limit=None):
    """Yield (com, star_pos[N,3], star_mass[N]) per subhalo given the DM/star
    stride `s_pm` and gas stride `s_gas`. Stops after `limit` subhalos if set.
    Raises ValueError on misalignment (absurd counts / truncation)."""
    size = Path(filename).stat().st_size
    seen = 0
    with open(filename, "rb") as f:
        while f.tell() < size:
            halo_hdr = f.read(_GALFIND_HDR_BYTES)
            if len(halo_hdr) < _GALFIND_HDR_BYTES:
                break
            nsub = int(np.frombuffer(halo_hdr, dtype="<i4", count=1)[0])
            if not (0 <= nsub <= _COUNT_CAP):
                raise ValueError("misaligned halo header")
            for _ in range(nsub):
                hdr = f.read(_GALFIND_HDR_BYTES)
                if len(hdr) < _GALFIND_HDR_BYTES:
                    raise ValueError(f"Truncated subhalo header: {filename}")
                c = np.frombuffer(hdr, dtype="<i4", count=6)
                ndm, ngas, nsink, nst = int(c[0]), int(c[1]), int(c[2]), int(c[3])
                if min(ndm, ngas, nsink, nst) < 0 or max(ndm, ngas, nsink, nst) > _COUNT_CAP:
                    raise ValueError("misaligned subhalo header")
                com = np.frombuffer(hdr[24:112], dtype="<f8", count=11)[5:8].copy()
                f.seek(ndm * s_pm + ngas * s_gas + nsink * _SINK_BYTES, 1)
                block = f.read(nst * s_pm)
                if len(block) < nst * s_pm:
                    raise ValueError(f"Truncated star block: {filename}")
                if nst:
                    rec = np.frombuffer(block, dtype=np.uint8).reshape(nst, s_pm)
                    pos = rec[:, _PM_POS_OFFSET:_PM_POS_OFFSET + 24].copy().view("<f8").reshape(nst, 3)
                    mass = rec[:, _PM_MASS_OFFSET:_PM_MASS_OFFSET + 8].copy().view("<f8").reshape(nst)
                else:
                    pos = np.empty((0, 3), dtype="<f8")
                    mass = np.empty(0, dtype="<f8")
                yield com, pos, mass
                seen += 1
                if limit is not None and seen >= limit:
                    return


def _detect_galfind_sizes(filename):
    """Return (s_pm, s_gas) for GALFIND.DATA by matching the file size to the
    GALCATALOG.LIST particle counts and validating star-mass sums vs mstar."""
    cat = _galcatalog_for(filename)
    if not cat.exists():
        raise ValueError(f"need {cat} alongside {filename} to size its records")
    hosts, galaxies, _ = read_galaxy_catalog(cat)
    n_halos, n_subs = len(hosts), len(galaxies)
    sum_dm = int(galaxies["npdm"].sum())
    sum_gas = int(galaxies["npgas"].sum())
    sum_sink = int(galaxies["npsink"].sum())
    sum_star = int(galaxies["npstar"].sum())

    size = Path(filename).stat().st_size
    base = (n_halos + n_subs) * _GALFIND_HDR_BYTES + sum_sink * _SINK_BYTES

    star_idx = [i for i in range(n_subs) if galaxies["npstar"][i] > 0]
    target = {i: float(galaxies["mstar"][i]) for i in star_idx[:30]}
    n_val = (max(target) + 1) if target else 0

    for s_pm in range(96, 520, 8):
        rem = size - base - s_pm * (sum_dm + sum_star)
        if sum_gas > 0:
            if rem <= 0 or rem % sum_gas != 0:
                continue
            s_gas = rem // sum_gas
            if not (96 <= s_gas <= 520):
                continue
        else:
            if rem != 0:
                continue
            s_gas = 104
        if _validate_sizes(filename, s_pm, s_gas, target, n_val):
            return int(s_pm), int(s_gas)
    raise ValueError(f"could not determine GALFIND.DATA record sizes for {filename}")


def _validate_sizes(filename, s_pm, s_gas, target, n_val, rtol=1e-3):
    if not target:
        return True
    try:
        for i, (_com, _pos, mass) in enumerate(
                _iter_galfind(filename, s_pm, s_gas, limit=n_val)):
            want = target.get(i)
            if want and want > 0 and abs(float(mass.sum()) - want) > rtol * want:
                return False
    except (ValueError, MemoryError, OverflowError):
        return False
    return True


def read_galfind_star_particles(filename, sizes=None):
    """Iterate (com_xyz, star_pos, star_mass) per subhalo, in GALCATALOG.LIST
    order. Record strides are auto-detected (sizes=(s_pm, s_gas) overrides).
    Positions and masses are in the catalogue's little-h code units."""
    s_pm, s_gas = sizes if sizes else _detect_galfind_sizes(filename)
    yield from _iter_galfind(filename, s_pm, s_gas)


def aperture_stellar_mass(filename, radius, box_size, centers=None, sizes=None):
    """Stellar mass within `radius` (code units) of each subhalo centre,
    periodic in `box_size`. Returns an array aligned with GALCATALOG order.

    radius, box_size : same (little-h, comoving) length units as the positions.
                       `radius` may be a sequence: the file is read once and an
                       (N, len(radius)) array is returned instead of (N,).
    centers          : optional (N, 3) centres; default = each subhalo's COM.
    sizes            : optional (s_pm, s_gas) record strides; auto-detected if None.
    """
    if sizes is None:
        sizes = _detect_galfind_sizes(filename)
    scalar = np.ndim(radius) == 0
    r2 = np.atleast_1d(np.asarray(radius, dtype=np.float64)) ** 2
    out = []
    for i, (com, pos, mass) in enumerate(read_galfind_star_particles(filename, sizes)):
        if len(mass) == 0:
            out.append(np.zeros(len(r2)))
            continue
        c = com if centers is None else centers[i]
        dx = pos - c
        dx -= box_size * np.round(dx / box_size)        # periodic minimum image
        d2 = (dx * dx).sum(axis=1)
        out.append([float(mass[d2 < r].sum()) for r in r2])
    out = np.asarray(out, dtype=np.float64).reshape(-1, len(r2))
    return out[:, 0] if scalar else out


def read_fof_catalog(filename):
    """Return the FoF header and all FoF halo records."""
    with open(filename, "rb") as stream:
        header_raw = np.fromfile(stream, dtype=FOF_HEADER_DTYPE, count=1)
        if len(header_raw) != 1:
            raise ValueError(f"Incomplete FoF header: {filename}")
        halos = np.fromfile(stream, dtype=FOF_HALO_DTYPE)

    header = {name: float(header_raw[0][name]) for name in FOF_HEADER_DTYPE.names}
    header["redshift"] = 1.0 / header["anow"] - 1.0
    return header, halos


def read_galaxy_catalog(filename):
    """Read GALCATALOG.LIST and return host summaries, galaxies, and host IDs.

    The returned host ID is the zero-based record number within GALCATALOG.LIST.
    It is not necessarily the row number in FoF_halo_cat, because PGALF can omit
    FoF halos that do not produce an accepted galaxy/subhalo catalog entry.
    """
    host_records = []
    galaxy_records = []
    host_ids = []

    with open(filename, "rb") as stream:
        host_id = 0
        while True:
            host = np.fromfile(stream, dtype=PGALF_HALO_DTYPE, count=1)
            if len(host) == 0:
                break
            if len(host) != 1:
                raise ValueError(f"Truncated host-halo record: {filename}")

            nsub = int(host[0]["nsub"])
            if nsub < 0:
                raise ValueError(f"Invalid nsub={nsub} in {filename}")

            galaxies = np.fromfile(
                stream, dtype=PGALF_SUBHALO_DTYPE, count=nsub
            )
            if len(galaxies) != nsub:
                raise ValueError(f"Truncated galaxy records: {filename}")

            host_records.append(host[0])
            galaxy_records.extend(galaxies)
            host_ids.extend([host_id] * nsub)
            host_id += 1

    hosts = np.asarray(host_records, dtype=PGALF_HALO_DTYPE)
    if galaxy_records:
        galaxies = np.asarray(galaxy_records, dtype=PGALF_SUBHALO_DTYPE)
    else:
        galaxies = np.empty(0, dtype=PGALF_SUBHALO_DTYPE)
    return hosts, galaxies, np.asarray(host_ids, dtype=np.int64)


def read_galaxy_centers(filename):
    """Read GALFIND.CENTER; each center is x, y, z, R."""
    centers = np.fromfile(filename, dtype=GALCENTER_DTYPE)
    expected_bytes = len(centers) * GALCENTER_DTYPE.itemsize
    actual_bytes = Path(filename).stat().st_size
    if expected_bytes != actual_bytes:
        raise ValueError(f"Truncated or incompatible center file: {filename}")
    return centers


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=int, help="Snapshot number, e.g. 274")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/gpfs/bkoh/eagle_pgalf/eagle"),
        help="PGALF run directory",
    )
    args = parser.parse_args()

    step = f"{args.snapshot:05d}"
    fof_dir = args.root / "FoF_Data" / f"FoF.{step}"
    fof_file = fof_dir / f"FoF_halo_cat.{step}"
    galaxy_file = fof_dir / f"GALCATALOG.LIST.{step}"
    center_file = fof_dir / f"GALFIND.CENTER.{step}"

    header, fof_halos = read_fof_catalog(fof_file)
    pgalf_hosts, galaxies, host_ids = read_galaxy_catalog(galaxy_file)

    print(f"snapshot={step}  z={header['redshift']:.6g}")
    print(f"FoF halos={len(fof_halos)}")
    if len(fof_halos):
        i = int(np.argmax(fof_halos["mass"]))
        print(
            "most massive FoF halo: "
            f"M={fof_halos[i]['mass']:.6e} Msun/h  "
            f"N={fof_halos[i]['np']}  "
            f"pos=({fof_halos[i]['x']:.5f}, "
            f"{fof_halos[i]['y']:.5f}, {fof_halos[i]['z']:.5f})"
        )

    print(f"PGALF host halos={len(pgalf_hosts)}  galaxies/subhalos={len(galaxies)}")
    if len(galaxies):
        i = int(np.argmax(galaxies["mstar"]))
        print(
            "most massive stellar galaxy: "
            f"Mstar={galaxies[i]['mstar']:.6e} Msun/h  "
            f"Mtot={galaxies[i]['totm']:.6e} Msun/h  "
            f"host={host_ids[i]}  "
            f"pos=({galaxies[i]['x']:.5f}, "
            f"{galaxies[i]['y']:.5f}, {galaxies[i]['z']:.5f})"
        )

    if center_file.exists():
        centers = read_galaxy_centers(center_file)
        print(f"refined galaxy centers={len(centers)}")


if __name__ == "__main__":
    main()