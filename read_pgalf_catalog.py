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


# GALFIND.DATA per-particle layout. Halo and subhalo headers are both
# 6*int32 + 11*float64 = 112 bytes. Particle record sizes (bytes):
_GALFIND_HDR_BYTES = 6 * 4 + 11 * 8        # 112
_GALFIND_DM_BYTES = 96                      # == star record size
_GALFIND_GAS_BYTES = 104
_GALFIND_SINK_BYTES = 168

# One star particle, matching read_PGalF_v1.py's STAR record.
STAR_RECORD_DTYPE = np.dtype(
    [
        ("pos", "<f8", (3,)),
        ("vel", "<f8", (3,)),
        ("mass", "<f8"),
        ("mass0", "<f8"),
        ("tp", "<f8"),
        ("zp", "<f8"),
        ("id", "<i4"),
        ("family", "i1"),
        ("tag", "i1"),
        ("_pad", "<i2"),
        ("level", "<i4"),
        ("_dum", "<i4"),
    ]
)
assert STAR_RECORD_DTYPE.itemsize == _GALFIND_DM_BYTES


def read_galfind_star_particles(filename):
    """Iterate GALFIND.DATA, yielding (com_xyz, star_pos, star_mass) per subhalo
    in the SAME order as GALCATALOG.LIST. DM/gas/sink blocks are skipped.

    com_xyz : the subhalo centre-of-mass stored in its header (code units)
    star_pos: (N, 3) star positions, star_mass: (N,) star masses
    Positions and masses are in the catalogue's little-h code units.
    """
    path = Path(filename)
    size = path.stat().st_size
    with open(path, "rb") as f:
        while f.tell() < size:
            halo_hdr = f.read(_GALFIND_HDR_BYTES)
            if len(halo_hdr) < _GALFIND_HDR_BYTES:
                break
            nsub = int(np.frombuffer(halo_hdr, dtype="<i4", count=1)[0])
            for _ in range(nsub):
                hdr = f.read(_GALFIND_HDR_BYTES)
                if len(hdr) < _GALFIND_HDR_BYTES:
                    raise ValueError(f"Truncated subhalo header: {filename}")
                counts = np.frombuffer(hdr, dtype="<i4", count=6)
                ndm, ngas, nsink, nst = (int(counts[0]), int(counts[1]),
                                         int(counts[2]), int(counts[3]))
                com = np.frombuffer(hdr[24:24 + 88], dtype="<f8", count=11)[5:8].copy()
                f.seek(ndm * _GALFIND_DM_BYTES
                       + ngas * _GALFIND_GAS_BYTES
                       + nsink * _GALFIND_SINK_BYTES, 1)
                stars = np.fromfile(f, dtype=STAR_RECORD_DTYPE, count=nst)
                if len(stars) != nst:
                    raise ValueError(f"Truncated star block: {filename}")
                yield com, stars["pos"], stars["mass"]


def aperture_stellar_mass(filename, radius, box_size, centers=None):
    """Stellar mass within `radius` (code units) of each subhalo centre,
    periodic in `box_size`. Returns an array aligned with GALCATALOG order.

    radius, box_size : same (little-h, comoving) length units as the positions.
    centers          : optional (N, 3) centres; default = each subhalo's COM.
    """
    out = []
    for i, (com, pos, mass) in enumerate(read_galfind_star_particles(filename)):
        if len(mass) == 0:
            out.append(0.0)
            continue
        c = com if centers is None else centers[i]
        dx = pos - c
        dx -= box_size * np.round(dx / box_size)        # periodic minimum image
        inside = (dx * dx).sum(axis=1) < radius * radius
        out.append(float(mass[inside].sum()))
    return np.asarray(out, dtype=np.float64)


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