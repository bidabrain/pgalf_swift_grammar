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