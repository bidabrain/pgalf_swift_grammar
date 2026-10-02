#!/bin/bash
#==============================================================
# run_pgalf.sh - PGalF pipeline driver for SWIFT single-file HDF5
#   SWIFT snapshot (flamingo_XXXX.hdf5)
#     -> NewDD -> opFoF -> NewGalFinder -> GalCenter
#==============================================================
set -euo pipefail

#--- Defaults (overridable via command-line flags) -----------
PGALF=/home/dbi224/opt/pgalf_swift   # PGalF install dir        (-p)
SIM=""                               # snapshot input dir       (-i)  required
OUT=""                               # output dir               (-o)  required
SNAPSPEC=""                          # snapshot number or range (-s)  required
PREFIX=flamingo                      # snapshot filename prefix (-x)  -> PREFIX_XXXX.hdf5
NSPLIT=64                            # NewDD number of slabs     (-n)
NP=64                                # MPI ranks for gfind       (-N)
MPIRUN="mpirun -np"                  # MPI launcher              (-m)  SLURM: "srun -n"
STAGES=all                           # stages: all | newdd,opfof,gfind,galcenter (-S)
ZMAX=""                              # skip snapshots with z > ZMAX (-z)  empty = no limit

usage() {
    echo "Usage: $0 -i <sim_dir> -o <out_dir> -s <snap|a-b|all> [-n nsplit] [-N np]"
    echo "          [-p pgalf_dir] [-x prefix] [-m \"mpirun -np\"] [-S stages] [-z zmax]"
    exit 1
}

# Read redshift of a SWIFT HDF5 snapshot (needs h5dump on PATH).
# Tries /Header/Redshift, /Cosmology/Redshift, then derives from Scale-factor.
# Echoes the redshift, or empty string if it cannot be determined.
get_redshift() {
    local file="$1" z a
    for attr in /Header/Redshift /Cosmology/Redshift; do
        z=$(h5dump -a "$attr" "$file" 2>/dev/null | awk '/\(0\):/{print $2; exit}')
        [ -n "$z" ] && { echo "$z"; return 0; }
    done
    a=$(h5dump -a /Header/Scale-factor "$file" 2>/dev/null | awk '/\(0\):/{print $2; exit}')
    [ -n "$a" ] && { awk "BEGIN{printf \"%.6f\", 1.0/$a - 1.0}"; return 0; }
    echo ""; return 1
}

#--- Parse arguments -----------------------------------------
while getopts "i:o:s:n:N:p:x:m:S:z:h" opt; do
    case $opt in
        i) SIM=$OPTARG ;;   o) OUT=$OPTARG ;;   s) SNAPSPEC=$OPTARG ;;
        n) NSPLIT=$OPTARG ;; N) NP=$OPTARG ;;   p) PGALF=$OPTARG ;;
        x) PREFIX=$OPTARG ;; m) MPIRUN=$OPTARG ;; S) STAGES=$OPTARG ;;
        z) ZMAX=$OPTARG ;;
        h|*) usage ;;
    esac
done

if [ -z "$SIM" ] || [ -z "$OUT" ] || [ -z "$SNAPSPEC" ]; then
    echo "Error: -i (input) -o (output) -s (snapshot) are required"; usage
fi

# Snapshot selection:
#   "all"   -> every PREFIX_XXXX.hdf5 found in SIM
#   "30-40" -> inclusive range
#   "37"    -> single snapshot
if [ "$SNAPSPEC" = all ]; then
    SNAPS=$(for f in "$SIM"/${PREFIX}_[0-9][0-9][0-9][0-9].hdf5; do
                [ -e "$f" ] || continue
                b=$(basename "$f" .hdf5); echo $((10#${b#${PREFIX}_}))
            done | sort -n -u)
    [ -z "$SNAPS" ] && { echo "No ${PREFIX}_XXXX.hdf5 found in $SIM"; exit 1; }
elif [[ "$SNAPSPEC" == *-* ]]; then
    A=${SNAPSPEC%-*}; B=${SNAPSPEC#*-}
    SNAPS=$(seq "$((10#$A))" "$((10#$B))")
else
    SNAPS=$((10#$SNAPSPEC))
fi

# Whether a given stage should run
want() { [ "$STAGES" = all ] || [[ ",$STAGES," == *",$1,"* ]]; }

mkdir -p "$OUT"
echo "=== PGalF: in=$SIM out=$OUT snaps=$(echo $SNAPS|tr '\n' ' ') nsplit=$NSPLIT np=$NP stages=$STAGES zmax=${ZMAX:-none} ==="

#--- Run the pipeline for each snapshot ----------------------
for SNAP in $SNAPS; do
    S4=$(printf "%04d" "$SNAP")
    SRC="$SIM/${PREFIX}_${S4}.hdf5"
    LINK="$OUT/snap_${S4}.hdf5"

    echo; echo "########## snapshot $SNAP (${S4}) ##########"
    if [ ! -f "$SRC" ]; then echo "!! snapshot not found: $SRC -- skipping"; continue; fi

    # Skip high-redshift snapshots (FoF percolates / crashes before structure forms)
    if [ -n "$ZMAX" ]; then
        Z=$(get_redshift "$SRC")
        if [ -z "$Z" ]; then
            echo "   warning: cannot read redshift from $SRC -- running anyway"
        elif awk "BEGIN{exit !($Z > $ZMAX)}"; then
            echo "   skip: z=$Z > zmax=$ZMAX"; continue
        else
            echo "   z=$Z (<= $ZMAX)"
        fi
    fi

    ln -sf "$SRC" "$LINK"      # symlink to the name NewDD expects: snap_XXXX.hdf5
    cd "$OUT"

    if want newdd;     then echo ">>> [1/4] NewDD";     "$PGALF/NewDD/newdd.exe" "$SNAP" "$NSPLIT"; fi
    if want opfof;     then echo ">>> [2/4] opFoF";      $MPIRUN "$NSPLIT" "$PGALF/opFoF/opfof.exe" "$SNAP" "$NSPLIT"; fi
    if want gfind;     then echo ">>> [3/4] gfind";      $MPIRUN "$NP" "$PGALF/NewGalFinder/gfind.exe" "$SNAP"; fi
    if want galcenter; then echo ">>> [4/4] galcenter"; "$PGALF/GalCenter/galcenter.exe" "$SNAP"; fi

    echo "########## snapshot $SNAP done ##########"
done
echo; echo "All done. Output in $OUT/FoF_Data/"

