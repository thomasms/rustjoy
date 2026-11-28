#!/bin/sh
# Run one PREPRO program with the Fortran build and with prustpro in two
# copies of an input directory, then compare every file the runs leave,
# the screen output and the exit status.
#
# usage: tools/prepro/compare.sh PROGRAM INPUT_DIR FORTRAN_BIN_DIR [WORK_DIR]
#
# Only the CPU times ("Total Execution Time ... Seconds") are masked:
# they differ between any two runs.
set -u
prog=$1
inp=$2
fbin=$3
work=${4:-${TMPDIR:-/tmp}/prustpro-compare}
here=$(cd "$(dirname "$0")/../.." && pwd)
rbin=${PRUSTPRO:-$here/target/release/prustpro}

rm -rf "$work/f" "$work/r"
mkdir -p "$work/f" "$work/r"
cp -r "$inp"/. "$work/f/"
cp -r "$inp"/. "$work/r/"

(cd "$work/f" && "$fbin/$prog" > ../f.stdout 2> ../f.stderr; echo $? > ../f.status)
(cd "$work/r" && "$rbin" "$prog" > ../r.stdout 2> ../r.stderr; echo $? > ../r.status)

mask() {
    sed -E 's/(Total Execution Time) +[0-9.*]+ Seconds/\1 X Seconds/' "$1"
}

bad=0
for x in status stdout; do
    if ! mask "$work/f.$x" | cmp -s - "$(mask "$work/r.$x" > "$work/r.$x.m"; echo "$work/r.$x.m")"; then
        echo "DIFF $prog $x"
        bad=1
    fi
done
for f in $(cd "$work/f" && find . -type f | sort); do
    if [ ! -f "$work/r/$f" ]; then
        echo "MISSING $prog $f"
        bad=1
        continue
    fi
    if ! cmp -s "$work/f/$f" "$work/r/$f"; then
        mask "$work/f/$f" > "$work/a.m"
        mask "$work/r/$f" > "$work/b.m"
        if ! cmp -s "$work/a.m" "$work/b.m"; then
            echo "DIFF $prog $f: $(cmp "$work/a.m" "$work/b.m" | head -1)"
            bad=1
        fi
    fi
done
for f in $(cd "$work/r" && find . -type f | sort); do
    [ -f "$work/f/$f" ] || { echo "EXTRA $prog $f"; bad=1; }
done
[ $bad -eq 0 ] && echo "SAME $prog"
exit $bad
