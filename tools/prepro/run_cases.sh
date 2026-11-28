#!/bin/sh
# Run tools/prepro/compare.sh on every case directory of one program.
#
# usage: tools/prepro/run_cases.sh PROGRAM CASES_DIR FORTRAN_BIN_DIR [JOBS]
#
# Prints one line per case (SAME or the differences) and a summary.
set -u
prog=$1
cases=$2
fbin=$3
jobs=${4:-4}
here=$(cd "$(dirname "$0")" && pwd)
work=${TMPDIR:-/tmp}/prustpro-cases-$prog
mkdir -p "$work"
ls "$cases/$prog" | xargs -P "$jobs" -I{} sh -c \
    '"$0/compare.sh" "$1" "$2/$1/$3" "$4" "$5/$3" > "$5/$3.result" 2>&1; rm -rf "$5/$3"; head -5 "$5/$3.result"' \
    "$here" "$prog" "$cases" {} "$fbin" "$work"
same=$(cat "$work"/*.result | grep -c '^SAME')
total=$(ls "$work"/*.result | wc -l)
echo "$prog: $same of $total cases the same"
