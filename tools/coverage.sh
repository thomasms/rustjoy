#!/bin/bash
# Run every corpus case with a gcov-instrumented NJOY build and report the
# line coverage of reconr.f90 and samm.f90.
# Usage: coverage.sh BUILD_COV_DIR CORPUS_DIR WORK_DIR
# WORK_DIR must not exist yet. Nothing is deleted: the counters (.gcda) are
# written under WORK_DIR/gcda through GCOV_PREFIX, not next to the build.
set -e
B=$(cd "$1" && pwd); C=$(cd "$2" && pwd); W=$3
SRC=$(cd "$B/../src" && pwd)
OBJ="$B/CMakeFiles/njoy_executable.dir/src"
if [ -e "$W" ]; then echo "$W exists; give a new directory" >&2; exit 1; fi
mkdir -p "$W/runs" "$W/gcda"
W=$(cd "$W" && pwd)
export GCOV_PREFIX="$W/gcda" GCOV_PREFIX_STRIP=0
for d in "$C"/*/; do
  n=$(basename "$d"); w="$W/runs/$n"; mkdir -p "$w"
  for t in "$d"/tapes/*; do [ -e "$t" ] && ln -s "$(readlink -f "$t")" "$w/$(basename "$t")"; done
  (cd "$w" && "$B/njoy" < "$d/pre.input" >/dev/null 2>&1; "$B/njoy" < "$d/run.input" >/dev/null 2>&1) || true
done
G="$W/gcda$OBJ"
for f in reconr samm; do
  cp "$OBJ/$f.f90.gcno" "$G/"
  (cd "$W" && gcov -o "$G" "$SRC/$f.f90" >/dev/null)
  # executable lines: count field not '-'; covered: not '#####'
  awk -F: '$1 !~ /-/ {t++; if ($1 !~ /#####/) c++} END {printf "%s: %d of %d executable lines run (%.1f%%)\n", FILENAME, c, t, 100*c/t}' "$W/$f.f90.gcov"
done
