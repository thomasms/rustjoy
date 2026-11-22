#!/bin/bash
# Run every corpus case with a gcov-instrumented NJOY build and report the
# line coverage of reconr.f90 and samm.f90.
# Usage: coverage.sh BUILD_COV_DIR CORPUS_DIR WORK_DIR
set -e
B=$1; C=$2; W=$3; SRC=$(cd "$B/../src" && pwd)
mkdir -p "$W"
find "$B" -name "*.gcda" -delete
for d in "$C"/*/; do
  n=$(basename "$d"); w="$W/$n"; rm -rf "$w"; mkdir -p "$w"
  cp -L "$d"/tapes/* "$w"/ 2>/dev/null || true
  (cd "$w" && "$B/njoy" < "$d/pre.input" >/dev/null 2>&1; "$B/njoy" < "$d/run.input" >/dev/null 2>&1) || true
  rm -rf "$w"
done
cd "$B"; B=$(pwd); cd "$B/CMakeFiles/njoy_executable.dir/src"
for f in reconr samm; do
  gcov -o . "$SRC/$f.f90" >/dev/null 2>&1
  mv -f "$f.f90.gcov" "$W/" 2>/dev/null || true
  # executable lines: count field not '-' ; covered: not '#####'
  awk -F: '$1 !~ /-/ {t++; if ($1 !~ /#####/) c++} END {printf "%s: %d of %d executable lines run (%.1f%%)\n", FILENAME, c, t, 100*c/t}' "$W/$f.f90.gcov"
done
