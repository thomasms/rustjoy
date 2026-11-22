#!/bin/bash
# Rerun the whole equivalence check from a clean NJOY2016 checkout.
#
# Usage: tools/run_all.sh NJOY2016_DIR WORK_DIR
#
# NJOY2016_DIR  a checkout of https://github.com/njoy/NJOY2016 (commit ab27c64)
# WORK_DIR      a new directory; it is created and nothing in it is deleted
#
# Builds, all inside WORK_DIR:
#   njoy-ref   the reference Fortran (NJOY CMake defaults: gfortran -O0)
#   njoy-cov   the same with gcov counters (--coverage)
#   njoy-fix   the reference plus tools/njoy-sigp5-init.patch, used only for
#              the synthetic Adler-Adler / hybrid cases (see README)
set -e
N=$(cd "$1" && pwd); W=$2
R=$(cd "$(dirname "$0")/.." && pwd)
if [ -e "$W" ]; then echo "$W exists; give a new directory" >&2; exit 1; fi
mkdir -p "$W"; W=$(cd "$W" && pwd)

build() {  # build DIR [cmake args...]
  local d=$1; shift
  mkdir -p "$d" && (cd "$d" && cmake "$@" > cmake.log 2>&1 && make -j"$(nproc)" njoy_executable > make.log 2>&1)
}

echo "== build the Fortran references"
build "$W/njoy-ref" "$N"
build "$W/njoy-cov" "$N" -DCMAKE_Fortran_FLAGS=--coverage
mkdir -p "$W/njoy-fix-src"
cp -r "$N/CMakeLists.txt" "$N/src" "$N/metaconfigure" "$W/njoy-fix-src/"
ln -s "$N/tests" "$W/njoy-fix-src/tests"
(cd "$W/njoy-fix-src" && patch -p1 < "$R/tools/njoy-sigp5-init.patch" > /dev/null)
build "$W/njoy-fix" "$W/njoy-fix-src"

echo "== build the Rust port"
(cd "$R" && CARGO_TARGET_DIR="$W/target" cargo build --release -q)
RUST="$W/target/release/njoy-reconr"

echo "== make the corpus"
C="$W/corpus"
python3 "$R/tools/build_corpus.py" "$N" "$C"
python3 "$R/tools/make_variants.py" "$C" 02_1_a 18_1_a 36_1_a 82_2_a 19_1_a 26_1_a 81_1_a 60_1_a 21_1
python3 "$R/tools/make_synthetic.py" "$N" "$C"

echo "== compare Fortran and Rust, byte for byte"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$C" "$W/compare" | tail -1
echo "   synthetic cases against the patched Fortran:"
python3 "$R/tools/compare.py" "$W/njoy-fix/njoy" "$RUST" "$C" "$W/compare-fix" \
  syn_aa7 syn_aa7_t syn_aa5 syn_aa6_t syn_hy syn_hy_t | tail -1

echo "== numeric kernels"
mkdir -p "$W/num"
(cd "$W/num" && gfortran -O0 -o numcheck "$N/src/locale.f90" "$N/src/mainio.f90" \
   "$N/src/util.f90" "$N/src/endf.f90" "$R/tools/numcheck.f90")
(cd "$R" && CARGO_TARGET_DIR="$W/target" cargo build --release -q --example numcheck)
python3 "$R/tools/numcheck_gen.py" > "$W/num/in.bin"
"$W/num/numcheck" < "$W/num/in.bin" > "$W/num/fortran.txt"
"$W/target/release/examples/numcheck" < "$W/num/in.bin" > "$W/num/rust.txt"
if cmp -s "$W/num/fortran.txt" "$W/num/rust.txt"; then
  echo "   identical: $(wc -l < "$W/num/rust.txt") values"
else
  echo "   DIFFER: see $W/num"; exit 1
fi

echo "== Fortran line coverage of the corpus"
bash "$R/tools/coverage.sh" "$W/njoy-cov" "$C" "$W/coverage"
