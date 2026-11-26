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
CB="$W/corpus-broadr"
python3 "$R/tools/build_corpus.py" "$N" "$CB" --module broadr
CV="$W/corpus-broadr-var"
mkdir -p "$CV"
BB="01_1 36_1_a 19_1 28_1_a 60_1_a 02_1 02_1_a"
for b in $BB; do cp -a "$CB/$b" "$CV/"; done
python3 "$R/tools/make_variants_broadr.py" "$CV" $BB --multi 02_1 1056
CH="$W/corpus-heatr"
python3 "$R/tools/build_corpus.py" "$N" "$CH" --module heatr
CHV="$W/corpus-heatr-var"
python3 "$R/tools/make_variants_heatr.py" "$N" "$CHV"
CM="$W/corpus-moder"
python3 "$R/tools/build_corpus.py" "$N" "$CM" --module moder
CMV="$W/corpus-moder-var"
python3 "$R/tools/make_variants_moder.py" "$N" "$CMV"
CT="$W/corpus-thermr"
python3 "$R/tools/build_corpus.py" "$N" "$CT" --module thermr
CTV="$W/corpus-thermr-var"
python3 "$R/tools/make_variants_thermr.py" "$N" "$CTV"
CU="$W/corpus-unresr"
python3 "$R/tools/build_corpus.py" "$N" "$CU" --module unresr
CUV="$W/corpus-unresr-var"
python3 "$R/tools/make_variants_unresr.py" "$N" "$CUV"

echo "== compare Fortran and Rust, byte for byte"
echo "   reconr:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$C" "$W/compare" | tail -1
echo "   broadr:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CB" "$W/compare-broadr" | tail -1
echo "   broadr variants:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CV" "$W/compare-broadr-var" | tail -1
echo "   heatr:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CH" "$W/compare-heatr" | tail -1
echo "   heatr variants:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CHV" "$W/compare-heatr-var" | tail -1
echo "   moder:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CM" "$W/compare-moder" | tail -1
echo "   moder variants:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CMV" "$W/compare-moder-var" | tail -1
echo "   thermr:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CT" "$W/compare-thermr" | tail -1
echo "   thermr variants:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CTV" "$W/compare-thermr-var" | tail -1
echo "   unresr:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CU" "$W/compare-unresr" | tail -1
echo "   unresr variants:"
python3 "$R/tools/compare.py" "$W/njoy-ref/njoy" "$RUST" "$CUV" "$W/compare-unresr-var" | tail -1
echo "   synthetic cases against the patched Fortran:"
python3 "$R/tools/compare.py" "$W/njoy-fix/njoy" "$RUST" "$C" "$W/compare-fix" \
  syn_aa7 syn_aa7_t syn_aa5 syn_aa6_t syn_hy syn_hy_t | tail -1

echo "== the parallel prefetch: 1 thread against the default ($(nproc) CPUs), byte for byte"
# compare.py runs its first program in the role of the reference; here that
# is the port with NJOY_THREADS=1 (no worker threads at all). The modules
# before each call still run with the Fortran program (COMPARE_SETUP).
printf '#!/bin/sh\nNJOY_THREADS=1 exec "%s"\n' "$RUST" > "$W/rust-1thread"
chmod +x "$W/rust-1thread"
echo "   reconr:"
COMPARE_SETUP="$W/njoy-ref/njoy" python3 "$R/tools/compare.py" "$W/rust-1thread" "$RUST" "$C" "$W/threads" | tail -1
echo "   broadr:"
COMPARE_SETUP="$W/njoy-ref/njoy" python3 "$R/tools/compare.py" "$W/rust-1thread" "$RUST" "$CB" "$W/threads-broadr" | tail -1
echo "   broadr variants:"
COMPARE_SETUP="$W/njoy-ref/njoy" python3 "$R/tools/compare.py" "$W/rust-1thread" "$RUST" "$CV" "$W/threads-broadr-var" | tail -1

echo "== numeric kernels"
echo "   allocation-free formatting and reading against the previous code:"
(cd "$R" && CARGO_TARGET_DIR="$W/target" cargo build --release -q --example ffmt_check)
"$W/target/release/examples/ffmt_check" | sed 's/^/   /'
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
sed 's/^   public broadr$/   public broadr, hnabb/' "$N/src/broadr.f90" > "$W/num/broadr_pub.f90"
(cd "$W/num" && gfortran -O0 -o numcheck_broadr "$N/src/locale.f90" "$N/src/mainio.f90" \
   "$N/src/phys.f90" "$N/src/util.f90" "$N/src/endf.f90" "$N/src/mathm.f90" \
   broadr_pub.f90 "$R/tools/numcheck_broadr.f90" 2> /dev/null)
(cd "$R" && CARGO_TARGET_DIR="$W/target" cargo build --release -q --example numcheck_broadr)
python3 "$R/tools/numcheck_gen.py" broadr > "$W/num/in_broadr.bin"
"$W/num/numcheck_broadr" < "$W/num/in_broadr.bin" > "$W/num/fortran_broadr.txt"
"$W/target/release/examples/numcheck_broadr" < "$W/num/in_broadr.bin" > "$W/num/rust_broadr.txt"
if cmp -s "$W/num/fortran_broadr.txt" "$W/num/rust_broadr.txt"; then
  echo "   erfc/hnabb identical: $(wc -l < "$W/num/rust_broadr.txt") inputs (6 values each)"
else
  echo "   erfc/hnabb DIFFER: see $W/num"; exit 1
fi

echo "== Fortran line coverage of the corpus"
bash "$R/tools/coverage.sh" "$W/njoy-cov" "$C" "$W/coverage"
CA="$W/corpus-broadr-all"
mkdir -p "$CA"
for d in "$CB"/* "$CV"/*; do [ -e "$CA/$(basename "$d")" ] || ln -s "$d" "$CA/$(basename "$d")"; done
bash "$R/tools/coverage.sh" "$W/njoy-cov" "$CA" "$W/coverage-broadr" broadr mathm
CHA="$W/corpus-heatr-all"
mkdir -p "$CHA"
for d in "$CH"/* "$CHV"/*; do [ -e "$CHA/$(basename "$d")" ] || ln -s "$d" "$CHA/$(basename "$d")"; done
bash "$R/tools/coverage.sh" "$W/njoy-cov" "$CHA" "$W/coverage-heatr" heatr
CMA="$W/corpus-moder-all"
mkdir -p "$CMA"
for d in "$CM"/* "$CMV"/*; do [ -e "$CMA/$(basename "$d")" ] || ln -s "$d" "$CMA/$(basename "$d")"; done
bash "$R/tools/coverage.sh" "$W/njoy-cov" "$CMA" "$W/coverage-moder" moder
CTA="$W/corpus-thermr-all"
mkdir -p "$CTA"
for d in "$CT"/* "$CTV"/*; do [ -e "$CTA/$(basename "$d")" ] || ln -s "$d" "$CTA/$(basename "$d")"; done
bash "$R/tools/coverage.sh" "$W/njoy-cov" "$CTA" "$W/coverage-thermr" thermr
CUA="$W/corpus-unresr-all"
mkdir -p "$CUA"
for d in "$CU"/* "$CUV"/*; do [ -e "$CUA/$(basename "$d")" ] || ln -s "$d" "$CUA/$(basename "$d")"; done
bash "$R/tools/coverage.sh" "$W/njoy-cov" "$CUA" "$W/coverage-unresr" unresr
