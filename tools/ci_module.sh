#!/bin/bash
# The equivalence check of tools/run_all.sh for one module, for CI (one job
# per module) or to run by hand. It stops with status 1 if any case differs,
# other than the known cases listed below.
#
# Usage: tools/ci_module.sh MODULE NJOY2016_DIR WORK_DIR
#
# MODULE        reconr, broadr, heatr, moder, thermr, unresr, gaspr, purr,
#               leapr, mixr, resxsr, dtfr, gaminr, groupr or powr
# NJOY2016_DIR  a checkout of https://github.com/njoy/NJOY2016 (commit ab27c64)
# WORK_DIR      a new directory; it is created and nothing in it is deleted
#
# Set NJOY_REF to an njoy program built with the NJOY CMake defaults to use
# it instead of building one in WORK_DIR.
set -e
M=$1
N=$(cd "$2" && pwd); W=$3
R=$(cd "$(dirname "$0")/.." && pwd)
if [ -e "$W" ]; then echo "$W exists; give a new directory" >&2; exit 1; fi
mkdir -p "$W"; W=$(cd "$W" && pwd)

build() {  # build DIR SRC [cmake args...]
  local d=$1; shift
  mkdir -p "$d" && (cd "$d" && cmake "$@" > cmake.log 2>&1 && make -j"$(nproc)" njoy_executable > make.log 2>&1) ||
    { echo "Fortran build failed in $d:" >&2; tail -30 "$d"/*.log >&2; exit 1; }
}
patched() {  # patched NAME PATCH: the reference with a patch, in WORK_DIR/NAME
  mkdir -p "$W/$1-src"
  cp -r "$N/CMakeLists.txt" "$N/src" "$N/metaconfigure" "$W/$1-src/"
  ln -s "$N/tests" "$W/$1-src/tests"
  (cd "$W/$1-src" && patch -p1 < "$R/tools/$2" > /dev/null)
  build "$W/$1" "$W/$1-src"
}

echo "== build"
if [ -n "$NJOY_REF" ]; then F=$NJOY_REF; else build "$W/njoy-ref" "$N"; F="$W/njoy-ref/njoy"; fi
(cd "$R" && CARGO_TARGET_DIR="$W/target" cargo build --release -q)
RUST="$W/target/release/rustjoy"
printf '#!/bin/sh\nNJOY_THREADS=1 exec "%s"\n' "$RUST" > "$W/rust-1thread"
chmod +x "$W/rust-1thread"

FAIL=0
# check NAME PROGRAM_A PROGRAM_B CORPUS [KNOWN...]: compare.py, then fail on
# any case that is not identical and is not in KNOWN
check() {
  local name=$1 a=$2 b=$3 corpus=$4; shift 4
  python3 "$R/tools/compare.py" "$a" "$b" "$corpus" "$W/c-$name" | tail -1 | sed "s/^/   $name: /"
  python3 - "$W/c-$name/results.json" "$@" <<'EOF' || FAIL=1
import json, sys
res = json.load(open(sys.argv[1])); known = set(sys.argv[2:])
bad = [r["case"] for r in res if r["status"] != "identical" and r["case"] not in known]
gone = sorted(known - {r["case"] for r in res if r["status"] != "identical"})
if gone:
    print("   note: known cases now identical: " + " ".join(gone))
if bad:
    print("   DIFFERS: " + " ".join(bad))
    sys.exit(1)
EOF
}
# threads NAME CORPUS: 1 thread against the default, byte for byte
threads() {
  COMPARE_SETUP="$F" check "threads-$1" "$W/rust-1thread" "$RUST" "$2"
}

echo "== corpus and compare: $M"
C="$W/corpus"
case $M in
  reconr)
    python3 "$R/tools/build_corpus.py" "$N" "$C"
    python3 "$R/tools/make_variants.py" "$C" 02_1_a 18_1_a 36_1_a 82_2_a 19_1_a 26_1_a 81_1_a 60_1_a 21_1
    python3 "$R/tools/make_synthetic.py" "$N" "$C"
    # these synthetic cases differ from stock NJOY2016 (see README); they
    # are compared with the patched Fortran below
    SYN="syn_aa7 syn_aa7_t syn_aa5 syn_aa6_t syn_hy syn_hy_t"
    check reconr "$F" "$RUST" "$C" $SYN
    patched njoy-fix njoy-sigp5-init.patch
    python3 "$R/tools/compare.py" "$W/njoy-fix/njoy" "$RUST" "$C" "$W/c-fix" $SYN | tail -1 | sed "s/^/   synthetic, patched Fortran: /"
    python3 -c "import json,sys; sys.exit(any(r['status']!='identical' for r in json.load(open('$W/c-fix/results.json'))))" || { echo "   DIFFERS"; FAIL=1; }
    threads reconr "$C"
    ;;
  broadr)
    python3 "$R/tools/build_corpus.py" "$N" "$C" --module broadr
    BB="01_1 36_1_a 19_1 28_1_a 60_1_a 02_1 02_1_a"
    mkdir -p "$W/corpus-var"
    for b in $BB; do cp -a "$C/$b" "$W/corpus-var/"; done
    python3 "$R/tools/make_variants_broadr.py" "$W/corpus-var" $BB --multi 02_1 1056
    check broadr "$F" "$RUST" "$C"
    check broadr-var "$F" "$RUST" "$W/corpus-var"
    threads broadr "$C"
    threads broadr-var "$W/corpus-var"
    ;;
  leapr)
    python3 "$R/tools/build_corpus.py" "$N" "$C" --module leapr
    python3 "$R/tools/make_variants_leapr.py" "$W/corpus-var"
    check leapr "$F" "$RUST" "$C"
    check leapr-var "$F" "$RUST" "$W/corpus-var"
    threads leapr "$C"
    threads leapr-var "$W/corpus-var"
    ;;
  dtfr)
    python3 "$R/tools/build_corpus.py" "$N" "$C" --module dtfr
    python3 "$R/tools/make_variants_dtfr.py" "$N" "$W/dtfr-work" "$W/corpus-var"
    check dtfr "$F" "$RUST" "$C"
    check dtfr-var "$F" "$RUST" "$W/corpus-var"
    ;;
  mixr|resxsr)
    python3 "$R/tools/make_variants_$M.py" "$N" "$W/corpus-var"
    check "$M-var" "$F" "$RUST" "$W/corpus-var"
    ;;
  powr)
    python3 "$R/tools/make_variants_powr.py" "$N" "$W/corpus-var"
    python3 "$R/tools/make_variants_powr.py" --cards "$N" "$W/corpus-cards"
    check powr-var "$F" "$RUST" "$W/corpus-var"
    patched njoy-pfix njoy-powr-cards.patch
    (cd "$R" && CARGO_TARGET_DIR="$W/target-cards" cargo build --release -q --features powr-cards)
    check powr-cards "$W/njoy-pfix/njoy" "$W/target-cards/release/rustjoy" "$W/corpus-cards"
    ;;
  heatr|moder|thermr|unresr|gaspr|purr|gaminr|groupr)
    python3 "$R/tools/build_corpus.py" "$N" "$C" --module "$M"
    python3 "$R/tools/make_variants_$M.py" "$N" "$W/corpus-var"
    # moder 10_2_a: a later Fortran module of test 10 stops, so there is
    # nothing to compare (see README)
    KNOWN=""; [ "$M" = moder ] && KNOWN="10_2_a"
    check "$M" "$F" "$RUST" "$C" $KNOWN
    check "$M-var" "$F" "$RUST" "$W/corpus-var"
    case $M in thermr|purr) threads "$M" "$C"; threads "$M-var" "$W/corpus-var";; heatr|groupr|gaminr) threads "$M" "$C";; esac
    ;;
  *) echo "unknown module: $M" >&2; exit 2;;
esac

if [ $FAIL -ne 0 ]; then echo "== $M: FAILED"; exit 1; fi
echo "== $M: all identical"
