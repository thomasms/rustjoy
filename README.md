# rustjoy

A Rust port of the **RECONR** and **BROADR** modules of
[NJOY2016](https://github.com/njoy/NJOY2016) (commit `ab27c64`, version
2016.80), with a harness that checks them against the Fortran program byte
for byte.

* RECONR reconstructs pointwise cross sections from ENDF resonance
  parameters and writes a PENDF tape.
* BROADR Doppler-broadens and thins those pointwise cross sections to one or
  more temperatures.

## Results: BROADR

| Check | Cases | Identical |
| --- | --- | --- |
| Every BROADR run in the NJOY2016 test suite, binary and ASCII forms | 108 | **108** |
| Variants: 3 temperatures, bootstrap, thinning only, negative `thnmax`, > 9999 K, 293.6 K thermal quantities, restart, two materials in one call, skipping a material | 51 | **51** |
| Numeric kernels `erfc` and `hnabb` (n = 0..4) | 1,000,152 inputs, 6 values each | **all** |

Fortran line coverage of the 159 cases: `broadr.f90` 1,073 of 1,135
executable lines run (**94.5%**); `erfc` 28 of its 30 lines. What is not run:
error exits, MT152 data with several σ₀ values or with resolved/unresolved
overlap (RECONR does not write either), MF1 sections longer than one page,
and the branches that need more than one page of points (about 10⁶ per
reaction set), among them the `bsigma` truncation warning.

Found during the comparison: the listing line with the broadened MT numbers
(`16i4,(/11x,16i4)`) has an empty line before each group of 16 after the
second, because format reversion ends a record and the `/` ends another.
The port now writes the same.

## Results: RECONR

| Check | Cases | Identical |
| --- | --- | --- |
| Every RECONR run in the NJOY2016 test suite, binary and ASCII forms | 144 | **144** |
| Variants with other options (tolerances, 293.6 K, user grid, cards) | 27 | **27** |
| Synthetic Adler-Adler / hybrid R-function data, against stock NJOY2016 | 6 | 1 (see below) |
| The same 6, against NJOY2016 with the one-line fix `tools/njoy-sigp5-init.patch` | 6 | **6** |
| Numeric kernels (`a11`, `a10`, `sigfig`, `x**n`, `e11.0` read, `E`/`F` output) | 2,009,128 values | **2,009,128** |

"Identical" means: every tape that the run writes, the listing file
`output`, standard output and the exit status are equal byte for byte
(masking only the bytes listed under "How equivalence is checked").

Fortran line coverage of all 177 cases (gcov):

* `reconr.f90`: 3,304 of 3,655 executable lines run (**90.4%**);
* `samm.f90`: 1,809 of 3,979 (45.5%). Of the routines that RECONR can
  reach and that are ported, 57.9% of lines run. The routines that are not
  ported (angular distributions and derivatives, see below) are most of the
  rest.

### Why 5 synthetic cases differ from stock NJOY2016

`csaa` (Adler-Adler) and `cshybr` (hybrid R-function) set only `sigp(1:4)`,
but `sigma` uses `sigp(1:nsig)` with `nsig=5`. In the Fortran, `sigp(5)` is an
uninitialised local array element, so it holds whatever an earlier call left
on the stack (in these runs: 0 for the first 757 calls, then 1.4e7). That
value goes into the convergence test in `resxs`, so the energy grid depends
on stack contents. This is undefined behaviour in NJOY2016 and cannot be
copied; the Rust port uses 0.

`tools/njoy-sigp5-init.patch` changes `do i=1,4` to `do i=1,nsig` in those
two routines. With it, all 6 synthetic cases are identical. No evaluation in
the NJOY2016 test suite uses these two formats (they are from old ENDF/B
versions), which is why synthetic data was needed to test them at all.

### Differences found and fixed during the comparison

* Fe-58 (`NLS=4`, unresolved case A): `unfac` sets nothing for l ≥ 3, so the
  Fortran reuses `vl`/`ps` from the previous (l, J) step, and `vl` is
  multiplied by √E again. That gives a non-physical MT152 table (14,222 b
  elastic at 350 keV), but it is what NJOY2016 writes, so the port now does
  the same. The same reuse is copied for `facts` (unset for l > 4).
* Earlier, before the first commit: CR-LF line ends in one ENDF file, the
  `nrtot` value that `resxs` returns by reference, and `x**n` with a constant
  `n` (gfortran -O0 calls `__powidf2`, it does not expand the product).

## What was ported

| Fortran source | Rust | Notes |
| --- | --- | --- |
| `reconr.f90` (5,747 lines) | `src/reconr.rs` | all routines |
| `samm.f90` (7,169 lines) | `src/samm.rs` | every routine RECONR can reach (LRF=7 R-matrix limited, incl. Coulomb functions and the complex symmetric solver) |
| `endf.f90` | `src/endf.rs` | ENDF tape I/O, `a11`, interpolation |
| `util.f90`, `phys.f90` | `src/util.rs` | `sigfig`, `a10`, scratch buffers, constants |
| `broadr.f90` (1,947 lines) | `src/broadr.rs` | all routines |
| `mathm.f90` | `src/mathm.rs` | `erfc` with `csevl` and `initds` (SLATEC) |
| `main.f90` | `src/main.rs` | banner and module loop (`reconr` and `broadr`) |
| gfortran runtime | `src/fio.rs`, `src/ffmt.rs`, `src/fmath.rs` | sequential units, formatted I/O, `__powidf2`, libm calls |

Not ported (unreachable from RECONR):

* code in `samm.f90` guarded by `Want_Partial_Derivs` or `Want_Angular_Dist`;
  RECONR sets both to false and nothing changes them;
* the SAMMY Reich-Moore (`rdsammy` mode 3) resonance loop: it needs
  `Want_SAMRML_RM`, which is false in RECONR. It stops with an error message
  if it is ever reached.

## Build and run

```bash
cargo build --release
```

The program is used like `njoy`: the input deck comes from standard input,
the listing goes to `output`, tapes are `tapeNN` files:

```bash
./target/release/njoy-reconr < input
```

Decks may contain `reconr`, `broadr` and `--` cards; any other module name
stops with `illegal module name`. For other modules in a chain, run them
with the Fortran `njoy` on the same tapes.

## How equivalence is checked

1. **Reference corpus** (`tools/build_corpus.py`): every RECONR call in the
   NJOY2016 test suite (72 calls) becomes a case, in two forms: the original
   binary-tape form (input made by Fortran `moder`) and an ASCII form that
   reads the ENDF file directly. 144 cases.
2. **Variants** (`tools/make_variants.py`): 3 extra decks for each of 9
   materials chosen to cover SLBW, MLBW, Reich-Moore, R-matrix limited,
   unresolved case A and no-resonance data: looser tolerance; 293.6 K
   (Doppler code, or the "not coded for temperature" error exits); tight
   tolerance with `errmax`/`errint`, 2 descriptive cards and 5 user grid
   points.
3. **Synthetic data** (`tools/make_synthetic.py`): the MF2 of Pu-241 is
   replaced by made-up Adler-Adler (LI = 5, 6, 7; 0 K and Doppler) and
   hybrid R-function data (0 K, and the temperature error exit).
4. **BROADR corpus** (`tools/build_corpus.py --module broadr`): every BROADR
   call in the test suite becomes a case. `pre.input` holds all the modules
   before the call (the Fortran runs them to make the PENDF input), and the
   `_a` form makes the whole chain ASCII. 108 cases.
   `tools/make_variants_broadr.py` adds the 51 option variants.
5. **Comparison** (`tools/compare.py`): each case is run with the Fortran
   `njoy` and with `njoy-reconr` in separate directories. Every tape written,
   the listing `output`, standard output and the exit status are compared
   **byte for byte**. Only these bytes are masked, because they differ between
   two runs of the Fortran program itself:
   * CPU-time fields (`f8.1`) and the date/time in the banner;
   * bytes 68–135 of the binary tape-label record and of each 136-byte
     Hollerith card in binary MF1/MT451 records. The Fortran fills them from
     uninitialised memory (`equivalence` of `character(4)` and `real(kr)`
     arrays). They never appear in ASCII tapes.
6. **Numeric kernels** (`tools/numcheck_gen.py`, `tools/numcheck.f90`,
   `examples/numcheck.rs`): 2,009,128 edge-case and seeded random values go
   through both languages and the outputs are compared byte for byte. For
   BROADR, `tools/numcheck_broadr.f90` and `examples/numcheck_broadr.rs` do
   the same for `erfc` and `hnabb` (`numcheck_gen.py broadr`).
7. **Coverage** (`tools/coverage.sh`): a gcov build of NJOY runs the corpus,
   to show which Fortran lines the comparison exercises.

To rerun everything (builds three Fortran copies, the port, the corpus, and
runs every check; needs gfortran, cmake, cargo and python3):

```bash
tools/run_all.sh /path/to/NJOY2016 /path/to/new-work-dir
```

The scripts never delete anything: they stop if an output directory already
exists.

## Matching gfortran exactly

The Fortran reference is built with NJOY's own CMake defaults (no flags, so
gfortran `-O0`, no FMA, no reordering). Points that the port copies on
purpose:

* `x**n` with an integer `n` (constant or not) calls libgcc `__powidf2`;
  `fmath::powi` is that algorithm. Rust's `powi` is not used.
* `exp`, `log`, `sin`, `atan`, `pow`, `cexp`, `cabs` … are called from the
  system libm, as gfortran does.
* complex multiply/divide use gfortran's Fortran rules (plain product, Smith's
  division);
* a few single-precision literals (`1./(k+1)` in `bigeta`, `float()`) are
  computed in `f32`;
* `SAVE`d locals (e.g. `elim` in `lunion`, `enow` in `csmlbw2`) keep their
  values between materials;
* outputs that a callee leaves unset (`vl`/`ps` from `unfac` for l ≥ 3,
  `se`/`pe` from `facts` for l > 4) keep the caller's previous values;
* by-reference arguments that the Fortran changes (`nrtot`, `e` in `csunr1`,
  `nchan` in `setr`) are changed here too;
* SAMMY arrays are flat column-major memory with the declared shapes, so
  routines that receive only the first slice (`scale3`, `yfour`) or read
  along the wrong index (`backgrdata` in `tab1io`) behave the same;
* CR-LF line ends are record ends, as in libgfortran;
* `initds` adds `abs(real(os(i)))`: each coefficient is rounded to single
  precision first, as the Fortran `real()` does;
* format reversion and trailing `/` produce the same empty records as
  libgfortran (checked against gfortran for the formats used).

## Known limits

* Behaviour that is undefined in the Fortran (reads of uninitialised
  variables or past array ends) cannot be copied. Where the Rust port must
  pick a value, it uses 0. Known places: `sigp(5)` after `csaa`/`cshybr`
  (above), `egridu(79)` in `rdf2u0`, `jj` in `s2sammy` for long LRF=7
  channel lists, `dss` in `genpsf` (l > 4), `igroup` in the `rdsammy` BW
  branch, `eeff` in `cshybr` when there are no competing channels but a
  background R-function or phase table is given, and the first (l, J) step
  when it already has l ≥ 3 (`unfac`) or l > 4 (`facts`).
* BROADR: `thinb` reads its local `tt` before setting it on every page
  after the first, `getunx` can read `loc1` unset, and `mt4i` and `ftev`
  can be read unset. The port keeps `tt` from the previous page and uses 0
  for the others. None of these paths is reached by the test cases (they
  need more than one page of points, or unusual MT152/MT51 data).
* BROADR's `ctev` is reset to 0 for every energy point, so the "thermal
  capture g-factor" line almost never appears. That is NJOY2016's
  behaviour, and the port copies it.
* In `cshybr` the `facts` outputs are not yet kept between calls for l > 4
  (as they are in the other routines); hybrid data with l > 4 is not tested.
* Only Linux x86-64 with glibc was tested. Another libm can change the last
  bit of transcendental functions in both programs differently.
