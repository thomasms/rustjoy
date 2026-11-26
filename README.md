# rustjoy

A Rust port of the **MODER**, **RECONR**, **BROADR**, **HEATR**, **THERMR**,
**UNRESR**, **GASPR** and **PURR** modules of
[NJOY2016](https://github.com/njoy/NJOY2016) (commit `ab27c64`, version
2016.80), with a harness that checks them against the Fortran program byte
for byte.

* MODER converts ENDF-type tapes between text and blocked-binary form, and
  selects materials from one or more tapes onto a new tape.
* RECONR reconstructs pointwise cross sections from ENDF resonance
  parameters and writes a PENDF tape.
* BROADR Doppler-broadens and thins those pointwise cross sections to one or
  more temperatures.
* HEATR computes heating (kerma) and radiation damage energy production
  from the PENDF and ENDF data and adds them to the PENDF tape (MT301 and
  the requested partial kermas).
* THERMR computes thermal neutron scattering cross sections and
  energy-angle distributions (free gas, S(α,β) tables, coherent and
  incoherent elastic) and adds them to the PENDF tape.
* UNRESR computes self-shielded effective cross sections in the unresolved
  resonance range (Bondarenko method) and adds them to the PENDF tape as
  MF2/MT152.
* GASPR adds gas production cross sections (MT203-207: protons,
  deuterons, tritons, helium-3 and alphas) to the PENDF tape.
* PURR samples resonance ladders in the unresolved range (Monte Carlo) to
  make probability tables (MF2/MT153) and Bondarenko cross sections
  (MF2/MT152) and adds them to the PENDF tape.

## Speed

The port gives the same bytes as before, and is faster:

Wall-clock seconds, best of 2 runs (Fortran: 1 run), on one 4-CPU Linux
machine. "Before" is the previous commit of this port.

| Case | Fortran (gfortran -O0) | Port before | Port now, 1 thread | Port now, 4 threads |
| --- | ---: | ---: | ---: | ---: |
| RECONR U-235 ENDF/B-VIII.0 (test 65), binary tapes | 68.0 | 40.0 | 14.5 | **5.8** |
| RECONR same, ASCII tapes | 101.2 | 53.3 | 22.4 | **15.1** |
| RECONR U-235 JENDL-3.3 (test 17) | 81.7 | 45.7 | 14.3 | **6.1** |
| RECONR Pu-239 (test 27), ASCII | 16.3 | 8.9 | 4.4 | **2.2** |
| BROADR U-235 (test 65), binary tapes | 26.4 | 13.5 | 13.4 | **5.6** |
| BROADR same, ASCII tapes | 50.0 | 28.4 | 21.7 | **13.3** |
| BROADR U-235 JENDL-3.3 (test 17) | 19.1 | 7.8 | 7.1 | **3.2** |
| BROADR test 15, ASCII | 15.4 | 6.9 | 6.1 | **2.6** |

Whole corpus (sum over all cases, from `run_all.sh`; 4 threads):

| Corpus | Fortran | Port before | Port now |
| --- | ---: | ---: | ---: |
| RECONR, 177 cases | 576 s | 363 s | **86 s** |
| BROADR, 108 cases | 355 s | 186 s | **75 s** |
| BROADR variants, 58 cases | 34 s | 17 s | **11 s** |

The ASCII cases gain less: reading and writing the text tapes runs on the
main thread.

What changed (none of it changes a floating-point operation or its order):

* **Reich-Moore resonances (`csrmat`).** For each energy, the Fortran scans
  the full resonance list once per J value and channel sign and skips the
  resonances that do not match. The port now builds, once per section, the
  list of resonances that each scan keeps (in the same order) and loops over
  only those. The sums are the same additions in the same order.
* **Fewer allocations:** `sigma` reuses its work arrays (zeroed each call,
  as a new allocation would be); Fortran `I`, `F`, `E` editing and the
  ENDF number reading build their text on the stack. The old routines are
  kept, and `examples/ffmt_check.rs` checks the new ones against them on
  110,998,901 random and edge-case inputs (all identical); `numcheck`
  (against gfortran) is unchanged.
* **Parallel prefetch** (`src/prefetch.rs`). RECONR `resxs` and BROADR
  `broadn` refine an energy grid panel by panel. Worker threads now refine
  the panels ahead of the main program, each with its own copy of the
  data, and keep every cross section they compute. The main program then
  runs the original sequential code unchanged; when it needs a cross
  section, it takes the worker's value if one exists for exactly that
  energy, else it computes it itself. A worker's value is used only if the
  call had no side effect: it wrote nothing, did not stop, and did not use
  state from an earlier call (for example the `csrmat` singular-matrix
  case). State that a call leaves behind is rebuilt before the next real
  call. BROADR's "broadening truncated" message is repeated from the stored
  value. So the output cannot depend on the thread count or timing; the
  check in `run_all.sh` compares 1 thread with the default, byte for byte,
  for every case. The prefetch is used for SLBW, MLBW (0 K), Reich-Moore and
  materials without resolved resonances; the other formats (MLBW with
  temperature, which keeps state between calls; Adler-Adler; hybrid;
  R-matrix limited through SAMMY) run as before.

The number of threads is the number of CPUs; set `NJOY_THREADS=1` to turn
the prefetch off, or another number to choose it. The machine used for the
times above has 4 virtual CPUs that give about 2.7 CPUs of real throughput,
so a machine with more cores gains more.

## Results: THERMR

| Check | Cases | Identical |
| --- | --- | --- |
| Every THERMR run in the NJOY2016 test suite, binary and ASCII forms | 30 | **30** |
| Variants: S(α,β) data from 5 thermal files (H in H2O, Al, D in LiD, H and Zr in ZrH), each in E-μ-E' and E-E'-μ form and with the full listing; free gas at several temperatures, bin counts and print levels; a temperature above 3000 K; built-in coherent elastic (graphite, Be, BeO) and incoherent elastic (icoh 11, 12, 13); 4 error exits | 33 | **33** |

Fortran line coverage of the 63 cases: `thermr.f90` 1,912 of 2,076
executable lines run (**92.1%**). What is not run: most error exits, the
continuation records of long thermal-file records, the second effective
temperature of mixed-moderator data (`gatef2`), and some limit branches of
`calcem` and `sig`.

THERMR also needed `acos` from libm, added to `src/fmath.rs`.

## Results: UNRESR

| Check | Cases | Identical |
| --- | --- | --- |
| Every UNRESR run in the NJOY2016 test suite, binary and ASCII forms | 10 | **10** |
| Variants: 15 materials (Pu-240, Fe-58, U-235, Zr-90, Mn-55, Ag-109, Sr-88, Am-242, Co-58, Ni-63, Kr-83, Ar-37, JENDL-3.3 U-238, and Al-27 and H-1 with no unresolved range), each with one temperature and two sigma-zero values, and with three temperatures, seven sigma-zero values and the full listing | 30 | **30** |

Fortran line coverage of the 40 cases: `unresr.f90` 985 of 1,044
executable lines run (**94.3%**). What is not run: error exits, the
continuation records of long lists, an energy-dependent scattering radius
inside the unresolved range (`nro=1`), and some branches of the `uw`
complex probability function.

## Results: GASPR

| Check | Cases | Identical |
| --- | --- | --- |
| Every GASPR run in the NJOY2016 test suite, binary and ASCII forms | 9 | **9** |
| Variants: every neutron ENDF file of the test resources and 6 more evaluations (JENDL-3.3, ENDF/B-IV and ENDF/B-V tapes), 28 of them with MF6/MT5 yields; for 6 materials also three temperatures, blocked-binary tapes, and a second GASPR on a tape that has gas data; 1 error exit | 60 | **60** |

Fortran line coverage of the 69 cases: `gaspr.f90` 790 of 955 executable
lines run (**82.7%**). Most lines not run are entries of the reaction
tables (MT152-200 and rare LR values) that no test file has; a script
wrote those Rust lines directly from the Fortran. Also not run: MF6 laws 6
and 7, and the gas production listing (`iprint` is fixed at 0).

## Results: PURR

| Check | Cases | Identical |
| --- | --- | --- |
| Every PURR run in the NJOY2016 test suite, binary and ASCII forms | 27 | **27** |
| Variants: 15 materials (as for UNRESR), each with six option sets: one or three temperatures, two or seven sigma-zero values, 15 to 20 bins, both print options, partial or total heating from HEATR, an old MT152 from UNRESR, and blocked-binary tapes; 3 error exits (too few bins, a second material, text input with binary output) | 93 | **93** |

Fortran line coverage of the 120 cases: `purr.f90` 1,729 of 1,823
executable lines run (**94.8%**). What is not run: error exits, an
energy-dependent scattering radius (`nro=1`), the fluctuation integral
with three widths, and code that fixed values make unreachable (`nmode`
is 0 and `nsamp` is 10000 in the Fortran). The random-number sequence
(`rann`) and the quicksort are ported as they are, so the sampled ladders
are the same.

## Results: MODER

| Check | Cases | Identical |
| --- | --- | --- |
| Every MODER run in the NJOY2016 test suite, binary and ASCII forms | 286 | **286** |
| Variants: every ENDF file of the test resources text→binary and binary→text (neutron, photonuclear, charged-particle, thermal scattering, photo-atomic, atomic relaxation, decay and covariance data, ENDF/B-III to VIII), and material selection from several text or binary tapes, a blank tape id, and the wrong-order error | 142 | **142** |

One more case of the test suite cannot be compared: in the all-ASCII form of
test 10 (`10_2_a`) a later Fortran module of the chain stops, so there is no
input tape for the MODER call (the binary form `10_2` is compared).

Fortran line coverage of the 428 cases: `moder.f90` 753 of 968 executable
lines run (**77.8%**). What is not run: error exits, the continuation
records of record types that no test file has longer than one page, and a
few rare formats (File 7 MT451, File 15 LF=2, tabulated File 14, some
File 32 forms). These use the same record routines as the rest.

MODER also needed `intgio` (INTG records of File 32), added to
`src/endf.rs`.

## Results: HEATR

| Check | Cases | Identical |
| --- | --- | --- |
| Every HEATR run in the NJOY2016 test suite, binary and ASCII forms | 27 | **27** |
| Variants: 46 materials from the test resources (H to Pu, light-ion incident data, ENDF/B-III to VIII, JEFF, JENDL, TENDL), each with 2 of 5 option sets: 14 partial kermas, the kinematic check with plot file, user and energy-dependent Q values, local photon deposition with a given damage energy, defaults | 92 | **92** |

Of the 92 variants, 10 stop with an NJOY error (for example data that HEATR
does not support); the error message, the tapes written up to that point and
the exit status are identical as well.

Fortran line coverage of the 119 cases: `heatr.f90` 3,614 of 3,994
executable lines run (**90.5%**); the test suite alone runs 81.9%. What is
not run: error exits, File 6 discrete lines (`h6dis` and the delta-function
part of `h6cm`), tabulated cm angular distributions in File 6, polynomial
ν̄ and MT455 data, recoil from lab-system File 6 data, and ENDF-5 charged
particle levels without angular data.

Found during the comparison:

* `nsyso` is unit 0 (`ERROR_UNIT`, connected to `output`), so the
  `tab1io(0,nsyso,...)` that should list the user's energy-dependent Q
  table writes nothing: only the line ` input q ----` appears.
* An `X` skip at the end of a record writes nothing in libgfortran (the
  `low`/`high` marker lines of the kinematic check).

## Results: BROADR

| Check | Cases | Identical |
| --- | --- | --- |
| Every BROADR run in the NJOY2016 test suite, binary and ASCII forms | 108 | **108** |
| Variants: 3 temperatures, bootstrap, thinning only, negative `thnmax`, > 9999 K, 293.6 K thermal quantities, restart, two materials in one call, skipping a material | 51 | **51** |
| Numeric kernels `erfc` and `hnabb` (n = 0..4) | 1,000,152 inputs, 6 values each | **all** |
| Port with 1 thread against 4 threads (all BROADR cases and variants) | 166 | **166** |

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
| Port with 1 thread against 4 threads (every RECONR case) | 177 | **177** |
| New allocation-free `I`/`F`/`E` editing and number reading against the previous code | 110,998,901 cases | **all** |

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
| `moder.f90` (1,714 lines) | `src/moder.rs` | all routines |
| `heatr.f90` (6,322 lines) | `src/heatr.rs` | all routines |
| `thermr.f90` (3,387 lines) | `src/thermr.rs` | all routines |
| `unresr.f90` (1,665 lines) | `src/unresr.rs` | all routines |
| `gaspr.f90` (1,153 lines) | `src/gaspr.rs` | all routines |
| `purr.f90` (2,919 lines) | `src/purr.rs` | all routines |
| `mathm.f90` | `src/mathm.rs` | `erfc` with `csevl` and `initds` (SLATEC), `legndr` |
| `endf.f90` (`skip6`, `intgio`) | `src/endf.rs` | added for HEATR and MODER |
| `main.f90` | `src/main.rs` | banner and module loop (`moder`, `reconr`, `broadr`, `heatr`, `thermr`, `unresr`, `gaspr`, `purr`) |
| — | `src/prefetch.rs` | worker threads for the parallel prefetch (no Fortran equivalent) |
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

Decks may contain `moder`, `reconr`, `broadr`, `heatr`, `thermr`, `unresr`,
`gaspr`, `purr` and `--` cards; any other module name
stops with `illegal module name`. For other modules in a chain, run them
with the Fortran `njoy` on the same tapes.

`NJOY_THREADS=n` sets the number of worker threads (default: the number of
CPUs; 1 turns them off). The output is the same for any value.

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
   **HEATR corpus**: the same for every HEATR call (27 cases), and
   `tools/make_variants_heatr.py` adds 92 variants (RECONR by the Fortran,
   then HEATR with one of five option sets).
   **MODER corpus**: the same for every MODER call (286 cases), and
   `tools/make_variants_moder.py` adds 142 variants (every resource file in
   both directions, and material selection).
   **THERMR corpus**: the same for every THERMR call (30 cases), and
   `tools/make_variants_thermr.py` adds 33 variants (RECONR and BROADR by
   the Fortran, then THERMR with S(α,β), free-gas or built-in data).
   **UNRESR corpus**: the same for every UNRESR call (10 cases), and
   `tools/make_variants_unresr.py` adds 30 variants (RECONR and BROADR by
   the Fortran, then UNRESR on 15 materials).
   **GASPR corpus**: the same for every GASPR call (9 cases), and
   `tools/make_variants_gaspr.py` adds 60 variants (RECONR and other
   modules by the Fortran, then GASPR).
   **PURR corpus**: the same for every PURR call (27 cases), and
   `tools/make_variants_purr.py` adds 93 variants (RECONR, BROADR and
   sometimes HEATR, UNRESR or MODER by the Fortran, then PURR).
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
7. **Threads**: `compare.py` also compares the port with `NJOY_THREADS=1`
   against the port with the default thread count, for every case.
8. **Coverage** (`tools/coverage.sh`): a gcov build of NJOY runs the corpus,
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
* HEATR: with a plot file (`nplot`) but without `iprint=2`, `hout` plots the
  minimum and maximum kinematic columns, which exist only with `iprint=2`:
  it reads elements of its local `c(30)` that are never set (the Fortran
  writes NaN or other leftover values). The port has 0 there, so these plot
  files differ; the variants use the plot file only with `iprint=2`. The
  test suite never asks for this.
* HEATR: `conbar` sets `t1`/`t2` (discrete-energy File 5 law, LF=3) only
  when it initializes a reaction, but they are not `save`d; the port keeps
  the values from that call. No test data uses LF=3.
* BROADR's `ctev` is reset to 0 for every energy point, so the "thermal
  capture g-factor" line almost never appears. That is NJOY2016's
  behaviour, and the port copies it.
* In `cshybr` the `facts` outputs are not yet kept between calls for l > 4
  (as they are in the other routines); hybrid data with l > 4 is not tested.
* Only Linux x86-64 with glibc was tested. Another libm can change the last
  bit of transcendental functions in both programs differently.
