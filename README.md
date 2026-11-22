# reconr-rs

A Rust port of the **RECONR** module of [NJOY2016](https://github.com/njoy/NJOY2016)
(commit `ab27c64`, version 2016.80), with a harness that checks it against the
Fortran program byte for byte.

RECONR reconstructs pointwise cross sections from ENDF resonance parameters
and writes a PENDF tape.

## What was ported

| Fortran source | Rust | Notes |
| --- | --- | --- |
| `reconr.f90` (5,747 lines) | `src/reconr.rs` | all routines |
| `samm.f90` (7,169 lines) | `src/samm.rs` | every routine RECONR can reach (LRF=7 R-matrix limited, incl. Coulomb functions and the complex symmetric solver) |
| `endf.f90` | `src/endf.rs` | ENDF tape I/O, `a11`, interpolation |
| `util.f90`, `phys.f90` | `src/util.rs` | `sigfig`, `a10`, scratch buffers, constants |
| `main.f90` | `src/main.rs` | banner and module loop (`reconr` only) |
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

Decks may contain `reconr` and `--` cards; any other module name stops with
`illegal module name`.

## How equivalence is checked

1. **Reference corpus** (`tools/build_corpus.py`): every RECONR call in the
   NJOY2016 test suite (72 calls) becomes a case, in two forms: the original
   binary-tape form (input made by Fortran `moder`) and an ASCII form that
   reads the ENDF file directly. 144 cases.
2. **Variants** (`tools/make_variants.py`): extra decks with other options
   (looser tolerance, 293.6 K, user grid points, descriptive cards).
3. **Comparison** (`tools/compare.py`): each case is run with the Fortran
   `njoy` and with `njoy-reconr` in separate directories. Every tape written,
   the listing `output`, standard output and the exit status are compared
   **byte for byte**. Only these bytes are masked, because they differ between
   two runs of the Fortran program itself:
   * CPU-time fields (`f8.1`) and the date/time in the banner;
   * bytes 68–135 of the binary tape-label record and of each 136-byte
     Hollerith card in binary MF1/MT451 records. The Fortran fills them from
   uninitialised memory (`equivalence` of `character(4)` and `real(kr)`
     arrays). They never appear in ASCII tapes.
4. **Numeric kernels** (`tools/numcheck.f90`, `examples/numcheck.rs`): `a11`,
   `a10`, `sigfig`, `x**n`, `e11.0` reads and `E`/`F` output are run on
   2,009,128 random and edge-case values in both languages and the outputs
   compared byte for byte.
5. **Coverage** (`tools/coverage.sh`): a gcov build of NJOY runs the corpus,
   to show which Fortran lines the comparison exercises.

To rerun everything:

```bash
python3 tools/build_corpus.py /path/to/NJOY2016 corpus
python3 tools/compare.py /path/to/NJOY2016/build/njoy target/release/njoy-reconr corpus work
```

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
* by-reference arguments that the Fortran changes (`nrtot`, `e` in `csunr1`,
  `nchan` in `setr`) are changed here too;
* SAMMY arrays are flat column-major memory with the declared shapes, so
  routines that receive only the first slice (`scale3`, `yfour`) or read
  along the wrong index (`backgrdata` in `tab1io`) behave the same;
* CR-LF line ends are record ends, as in libgfortran.

## Known limits

* Behaviour that is undefined in the Fortran (reads of uninitialised
  variables or past array ends) cannot be copied. Where the Rust port must
  pick a value, it uses 0. Known places: `egridu(79)` in `rdf2u0`, `jj` in
  `s2sammy` for long LRF=7 channel lists, `dss` in `genpsf` (l > 4),
  `igroup` in the `rdsammy` BW branch.
* Only Linux x86-64 with glibc was tested. Another libm can change the last
  bit of transcendental functions in both programs differently.
