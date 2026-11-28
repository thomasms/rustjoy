# rustjoy

A Rust port of the **MODER**, **RECONR**, **BROADR**, **HEATR**, **THERMR**,
**UNRESR**, **GASPR**, **PURR**, **LEAPR**, **MIXR**, **RESXSR**, **DTFR**,
**GAMINR**, **GROUPR** and **POWR** modules of [NJOY2016](https://github.com/njoy/NJOY2016) (commit `ab27c64`,
version 2016.80), with a harness that checks them against the Fortran
program byte for byte.

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
* LEAPR computes the thermal scattering law S(α,β) from a phonon frequency
  spectrum (phonon expansion, with free-gas or diffusive translation,
  discrete oscillators, cold hydrogen and deuterium, and Skold
  coherence), with incoherent or coherent (Bragg edge) elastic data, and
  writes it as an ENDF-6 File 7 tape. It needs no input tape.
* MIXR makes a PENDF tape of mixed cross sections: each output MT is a
  weighted sum of that MT from several materials, on the union energy
  grid (for example, elemental cross sections).
* RESXSR writes a RESXS file (binary, CCCC-like): the elastic, fission
  and capture cross sections of several materials and temperatures on
  one thinned energy grid, for hyper-fine flux calculations.
* DTFR converts multigroup data from a GENDF tape (GROUPR or GAMINR
  output) into DTF-format transport tables (text or CLAW format), and can
  write plot commands for PLOTR (group values and PENDF curves, and
  isometric plots of the scattering and photon production tables).
* GAMINR makes multigroup photon interaction cross sections, scattering
  matrices and heating kerma factors from ENDF photon data (MF23 cross
  sections and MF27 form factors) and a PENDF tape, and writes them on a
  GENDF tape.
* GROUPR makes multigroup cross sections from a PENDF tape: self-shielded
  (Bondarenko) cross sections and ratios, neutron scattering matrices
  (Files 4, 5 and 6, thermal MF6 data), fission and delayed-neutron
  spectra, photon and charged-particle production matrices, and nuclide
  production, with 36 neutron and 10 photon group structures, 12 weight
  functions and a flux calculator. It writes a GENDF tape.
* POWR formats multigroup data from a GENDF tape as input for the
  EPRI-CELL codes GAMTAP (fast library, 68 groups, lib=1) and LIBRAR
  (thermal library, 35 groups, lib=2), and for the EPRI-CPM code CLIB
  (69 groups, lib=3). A 185-group GENDF tape is collapsed to these
  structures.

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
  kept, and `tools/checks/ffmt_check.rs` checks the new ones against them on
  130,997,261 random and edge-case inputs (all identical); `numcheck`
  (against gfortran) is unchanged.
* **Parallel prefetch** (`src/rustjoy/prefetch.rs`). RECONR `resxs` and BROADR
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

### Second round: PURR, THERMR, LEAPR and ENDF output

Wall-clock seconds, best of 2 runs, on the same 4-CPU machine. "Before" is
commit `68b70b1` (POWR port); "now" is this version, with 4 threads. The
tapes are identical, byte for byte; only the times in the listing change.

| Case | Before | Now | Speed-up |
| --- | ---: | ---: | ---: |
| `examples/pu241-unresr` (UNRESR, HEATR, PURR) | 15.9 | **4.3** | 3.7x |
| `examples/zrh-thermr` (THERMR, ZrH) | 6.9 | **1.9** | 3.6x |
| THERMR test 25 (3 temperatures) | 73.8 | **7.1** | 10.4x |
| LEAPR test 80 | 43.0 | **16.7** | 2.6x |
| `examples/sn119-purr` | 0.6 | **0.4** | 1.5x |

What changed (again, no floating-point operation or its order changes):

* **PURR.** The ladders of resonances are still drawn one after the other,
  so the random numbers come in the same order. Then the energy points are
  divided into blocks, one block per thread. Each energy point gets the same
  additions in the same order as before (ladder, then temperature, then
  resonance). The cross-section table is now calculated when it is needed,
  and the search in it is the same as in the Fortran.
* **THERMR.** The search in the alpha and beta grids is a binary search when
  the grid increases (else the original loop). In `calcem`, worker threads
  calculate the energies ahead of the main program, for both the
  incoherent (`iform` 0) and the equi-probable angle (`iform` 1) forms. The
  only value that goes from one energy to the next (`jnz`) is returned by
  each worker and checked by the main program. The main program uses a
  worker's result only if the worker wrote nothing and did not stop; else
  it calculates the energy itself.
* **LEAPR (`trans`).** The Fortran can read values that an earlier beta or
  alpha left in the `sb` work array, when `sbfill` fills fewer points than
  the convolution reads. The port now finds first, for each beta, the step
  and the number of points that `sbfill` uses. Then worker threads
  calculate the betas, and each one rebuilds exactly the values that the
  Fortran reads, old values included. At the end of each alpha, the main
  `sb` array is rebuilt to its Fortran state. Small cases (fewer than
  200,000 points) use the original serial code.
* **ENDF output.** The 11-column numbers are now made from exact integer
  arithmetic, with no text on the heap, and each 80-column line is built in
  one buffer. `ffmt_check` compares the new code with the old one on 20
  million extra cases.

The check in `run_all.sh` compares 1 thread with 4 threads, byte for byte,
for every THERMR, PURR and LEAPR case, and for their variants.

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

THERMR also needed `acos` from libm, added to `src/rustjoy/fmath.rs`.

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

## Results: LEAPR

| Check | Cases | Identical |
| --- | --- | --- |
| Every LEAPR run in the NJOY2016 test suite (tests 09, 22, 23, 33 twice, 80), and the ASCII form of test 09 | 7 | **7** |
| Variants: H in H2O with free gas, diffusion and oscillators, the three print levels, several temperatures (with reuse of the parameters), `lat=1`, `isabt=1`, `ilog=1`, a pure solid (incoherent elastic), phonon order 260, free, diffusive and mixed secondary scatterers, 320 α values, large α values; coherent elastic for every lattice (`iel` 1 to 6), with one, two and three temperatures; ortho and para hydrogen and deuterium (`ncold` 1 to 4) with each `isabt`/`ilog` pair, for one and two temperatures; both s(κ) options (Vineyard, Skold); 2 error exits | 66 | **66** |

Fortran line coverage of the 73 cases: `leapr.f90` 1,990 of 2,038
executable lines run (**97.6%**). What is not run: error exits, code that
fixed values make unreachable (the moment check in `stable` has
`icheck=0`, the free-gas branch in `coldh` has `ifree=0`, and the extra
temperature points in `endout` need `ntf > ntempr`), a zero S(α,β) in
`terps`, large arguments in `sjbes`, and one branch of the hexagonal
Bragg-edge search.

Test 80 (200 α by 1,325 β values with diffusion) takes 154 s with the
Fortran and 58 s with the port.

## Results: MIXR

| Check | Cases | Identical |
| --- | --- | --- |
| No NJOY2016 test calls MIXR. Variants: two to ten materials and up to twenty MTs on PENDF tapes; an MT missing from one tape or from all tapes; 100,000-point sections (many pages); raw ENDF tapes (ENDF/B-IV, -V and -VI formats, incident protons, deuterons, helium-3 and alphas); the second and third temperature of BROADR tapes; blocked-binary tapes; a negative weight; 2 error exits | 28 | **28** |

Fortran line coverage of the 28 cases: `mixr.f90` 307 of 309 executable
lines run (**99.4%**). The two lines not run are error exits that the
code cannot reach (the array limit and an uninitialised `gety`).

## Results: RESXSR

| Check | Cases | Identical |
| --- | --- | --- |
| No NJOY2016 test calls RESXSR. Variants: U-235 with 6 energy ranges, thinning tolerances (0 to 0.1) and 0 to 10 text cards; Fe-56 (no fission); U-238; one to three temperatures of a BROADR tape (`maxt` 1, 2, 3 and 9); two materials on two tapes; a blocked-binary PENDF tape; 1e-5 eV to 20 MeV (many records); two RESXSR runs in one deck; 2 error exits | 18 | **18** |

Fortran line coverage of the 18 cases: `resxsr.f90` 211 of 212
executable lines run (**99.5%**). The line not run sets the ENDF-5 format
flag; no test file gives an ENDF-5 PENDF tape.

## Results: DTFR

| Check | Cases | Identical |
| --- | --- | --- |
| The DTFR run of the NJOY2016 test suite (test 03, on GAMINR output) | 1 | **1** |
| Variants on GROUPR output of tests 01, 02, 07, 08 and 30 and of U-235 with delayed neutrons (MT455, MF5): edits, thermal tables (MT229/230), upscatter, Legendre orders, photon production tables, CLAW format, plot commands (none, 1 and 4 plots per frame, plots in the listing with `nplot=0`), up to 11 MTs in one edit position, temperature and sigma-zero choice, a material that is not on the tape, binary tapes, two materials in one run; U photoatomic data (ENDF-6) through GAMINR; 5 error exits | 32 | **32** |

Fortran line coverage of the 33 cases: `dtfr.f90` 902 of 912 executable
lines run (**98.9%**). What is not run: the three storage-limit errors,
the MT19-after-MT18 and delayed-nubar errors, a GENDF header longer than
one page, and the zero padding of a CLAW photon table (the number of
neutron by photon groups is always a multiple of 6 here).

## Results: GAMINR

| Check | Cases | Identical |
| --- | --- | --- |
| Every GAMINR run in the NJOY2016 test suite (tests 03 and 30), binary and ASCII forms | 3 | **3** |
| Variants: all 10 group structures (one read-in structure with 320 groups, so that records are longer than one page, and one that ends below the pair threshold), the three weight options, Legendre orders 0 to 8, both print options, explicit reaction lists, no output tape, binary output, two GAMINR runs in one deck, an input GAMINR tape (material lower, equal and higher), H, Fe and U together (the incoherent matrix of the lowest Z is used again for higher Z), H and U photoatomic data (ENDF-6, MT522/525); 5 error exits and one end-of-file stop | 33 | **33** |

Fortran line coverage of the 36 cases: `gaminr.f90` 806 of 814 executable
lines run (**99.0%**). What is not run: the ENDF-5 format flag (no ENDF-5
photon file), the form-factor storage error, one exit of the incoherent
panel loop, and code that is not reachable (6-point quadrature in
`gpanel`, more than one page of total heating values for one group).

## Results: GROUPR

| Check | Cases | Identical |
| --- | --- | --- |
| Every GROUPR run that `build_corpus.py` finds in the NJOY2016 test suite (20 tests), binary and ASCII forms | 48 | **48** |
| Variants: all 35 built-in neutron group structures and all 10 gamma structures (two read in), all 12 weight functions, the flux calculator (iwt<0, with admixed and external moderators) and a flux read back from `ninwt` (iwt=0), Legendre orders 0 to 8, both print and smoothing options, two temperatures, sigma zeroes with and without unresolved data, the automatic lists (`3/`, `6/`, `10/`, `13/`, `16/`, `21/` to `26/`, `31/`, `35/`), extended nuclide production cards (`-1/`), MT251-253 and MT257-259, fission spectra, nubar and delayed neutrons (U-235 ENDF/B-V and the test 76 deck), File 6 laws 1, 2, 3, 6 and 7 and the discrete capture gamma of H-1, particle and nuclide production (Be, H-2, B-10, Al, Ca, Fe, Co, Ni, Pu-241), thermal MF6 data from THERMR (H in ZrH, Al), incident protons, deuterons and alphas, incident photons, gout tapes (replace a material, add MTs, a material below the first one on the tape, two materials); 12 error exits, and 8 runs that stop part way with the same error in both programs (data that GROUPR cannot handle, a wrong MT range, and the gout-tape options of Known limits) | 115 | **115** |

Fortran line coverage of the 163 cases: `groupr.f90` 5,571 of 6,236 executable
lines run (**89.3%**). What is not run: most error exits; the Madland-Nix
spectrum (not ported, see above); discrete lines in centre-of-mass File 6
data (`f6dis`) and tabulated angular distributions in File 6 law 1
(LANG=11-15); the Coulomb amplitude expansion of charged-particle elastic
data in `getco`; ENDF-5 charged-particle MT numbers; the temperature
search of `getfwt` for a second temperature; the recoil branch of
`gam102`; and the unresolved sigma-zero warnings.

## Results: POWR

| Check | Cases | Identical |
| --- | --- | --- |
| No NJOY2016 test uses POWR. Variants: lib=1 (carbon, U-235 with fission and delayed neutrons, Be-9, Pu-238 with unresolved self-shielding; one and two temperatures, up to three sigma zeroes, f-factors on and off, a reference temperature that is not the first, the reference sigma zero and sigma-zero limit, two materials, two POWR runs in one deck, 185-group collapse, ASCII GENDF, negative `nout`), lib=2 (carbon, U-235, H-1 free gas from THERMR; transport correction, two temperatures and two materials, 185-group collapse, two runs), lib=3 (burnup data, modes 1 and 2, 185-group collapse, two materials, a negative mat, a lib=3 run after a lib=1 run); 15 error exits | 44 | **44** |
| lib=3 against NJOY2016 with `tools/njoy-powr-cards.patch`, the port built with `--features powr-cards` (see below): subfiles 2 to 5 and the CLIB output, with materials from the tape (nina 0 and 1) and from the input (nina 2 and 3), modes 0, 1 and 2, burnup data, one and two temperatures, a reference sigma zero, no resonance fission; one of the 13 stops with the same runtime error in both programs | 13 | **13** |

Fortran line coverage: the 44 cases run 1,568 of 3,021 executable lines of
`powr.f90` (51.9%); the unpatched program cannot get past subfile 2 of
lib=3 (see below). Counting also the lines that the 13 patched cases run
(`tools/gcov_union.py` maps them back through the patch), 2,658 of the
3,021 lines (**88.0%**) are run by a compared case. What is not run:
`sfile6` and `sf6out` after their first records and subfile 6 of the CLIB
output (`sfile6` never finds a P1 matrix, see Known limits), the
absorption-only input of lib=1 (unreachable, see Known limits), and
storage and input errors.

### Why POWR lib=3 is also checked against a patched NJOY2016

In stock NJOY2016, lib=3 reads card 5 of each material into local
variables but then takes the values from the scratch array `scr` (zero in
the first POWR run), and it never stores the resonance lambdas. Every
material is then treated as `nina=0` with sigma zero 0, and `sfile2` stops
("reference sigma zero not found"), so subfiles 2 to 6 cannot be compared
with real input. `tools/njoy-powr-cards.patch` makes three changes for the
comparison only: it stores the card-5 values and the lambdas where the
code reads them back; it allocates `tempt` in `sfile4` also when there are
resonance materials (otherwise an unallocated array is written); and it
changes the continuation format of the P0 matrix print in `sf4out` from
`(9x,i4,2x,1p,5e12.4)` to `(14x,1p,5e12.4)` (the `i4` gets a real number,
which stops the program for every matrix row longer than 5 groups). The
Rust feature `powr-cards` makes the same three changes; the default build
has none of them and matches stock NJOY2016.

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
`src/rustjoy/endf.rs`.

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
| `reconr.f90` (5,747 lines) | `src/rustjoy/reconr.rs` | all routines |
| `samm.f90` (7,169 lines) | `src/rustjoy/samm.rs` | every routine RECONR can reach (LRF=7 R-matrix limited, incl. Coulomb functions and the complex symmetric solver) |
| `endf.f90` | `src/rustjoy/endf.rs` | ENDF tape I/O, `a11`, interpolation |
| `util.f90`, `phys.f90` | `src/rustjoy/util.rs` | `sigfig`, `a10`, scratch buffers, constants |
| `broadr.f90` (1,947 lines) | `src/rustjoy/broadr.rs` | all routines |
| `moder.f90` (1,714 lines) | `src/rustjoy/moder.rs` | all routines |
| `heatr.f90` (6,322 lines) | `src/rustjoy/heatr.rs` | all routines |
| `thermr.f90` (3,387 lines) | `src/rustjoy/thermr.rs` | all routines |
| `unresr.f90` (1,665 lines) | `src/rustjoy/unresr.rs` | all routines |
| `gaspr.f90` (1,153 lines) | `src/rustjoy/gaspr.rs` | all routines |
| `purr.f90` (2,919 lines) | `src/rustjoy/purr.rs` | all routines |
| `leapr.f90` (3,625 lines) | `src/rustjoy/leapr.rs` | all routines |
| `mixr.f90` (533 lines) | `src/rustjoy/mixr.rs` | all routines |
| `resxsr.f90` (505 lines) | `src/rustjoy/resxsr.rs` | all routines |
| `dtfr.f90` (1,510 lines) | `src/rustjoy/dtfr.rs` | all routines |
| `gaminr.f90` (1,517 lines) | `src/rustjoy/gaminr.rs` | all routines |
| `groupr.f90` (12,690 lines) | `src/rustjoy/groupr/` (`mod.rs`, `weight.rs`, `core.rs`, `mf6.rs`, `dist.rs`, `conver.rs`) and `src/rustjoy/groupr_tables.rs` (group structures and weight functions, made from the Fortran by a script) | all routines; the Madland-Nix spectrum (MF5 LF=12 in `anased`, which needs `e1` and `gami`) stops with a message |
| `powr.f90` (4,902 lines) | `src/rustjoy/powr/` (`mod.rs`: main routine and the fast library; `therm.rs`; `cpm.rs` and `cpm4.rs`: the CPM library) | all routines |
| `mathm.f90` | `src/rustjoy/mathm.rs` | `erfc` with `csevl` and `initds` (SLATEC), `legndr` |
| `endf.f90` (`skip6`, `intgio`) | `src/rustjoy/endf.rs` | added for HEATR and MODER |
| `main.f90` | `src/rustjoy/main.rs` | banner and module loop (`moder`, `reconr`, `broadr`, `heatr`, `thermr`, `unresr`, `gaspr`, `purr`, `leapr`, `mixr`, `resxsr`, `dtfr`, `gaminr`, `groupr`, `powr`) |
| — | `src/rustjoy/prefetch.rs` | worker threads for the parallel prefetch (no Fortran equivalent) |
| gfortran runtime | `src/rustjoy/fio.rs`, `src/common/ffmt.rs`, `src/rustjoy/fmath.rs` | sequential units, formatted I/O, `__powidf2`, libm calls |

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
./target/release/rustjoy < input
```

Decks may contain `moder`, `reconr`, `broadr`, `heatr`, `thermr`, `unresr`,
`gaspr`, `purr`, `leapr`, `mixr`, `resxsr`, `dtfr`, `gaminr`, `groupr`, `powr` and `--` cards; any other module name
stops with `illegal module name`. For other modules in a chain, run them
with the Fortran `njoy` on the same tapes.

`NJOY_THREADS=n` sets the number of worker threads (default: the number of
CPUs; 1 turns them off). The output is the same for any value.

[`examples/`](examples) has 11 input decks with their ENDF files, which
use all 14 ported modules, and tells how to run them.

## Tests

`tests/examples.rs` runs each deck in `examples/` with rustjoy and
compares the results with NJOY2016 byte for byte: the exit status,
standard output, the listing `output` and every tape that the deck writes.

```bash
cargo test --release
```

This needs no Fortran: it compares with the NJOY2016 results stored in
`examples/*/reference.txt` (a SHA-256 checksum of each output file). The
references were made with NJOY2016 commit `ab27c64`, built with its CMake
defaults (gfortran -O0) on Linux x86-64. On another platform, a different
libm can change the last bit of some results.

To run NJOY2016 on the same decks and compare the files directly (a
difference is then shown with the first line that differs), give the path
to an `njoy` program:

```bash
NJOY=/path/to/njoy cargo test --release --test examples
```

The debug build (`cargo test`) gives the same results, but it is about
six times slower. After a test passes, its work directory in `target/tmp` is
removed; after a failure, it is kept for inspection.

## How equivalence is checked

The example tests above are a small, quick check. The full check below
runs thousands of cases made from the NJOY2016 test suite.

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
   **LEAPR corpus**: the same for every LEAPR call (7 cases), and
   `tools/make_variants_leapr.py` adds 66 decks (LEAPR alone, no input
   tapes).
   **MIXR corpus**: no test calls MIXR, so
   `tools/make_variants_mixr.py` makes all 28 cases (RECONR and
   sometimes BROADR or MODER by the Fortran, then MIXR).
   **RESXSR corpus**: no test calls RESXSR either, so
   `tools/make_variants_resxsr.py` makes all 18 cases (RECONR and
   sometimes BROADR or MODER by the Fortran, then RESXSR).
   **DTFR corpus**: the same for the one DTFR call (test 03), and
   `tools/make_variants_dtfr.py` adds 32 variants (the module chain of a
   test up to its GROUPR run, or RECONR, BROADR, GROUPR or GAMINR, by the
   Fortran, then DTFR).
   **GAMINR corpus**: the GAMINR calls of tests 03 and 30, and
   `tools/make_variants_gaminr.py` adds 33 variants (RECONR, and for the
   input-tape cases a first GAMINR run, by the Fortran, then GAMINR).
   **GROUPR corpus**: the GROUPR calls of the test suite (the module
   chain up to each call is run by the Fortran), and
   `tools/make_variants_groupr.py` adds 115 variants (MODER, RECONR,
   BROADR and sometimes UNRESR, THERMR or a first GROUPR run by the
   Fortran, then GROUPR).
   **POWR corpus**: no test calls POWR, so `tools/make_variants_powr.py`
   makes all 44 cases (MODER, RECONR, BROADR and sometimes UNRESR or
   THERMR, then GROUPR for the 68-, 35-, 69- and 185-group structures, by
   the Fortran, then POWR). With `--cards` it makes the 13 lib=3 cases
   for the patched Fortran and the `powr-cards` build.
5. **Comparison** (`tools/compare.py`): each case is run with the Fortran
   `njoy` and with `rustjoy` in separate directories. Every tape written,
   the listing `output`, standard output and the exit status are compared
   **byte for byte**. Only these bytes are masked, because they differ between
   two runs of the Fortran program itself:
   * CPU-time fields (`f8.1`) and the date/time in the banner;
   * bytes 68–135 of the binary tape-label record and of each 136-byte
     Hollerith card in binary MF1/MT451 records. The Fortran fills them from
     uninitialised memory (`equivalence` of `character(4)` and `real(kr)`
     arrays). They never appear in ASCII tapes.
6. **Numeric kernels** (`tools/numcheck_gen.py`, `tools/numcheck.f90`,
   `tools/checks/numcheck.rs`): 2,009,128 edge-case and seeded random values go
   through both languages and the outputs are compared byte for byte. For
   BROADR, `tools/numcheck_broadr.f90` and `tools/checks/numcheck_broadr.rs` do
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
* LEAPR, coherent elastic output: when the number of Bragg edges kept
  (`jmax`) is a multiple of `npage/2` (153 for the default page), the
  first temperature fills the TAB1 record exactly, and the next edges are
  written to `scr(-1)` and `scr(0)`, before the start of the array. The
  Fortran then aborts (`munmap_chunk(): invalid pointer`) at the
  `deallocate` statements after `endout`, and loses the end of its listing
  and of the tape (the MEND and TEND records). The port keeps those two
  words (the Fortran reads one back as the running sum) and finishes
  normally: its tape is the same as the Fortran tape for every line the
  Fortran wrote. One variant hit this with BeO; it now uses graphite.
* LEAPR, `isym=3` with `ilog=1`: the Fortran tests `ssp(-nbeta+1,j,nt)`
  (a value nbeta words before `ssp(1,j,nt)`). The port reads the same
  word; for j=1 at the first temperature that word is before the array,
  and the port takes it as positive.
* MIXR reads each input tape forward only, as the Fortran does. An MT
  that comes after a larger MT in the input list is reported as "not
  present". When an MT is missing from a tape, the Fortran goes back to
  the start of MF3 and then skips one section, so the first MF3 section
  (normally MT1) of that tape is not found for the next MT. The port does
  the same.
* RESXSR reads the descriptive text cards but does not write them: the
  Fortran fills the set-hollerith record with blanks. The port does the
  same.
* RESXSR with an ENDF/B-IV tape sets the number of words per scratch
  record to `n2h` of the second MF1 card. For the test file t404 that is
  0, and the Fortran stops with a floating-point exception (division by
  zero in `loada`). The port stops with a panic at the same point. A
  RECONR PENDF tape of t404 makes the Fortran stop earlier, with a bad
  integer read; the port gives the same output there.
* RESXSR thinning keeps up to 50 points in `stack(50,31)` with no check.
  If more than 50 points in a row are inside the tolerance, the Fortran
  writes past the array. The port uses a larger array, so it can differ
  from the Fortran in that case. No case here keeps more than 50.
* DTFR with `nplot=0` writes its plot commands to unit 0, which is the
  listing (`nsyso` is `ERROR_UNIT`, connected to `output`). The port
  writes them to the listing too.
* GAMINR skips every record of material 1 on an input GAMINR tape
  (`ngam1`), because it takes `math=1` as the tape header. It then reads
  to the end of the file and stops with a runtime error. The port does the
  same (case `in_mat1`).
* GROUPR's options to add MTs to a material on an input gout tape do not
  work in NJOY2016: `matb<0` on card 2 is not made positive (findf fails);
  `matb<0` on card 10 skips the old material on the input tape and then
  reads past the end of it (runtime error); `matb<0` with no input gout
  tape takes the tape label as the first material. The port does the same
  (cases `out_add1`, `out_add`, `out_all`).
* GROUPR with an incident charged particle and a request for MF3 MT2
  (for example the automatic `3/` list) makes the Fortran `getmf6` clear
  one column of the feed-function array for each group, but that array
  has only 2 words for MF3. This writes past the array, corrupts the heap,
  and the Fortran stops with `free(): invalid pointer`. The port does not
  copy this; no case here uses it.
* GROUPR with incident photons (photonuclear data) above the top neutron
  group: the `jgmax` search in `getmf6` reads past the end of the group
  bounds. The port stops the search at the end of the array. The value is
  used only for an upscatter correction that does not apply to these
  cases.
* GROUPR reads three locals that are not set on every path: `ehigh` in
  `panel` (a jump to label 115 from the lower boundary), `lcd` in the
  normal entry of `getgfl`, and `x1`, `x2` in the Kalbach branch of
  `f6ddx` outside the tabulated range. The port uses the value of the last
  call, 1, and the last value. No case here shows a difference.
* GROUPR's Madland-Nix fission spectrum (MF5 LF=12) is not ported (no
  test file has it); the port stops with a message there.
* POWR lib=1: a negative `matd` (absorption read from the input) ends the
  material loop, because the loop test `matd.gt.0` comes before the sign
  is used. The absorption-only code is ported but never runs, in both
  programs.
* POWR lib=3 takes the card-5 values from `scr` (see Results: POWR). In a
  later POWR run of the same deck they are what the earlier run left
  there; the port keeps `scr` between runs, so case `c_after_fast` is
  identical. After a lib=2 run those values send `pinit` past its
  material search and the Fortran prints its group-bound array without
  setting it; that is not copied (no case uses it).
* POWR lib=3, undefined behaviour that is not copied: mode 0 with
  position 0 (from `scr`) uses `iflgs(0)` and `idnto(0)`, outside the
  arrays (the port keeps one word for `iflgs(0)` and uses 0 for
  `idnto(0)`); `indf` is printed but set only when a material has
  `ires=1`; the group bounds are unset when the first material is not on
  the GENDF tape; `sfile3` reads `flxr` and `flux` before their start when
  the first temperature has a zero flux (the port uses 0); `sfile4` writes
  the unallocated `tempt` when there are resonance materials (the Fortran
  crashes, the port panics). No compared case reaches these.
* POWR runtime errors that are copied (same message, exit status 2, and,
  as in libgfortran, the record being written is lost): the P0 matrix
  print in `sf4out` for rows longer than 5 groups, the `i6` edit of a real
  in the burnup input (`if5=2`), a second lib=3 run in one deck (its
  arrays are still allocated), and lib=3 with no `nina` 0 or 3 material and
  no `ires=1` (deallocation of `glamt`, which was not allocated).
* POWR `sfile6` jumps to `tosend` after the first record of every MF6
  section, so it never finds a P1 matrix and stops with "no required p1
  matrices found" for a material with `nina=0` and `ip1opt=0`; `sfile3`
  writes resonance integrals only for `nina>0`, so CLIB output with a
  `nina=0`, `ires=1` material fails with a runtime error. The port does
  the same.
* BROADR's `ctev` is reset to 0 for every energy point, so the "thermal
  capture g-factor" line almost never appears. That is NJOY2016's
  behaviour, and the port copies it.
* In `cshybr` the `facts` outputs are not yet kept between calls for l > 4
  (as they are in the other routines); hybrid data with l > 4 is not tested.
* Only Linux x86-64 with glibc was tested. Another libm can change the last
  bit of transcendental functions in both programs differently.

## License

rustjoy is released under the [MIT license](LICENSE). It is a modified
version of NJOY2016 and is not endorsed by LANL; the parts derived from
NJOY2016 also keep the [NJOY2016 license](LICENSE-NJOY2016).
