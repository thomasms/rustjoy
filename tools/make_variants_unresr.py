#!/usr/bin/env python3
"""Make extra UNRESR cases with other materials and input options.

Each case runs the Fortran RECONR and BROADR on a neutron ENDF file from
the NJOY2016 test resources (pre.input, loose tolerances to keep it fast)
and then UNRESR on the PENDF tape (run.input). The materials cover
unresolved formats LRF=1 and LRF=2, energy-dependent scattering radius,
several isotopes, and materials with no unresolved range (copy path).
For each material there are cases with:

  _a   one temperature, two sigma-zero values, iprint=0
  _b   three temperatures, seven sigma-zero values, iprint=1

Usage: make_variants_unresr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

MATS = [
    # name, ENDF file, mat
    ("pu240", "n-094_Pu_240-ENDF8.0.endf", 9440),
    ("fe58", "n-026_Fe_058-ENDF8.0-Beta4.endf", 2637),
    ("u235", "n-092_U_235-ENDF8.0.endf", 9228),
    ("zr90", "n-040_Zr_090-ENDF8.0.endf", 4025),
    ("mn55", "n-025_Mn_055-ENDF8.0.endf", 2525),
    ("ag109", "n-047_Ag_109-ENDF8.0.endf", 4731),
    ("sr88", "n-038_Sr_088-ENDF8.1.endf", 3837),
    ("am242", "n-095_Am_242-ENDF8.0.endf", 9546),
    ("co58", "n-027_Co_058-ENDF8.0.endf", 2722),
    ("ni63", "n-028_Ni_063-ENDF8.0.endf", 2840),
    ("kr83", "n-036_Kr_083-ENDF8.0.endf", 3640),
    ("ar37", "n-018_Ar_37-tendl2023.endf", 1828),
    ("u238j", "J33U238", 9237),
    ("al27", "n-013_Al_027-ENDF8.0.endf", 1325),
    ("h1", "n-001_H_001-ENDF8.0-Beta6.endf", 125),
]


def pre(mat, temps):
    t = " ".join("%g" % x for x in temps)
    return ("reconr\n 20 21/\n 'unresr variant'/\n %d 0/\n .01/\n 0/\n"
            "broadr\n 20 21 22/\n %d %d 0 0 0./\n .01/\n %s/\n 0/\nstop\n" % (mat, mat, len(temps), t))


def unresr(mat, temps, sigz, iprint):
    t = " ".join("%g" % x for x in temps)
    s = " ".join("%g" % x for x in sigz)
    return ("unresr\n 20 22 23/\n %d %d %d %d/\n %s/\n %s/\n 0/\nstop\n"
            % (mat, len(temps), len(sigz), iprint, t, s))


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    n = 0

    def case(name, f, pre_in, run_in):
        nonlocal n
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        os.symlink(os.path.join(res, f), os.path.join(d, "tapes", "tape20"))
        with open(os.path.join(d, "pre.input"), "w") as fo:
            fo.write(pre_in)
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(run_in)
        n += 1

    for name, f, mat in MATS:
        if not os.path.exists(os.path.join(res, f)):
            continue
        case(name + "_a", f, pre(mat, [293.6]), unresr(mat, [293.6], [1e10, 100.], 0))
        t3 = [300., 900., 2100.]
        case(name + "_b", f, pre(mat, t3),
             unresr(mat, t3, [1e10, 1e5, 1e4, 1000., 100., 10., 1.], 1))
    print("%d cases in %s" % (n, corpus))


if __name__ == "__main__":
    main()
