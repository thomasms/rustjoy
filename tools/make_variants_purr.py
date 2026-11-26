#!/usr/bin/env python3
"""Make extra PURR cases with other materials and input options.

Each case runs the Fortran RECONR and BROADR (and for some cases HEATR,
UNRESR or MODER) on a neutron ENDF file from the NJOY2016 test resources
(pre.input, loose tolerances to keep it fast) and then PURR on the PENDF
tape (run.input). The materials cover unresolved formats LRF=1 (with and
without energy-dependent fission widths) and LRF=2, energy-dependent
scattering radius, LSSF=0 and 1, competition, and materials with no
unresolved range (copy path). For each material there are cases with:

  _a     one temperature, two sigma-zero values, 15 bins, iprint=0, 4 points
  _b     three temperatures, seven sigma-zero values, 20 bins, iprint=1
  _h     HEATR partial heating (302, 318, 402) on the PENDF tape
  _h1    HEATR total heating only (301)
  _u     UNRESR first (the tape has an old MT152)
  _bin   blocked-binary PENDF tapes (MODER)

and a few cases that stop with an error:

  err_nbin   fewer than 15 bins
  err_2mat   a second material (the Fortran stops: array already allocated)
  err_mode   text input and binary output

Usage: make_variants_purr.py NJOY2016_DIR CORPUS_DIR
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


def pre(mat, temps, nout=22):
    t = " ".join("%g" % x for x in temps)
    return ("reconr\n 20 21/\n 'purr variant'/\n %d 0/\n .01/\n 0/\n"
            "broadr\n 20 21 %d/\n %d %d 0 0 0./\n .01/\n %s/\n 0/\n" % (mat, nout, mat, len(temps), t))


def heatr(mat, mts, nin=22, nout=23):
    return "heatr\n 20 %d %d/\n %d %d/\n %s/\n" % (nin, nout, mat, len(mts), " ".join(str(m) for m in mts))


def purr(mat, temps, sigz, nbin, nladr, iprint, nunx, units=(20, 22, 23), more=""):
    t = " ".join("%g" % x for x in temps)
    s = " ".join("%g" % x for x in sigz)
    return ("purr\n %d %d %d/\n %d %d %d %d %d %d %d/\n %s/\n %s/\n%s 0/\nstop\n"
            % (units + (mat, len(temps), len(sigz), nbin, nladr, iprint, nunx, t, s, more)))


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

    s2 = [1e10, 100.]
    s7 = [1e10, 1e5, 1e4, 1000., 100., 10., 1.]
    t1 = [293.6]
    t3 = [300., 900., 2100.]
    for name, f, mat in MATS:
        if not os.path.exists(os.path.join(res, f)):
            continue
        case(name + "_a", f, pre(mat, t1) + "stop\n", purr(mat, t1, s2, 15, 2, 0, 4))
        case(name + "_b", f, pre(mat, t3) + "stop\n", purr(mat, t3, s7, 20, 3, 1, 0))
        case(name + "_h", f, pre(mat, t1) + heatr(mat, [302, 318, 402]) + "stop\n",
             purr(mat, t1, s2, 16, 2, 1, 5, (20, 23, 24)))
        case(name + "_h1", f, pre(mat, t1) + heatr(mat, [444]) + "stop\n",
             purr(mat, t1, s2, 16, 2, 0, 3, (20, 23, 24)))
        case(name + "_u", f, pre(mat, t1) + "unresr\n 20 22 23/\n %d 1 2 0/\n 293.6/\n 1e10 100/\n 0/\nstop\n" % mat,
             purr(mat, t1, s2, 15, 2, 0, 3, (20, 23, 24)))
        case(name + "_bin", f, pre(mat, t1) + "moder\n 22 -23/\nstop\n",
             purr(mat, t1, s2, 15, 2, 0, 3, (20, -23, -24)))
    f, mat = "n-094_Pu_240-ENDF8.0.endf", 9440
    case("err_nbin", f, pre(mat, t1) + "stop\n", purr(mat, t1, s2, 10, 2, 0, 3))
    case("err_2mat", f, pre(mat, t1) + "stop\n",
         purr(mat, t1, s2, 15, 2, 0, 3, more=" %d 1 2 15 2 0 3/\n 293.6/\n 1e10 100/\n" % mat))
    case("err_mode", f, pre(mat, t1) + "stop\n", purr(mat, t1, s2, 15, 2, 0, 3, (20, 22, -23)))
    print("%d cases in %s" % (n, corpus))


if __name__ == "__main__":
    main()
