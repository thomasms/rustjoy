#!/usr/bin/env python3
"""Make extra GASPR cases with other materials and input options.

Each case runs the Fortran RECONR (and for some cases BROADR, MODER or a
first GASPR) on a neutron ENDF file from the NJOY2016 test resources
(pre.input, loose tolerance to keep it fast) and then GASPR on the PENDF
tape (run.input). GASPR reads the first material of the ENDF tape, so the
first material is used for the multi-material tapes (t404, t511).

  *_g      every neutron file: one temperature, text tapes
  *_t      three temperatures (BROADR)
  *_bin    blocked-binary PENDF tapes (MODER)
  *_again  a second GASPR on a tape that has gas data (old data deleted)
  err_mode text input and binary output (error exit)

Usage: make_variants_gaspr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

MORE = ["J33Pu239", "J33U235", "J33U238", "e6pu241c", "eni61", "t404", "t511"]
SOME = ["n-005_B_010-ENDF8.0.endf", "n-026_Fe_056-JEFF3.3.endf", "n-072_Hf_177-TENDL2021.endf",
        "n-092_U_235-ENDF8.0.endf", "n-082_Pb_208-TENDL2021.endf", "eni61"]


def first_mat(path):
    with open(path, errors="replace") as f:
        f.readline()
        return int(f.readline()[66:70])


def reconr(mat, nout=21):
    return "reconr\n 20 %d/\n 'gaspr variant'/\n %d 0/\n .01/\n 0/\n" % (nout, mat)


def gaspr(nendf, nin, nout):
    return "gaspr\n %d %d %d/\n" % (nendf, nin, nout)


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
            fo.write(pre_in + "stop\n")
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(run_in + "stop\n")
        n += 1

    files = sorted(f for f in os.listdir(res) if f.startswith("n-")) + MORE
    for f in files:
        mat = first_mat(os.path.join(res, f))
        name = f.split("-ENDF")[0].split("-tendl")[0].split("-TENDL")[0].split("-JEFF")[0]
        name = name.split("-IRDFF")[0].split("-beta")[0].replace("-", "_")
        if f.endswith("JEFF3.3.endf"):
            name += "j"
        elif f.endswith("TENDL19.endf"):
            name += "t"
        elif f.endswith("ENDF7.1.endf"):
            name += "7"
        case(name + "_g", f, reconr(mat), gaspr(20, 21, 22))
        if f in SOME:
            case(name + "_t", f, reconr(mat, 22) +
                 "broadr\n 20 22 21/\n %d 3 0 0 0./\n .01/\n 300 900 1800/\n 0/\n" % mat,
                 gaspr(20, 21, 22))
            case(name + "_bin", f, reconr(mat) + "moder\n 21 -23/\n", gaspr(20, -23, -24))
            case(name + "_again", f, reconr(mat, 22) + gaspr(20, 22, 21), gaspr(20, 21, 23))
    f = "n-026_Fe_058-ENDF8.0-Beta4.endf"
    case("err_mode", f, reconr(first_mat(os.path.join(res, f))), gaspr(20, 21, -22))
    print("%d cases in %s" % (n, corpus))


if __name__ == "__main__":
    main()
