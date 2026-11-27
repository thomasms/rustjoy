#!/usr/bin/env python3
"""Make RESXSR cases. No NJOY2016 test deck calls resxsr, so all cases
are here.

Each case links ENDF files from the NJOY2016 test resources (tape20,
tape30, ...), runs the Fortran RECONR (and for some cases BROADR or
MODER) in pre.input to make PENDF tapes, and then RESXSR in run.input.
The RESXS file is a binary tape.

  u5_*      U-235 (fissile): energy ranges, thinning tolerances, text cards
  fe        Fe-56 (no fission: two reactions)
  temps_*   BROADR tapes with three temperatures; maxt 1, 2, 3 and 9
  two       two materials on two tapes, with different temperatures
  bin       a blocked-binary PENDF tape (MODER)
  big       1e-5 eV to 20 MeV (many records)
  err_endf4 the RECONR PENDF tape of the ENDF/B-IV file t404 (the Fortran
            stops with a bad integer read)
  again     two RESXSR runs in one deck, on two copies of the PENDF tape
  err_again two RESXSR runs on the same unit (the Fortran does not rewind
            it and stops at the end of the file)

Usage: make_variants_resxsr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

UNITS = [20, 30, 40, 50]


def first_mat(path):
    with open(path, errors="replace") as f:
        f.readline()
        return int(f.readline()[66:70])


def reconr(nin, nout, mats, err=".005"):
    s = "reconr\n %d %d/\n 'resxsr variant'/\n" % (nin, nout)
    for m in mats:
        s += " %d 0/\n %s/\n" % (m, err)
    return s + " 0/\n"


def broadr(nendf, nin, nout, mat, temps):
    return "broadr\n %d %d %d/\n %d %d 0 0 0./\n .005/\n %s/\n 0/\n" % (
        nendf, nin, nout, mat, len(temps), " ".join(str(t) for t in temps))


def resxsr(nout, mats, maxt, efirst, elast, eps, holl=("resxsr variant",), uid="njoy user",
           ivers=1):
    s = "resxsr\n %d/\n %d %d %d %s %s %s/\n '%s' %d/\n" % (
        nout, len(mats), maxt, len(holl), efirst, elast, eps, uid, ivers)
    for h in holl:
        s += " '%s'/\n" % h
    for h, m, u in mats:
        s += " '%s' %d %d/\n" % (h, m, u)
    return s


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    n = 0

    def mat(f):
        return first_mat(os.path.join(res, f))

    def case(name, files, pre_in, run_in):
        nonlocal n
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        for u, f in zip(UNITS, files):
            os.symlink(os.path.join(res, f), os.path.join(d, "tapes", "tape%d" % u))
        with open(os.path.join(d, "pre.input"), "w") as fo:
            fo.write(pre_in + "stop\n")
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(run_in + "stop\n")
        n += 1

    u5 = "n-092_U_235-ENDF8.0.endf"
    u8 = "J33U238"
    fe = "n-026_Fe_056-JEFF3.3.endf"
    pu = "n-094_Pu_239-ENDF8.0-Beta6.endf"
    m5, m8, mfe, mpu = mat(u5), mat(u8), mat(fe), mat(pu)

    pre = reconr(20, 21, [m5])
    for name, ef, el, eps, holl in [
            ("u5_1", 1, 100, ".001", ("one text card",)),
            ("u5_2", 10, 2000, ".005", ("first card", "second card", "third card")),
            ("u5_3", ".5", 50, "0.", ()),
            ("u5_4", 1, 1000, ".01", ("x" * 72,)),
            ("u5_5", 2, 30, ".1", ("a",)),
            ("u5_6", 100, 2250, ".002", ("b",) * 10)]:
        case(name, [u5], pre, resxsr(31, [("u235", m5, 21)], 1, ef, el, eps, holl,
                                     uid="user %s" % name, ivers=len(holl)))
    pre = reconr(20, 21, [mfe])
    case("fe", [fe], pre, resxsr(31, [("fe56", mfe, 21)], 1, 1000, 100000, ".001"))

    pre = reconr(20, 22, [m5]) + broadr(20, 22, 21, m5, [300, 900, 2400])
    for maxt in (1, 2, 3, 9):
        case("temps_%d" % maxt, [u5], pre,
             resxsr(31, [("u235", m5, 21)], maxt, 1, 500, ".002"))

    pre = (reconr(20, 22, [m5]) + broadr(20, 22, 21, m5, [300, 600]) +
           reconr(30, 32, [mpu]) + broadr(30, 32, 31, mpu, [300, 900]))
    case("two", [u5, pu], pre,
         resxsr(41, [("u235", m5, 21), ("pu239xyz", mpu, 31)], 2, 1, 300, ".001",
                ("two materials",), uid="abcdefghijkl", ivers=7))

    pre = reconr(20, 21, [m5]) + "moder\n 21 -23/\n"
    case("bin", [u5], pre, resxsr(31, [("u235", m5, -23)], 1, 1, 200, ".001"))

    pre = reconr(20, 21, [m5], ".01")
    case("big", [u5], pre, resxsr(31, [("u235", m5, 21)], 1, "1.e-5", "2.e7", ".001"))

    m4 = mat("t404")
    pre = reconr(20, 21, [m4])
    case("err_endf4", ["t404"], pre, resxsr(31, [("t404", m4, 21)], 1, 1, 1000, ".001"))

    pre = reconr(20, 21, [m8])
    run = (resxsr(31, [("u238", m8, 21)], 1, 1, 1000, ".001") +
           resxsr(32, [("u238b", m8, 21)], 1, 5, 500, ".005"))
    case("err_again", [u8], pre, run)
    run = (resxsr(31, [("u238", m8, 21)], 1, 1, 1000, ".001") +
           resxsr(32, [("u238b", m8, 22)], 1, 5, 500, ".005"))
    case("again", [u8], pre + "moder\n 21 22/\n", run)

    pre = reconr(20, 21, [m8])
    case("u8", [u8], pre, resxsr(31, [("u238", m8, 21)], 1, 1, 10000, ".001"))
    print("%d cases in %s" % (n, corpus))


if __name__ == "__main__":
    main()
