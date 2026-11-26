#!/usr/bin/env python3
"""Make MIXR cases. No NJOY2016 test deck calls mixr, so all cases are here.

Each case links one or more ENDF files from the NJOY2016 test resources
(tape20, tape30, tape40, ...), runs the Fortran RECONR (and for some cases
BROADR or MODER) in pre.input to make PENDF tapes, and then MIXR in
run.input. Some cases give MIXR the ENDF tapes directly.

  pair_*     two materials, a few MTs, positive weights
  miss_*     an MT that is not on every tape (message; the Fortran then
             skips the first MF3 section of that tape for the next MT)
  none       an MT that is on no tape (no TAB1 record)
  big_*      many points (pages in gety and in the output TAB1)
  endf_*     raw ENDF tapes (interpolation ranges, other versions, nsub)
  temp_*     BROADR tapes; MIXR picks the second or third temperature
  bin        binary input and output tapes (MODER)
  max        ten materials and twenty MTs
  neg        a negative weight
  err_mat    a material that is not on the tape (error exit)
  err_temp   a temperature that is not on the tape (error exit)

Usage: make_variants_mixr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

UNITS = [20, 30, 40, 50, 60, 70, 80, 25, 35, 45]


def first_mat(path):
    with open(path, errors="replace") as f:
        f.readline()
        return int(f.readline()[66:70])


def reconr(nin, nout, mat, err=".01"):
    return "reconr\n %d %d/\n 'mixr variant'/\n %d 0/\n %s/\n 0/\n" % (nin, nout, mat, err)


def mixr(nout, nins, mts, mw, temp, matd, za, awr, text="mixed cross sections"):
    s = "mixr\n %d %s/\n" % (nout, " ".join(str(x) for x in nins))
    s += " %s/\n" % " ".join(str(x) for x in mts)
    s += " %s/\n" % " ".join("%d %s" % (m, w) for m, w in mw)
    s += " %s/\n %d %s %s/\n '%s'/\n" % (temp, matd, za, awr, text)
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

    def pendf(files, err=".01"):
        """reconr each file u -> u+1; return (pre, pendf units, mats)."""
        pre, outs, mats = "", [], []
        for u, f in zip(UNITS, files):
            m = mat(f)
            pre += reconr(u, u + 1, m, err)
            outs.append(u + 1)
            mats.append(m)
        return pre, outs, mats

    h1 = "n-001_H_001-ENDF8.0-Beta6.endf"
    b10 = "n-005_B_010-ENDF8.0.endf"
    al = "n-013_Al_027-ENDF8.0.endf"
    fe56 = "n-026_Fe_056-JEFF3.3.endf"
    fe58 = "n-026_Fe_058-ENDF8.0-Beta4.endf"
    u235 = "n-092_U_235-ENDF8.0.endf"
    pu239 = "n-094_Pu_239-ENDF8.0-Beta6.endf"
    pu240 = "n-094_Pu_240-ENDF8.0.endf"
    zr90 = "n-040_Zr_090-ENDF8.0.endf"
    be9 = "n-004_Be_009-ENDF8.0.endf"

    # two materials
    for name, fs, mts, w in [
            ("pair_hb", [h1, b10], [1, 2, 102], ["2.0", "0.5"]),
            ("pair_fe", [fe56, fe58], [1, 2, 4, 102], [".9175", ".0028"]),
            ("pair_al_be", [al, be9], [1, 2, 16, 102, 103, 107], ["1", "1"]),
            ("pair_pu", [pu239, pu240], [1, 2, 18, 102], [".94", ".06"])]:
        pre, outs, mats = pendf(fs)
        case(name, fs, pre, mixr(90, outs, mts, zip(mats, w), 0, 9999, 1000., 10.,
                                 "mix of " + name))

    # an MT not on every tape
    fs = [u235, fe56]
    pre, outs, mats = pendf(fs)
    case("miss_18", fs, pre, mixr(90, outs, [1, 18, 2, 102], zip(mats, ["1", "1"]),
                                   0, 9228, 92235., 233.))
    case("miss_first", fs, pre, mixr(90, outs, [18, 1, 102], zip(mats, ["1", "1"]),
                                      0, 9228, 92235., 233.))
    fs = [h1, b10]
    pre, outs, mats = pendf(fs)
    case("none", fs, pre, mixr(90, outs, [1, 999, 2], zip(mats, ["1", "1"]), 0, 1, 1., 1.))

    # many points
    fs = [u235, pu239]
    pre, outs, mats = pendf(fs, ".001")
    case("big_upu", fs, pre, mixr(90, outs, [1, 2, 18, 102], zip(mats, [".5", ".5"]),
                                  0, 9200, 92000., 235.))
    fs = [u235]
    pre, outs, mats = pendf(fs, ".001")
    case("big_u", fs, pre, mixr(90, outs, [1, 102], [(mats[0], "1.0")], 0, 9228, 92235., 233.))

    # raw ENDF tapes
    for name, fs, mts in [
            ("endf_fe", [fe56, fe58], [1, 2, 102]),
            ("endf_zr_al", [zr90, al], [1, 2, 4, 102]),
            ("endf_hb", [h1, b10], [1, 2, 102]),
            ("endf_t404", ["t404"], [1, 2, 102]),
            ("endf_t511", ["t511"], [1, 2, 102]),
            ("endf_p", ["p-001_H_001-ENDF8.0.endf"], [2]),
            ("endf_d", ["d-001_H_002-ENDF8.0.endf"], [2]),
            ("endf_he3", ["d-002_He_003.endf"], [2]),
            ("endf_a", ["a-002_He_004-ENDF8.0.endf"], [2]),
            ("endf_p2", ["p-001_H_002-ENDF8.0.endf"], [2]),
            ("endf_p3", ["p-001_H_003-ENDF8.0.endf"], [2])]:
        mats = [mat(f) for f in fs]
        case(name, fs, "", mixr(90, UNITS[:len(fs)], mts, [(m, "1.") for m in mats],
                                0, 100, 1., 1.))

    # temperatures
    for name, temp in [("temp_900", 900), ("temp_1800", 1800), ("temp_300", 300.4)]:
        fs = [fe56, al]
        pre = ""
        outs, mats = [], []
        for u, f in zip(UNITS, fs):
            m = mat(f)
            pre += reconr(u, u + 2, m)
            pre += "broadr\n %d %d %d/\n %d 3 0 0 0./\n .01/\n 300 900 1800/\n 0/\n" % (
                u, u + 2, u + 1, m)
            outs.append(u + 1)
            mats.append(m)
        case(name, fs, pre, mixr(90, outs, [1, 2, 102], zip(mats, [".5", ".5"]), temp,
                                 2600, 26000., 55.))

    # binary tapes
    fs = [fe56, al]
    pre, outs, mats = pendf(fs)
    pre += "moder\n 21 -23/\nmoder\n 31 -33/\n"
    case("bin", fs, pre, mixr(-90, [-23, -33], [1, 2, 102], zip(mats, ["1", "2"]), 0,
                              2600, 26000., 55.))

    # ten materials, twenty MTs
    fs = [h1, b10, al, fe56, fe58, zr90, be9, u235, pu239, pu240]
    pre, outs, mats = pendf(fs, ".05")
    mts = [1, 2, 4, 16, 17, 18, 22, 28, 51, 52, 53, 91, 102, 103, 104, 105, 106, 107, 111, 600]
    case("max", fs, pre, mixr(90, outs, mts, zip(mats, [".1"] * 10), 0, 1, 1., 1.,
                              "ten materials, twenty mts, a long text " + "x" * 40))

    # a negative weight
    fs = [fe56, fe58]
    pre, outs, mats = pendf(fs)
    case("neg", fs, pre, mixr(90, outs, [1, 2], zip(mats, ["1.", "-1."]), 0, 2600, 26000., 55.))

    # errors
    fs = [h1]
    pre, outs, mats = pendf(fs)
    case("err_mat", fs, pre, mixr(90, outs, [1], [(9999, "1.")], 0, 1, 1., 1.))
    case("err_temp", fs, pre, mixr(90, outs, [1], [(mats[0], "1.")], 293.6, 1, 1., 1.))
    print("%d cases in %s" % (n, corpus))


if __name__ == "__main__":
    main()
