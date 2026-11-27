#!/usr/bin/env python3
"""Make extra GAMINR cases.

The test suite calls GAMINR twice (tests 03 and 30). Here pre.input makes
the PENDF tape (RECONR, run by the Fortran) and run.input is GAMINR alone.

  g_*      H, Fe and U photon data of test 03 (gam23 + gam27, ENDF-4/5 style):
           each group structure and weight option, Legendre orders up to
           7, print options, explicit reaction lists, tape options
  t_*      GAMINR twice in one run (module values kept from run to run)
  in_*     an input GAMINR tape (ngam1): material lower, equal, higher;
           in_mat1 stops at end of file in both codes (mat 1 is taken
           as the tape header)
  pho_*    U and H photoatomic data (ENDF-6, MT522/525)
  err_*    error exits

Usage: make_variants_gaminr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

PRE03 = """reconr
 30 31/
 'pendf tape for photon interaction cross sections from dlc7e'/
 1 1 0/
 .001/
 '1-hydrogen'/
 92 1 0/
 .001/
 '92-uranium'/
 0/
"""

PRE3M = """reconr
 30 31/
 'pendf tape for h, fe and u'/
 1 1 0/
 .001/
 '1-hydrogen'/
 26 1 0/
 .001/
 '26-iron'/
 92 1 0/
 .001/
 '92-uranium'/
 0/
"""

PREU = """reconr
 20 21/
 'u photoatomic'/
 9200 0 0/
 .001/
 0/
"""

PREH = """reconr
 20 21/
 'h photoatomic'/
 100 0 0/
 .001/
 0/
"""

WT = "0. 0. 0 0 1 4 4 2 1.e3 1. 1.e5 .5 1.e7 2. 1.e9 1./\n"
GR = " 8/\n 1.e4 5.e4 1.e5 2.e5 3.e5 5.e5 7.e5 9.e5 1.e6/\n"


def many_groups(ng):
    """A read-in structure of ng log-spaced groups from 10 keV to 20 MeV."""
    e = ["%.6e" % (1.e4 * (2.e3 ** (i / ng))) for i in range(ng + 1)]
    lines = [" " + " ".join(e[i:i + 6]) for i in range(0, len(e), 6)]
    return " %d/\n" % ng + "\n".join(lines) + "/\n"


def gam(units, card2, mats, title="photon test", groups="", wt=""):
    """units: card 1; card2: 'igg iwt lord iprint' (mat put first);
    mats: list of (mat, [reaction cards]) where an empty list means -1."""
    s = "gaminr\n %s/\n %d %s/\n '%s'/\n" % (units, mats[0][0], card2, title)
    s += groups + wt
    for k, (m, cards) in enumerate(mats):
        if k > 0:
            s += " %d/\n" % m
        if cards:
            for cd in cards:
                s += " %s/\n" % cd
            s += " 0/\n"
        else:
            s += " -1/\n"
    return s + " 0/\n"


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    n = 0

    def case(name, pre, run_in, tapes):
        nonlocal n
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        for u, f in tapes:
            os.symlink(os.path.join(res, f), os.path.join(d, "tapes", "tape%d" % u))
        with open(os.path.join(d, "pre.input"), "w") as fo:
            fo.write(pre + "stop\n")
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(run_in + "stop\n")
        n += 1

    t03 = [(30, "gam23"), (32, "gam27")]
    tu = [(20, "photoat-092_U_000-ENDF8.0.endf")]
    th = [(20, "photoat-001_H_000-ENDF8.0.endf")]
    both = [(1, []), (92, [])]

    # group structures, weight options and Legendre orders
    case("g_igg2", PRE03, gam("32 31 0 33", "2 3 3 1", both), t03)
    case("g_igg4", PRE03, gam("32 31 0 33", "4 2 0 1", both), t03)
    case("g_igg5", PRE03, gam("32 31 0 33", "5 1 2 1", both, wt=WT), t03)
    case("g_igg6", PRE03, gam("32 31 0 33", "6 3 5 1", both), t03)
    case("g_igg7", PRE03, gam("32 31 0 33", "7 3 7 1", both), t03)
    case("g_igg8", PRE03, gam("32 31 0 33", "8 2 1 1", [(92, [])]), t03)
    case("g_igg9", PRE03, gam("32 31 0 33", "9 3 4 1", [(1, [])]), t03)
    case("g_igg10", PRE03, gam("32 31 0 33", "10 3 6 1", both), t03)
    # read-in groups that end below the pair threshold, explicit reactions
    rl = ["23 501 'total'", "26 502 'coherent'", "26 504 'incoh'", "26 516 'pair'",
          "23 516 'pair'", "23 602 'abs'", "23 621 'heat'"]
    case("g_read", PRE03, gam("32 31 0 33", "1 3 3 1", [(1, rl), (92, rl)], groups=GR), t03)
    # more than npage words in the head record and in the matrix records
    case("g_many", PRE03, gam("32 31 0 33", "1 3 1 1", [(1, [])], groups=many_groups(320)), t03)
    case("g_iprint0", PRE03, gam("32 31 0 33", "3 3 3 0", both), t03)
    case("g_nogam", PRE03, gam("32 31 0 0", "3 2 2 1", both), t03)
    case("g_bin", PRE03, gam("32 31 0 -33", "3 3 3 1", both), t03)
    case("g_list", PRE03, gam("32 31 0 33", "3 3 2 1",
                              [(1, ["23 504", "26 504 'incoherent'", "23 621"]), (92, [])]), t03)
    case("g_three", PRE3M, gam("32 31 0 33", "3 3 3 1", [(1, []), (26, []), (92, [])]), t03)
    case("g_three2", PRE3M, gam("32 31 0 33", "4 3 5 1", [(26, []), (92, [])]), t03)
    case("g_u", PRE03, gam("32 31 0 33", "3 3 4 1", [(92, [])]), t03)
    # two runs in one deck
    case("t_two", PRE03, gam("32 31 0 33", "3 3 3 1", [(1, [])])
         + gam("32 31 0 -34", "5 2 2 1", [(92, [])]), t03)
    case("t_bin1", PRE03, gam("32 31 0 -33", "3 3 1 1", [(1, [])])
         + gam("32 31 0 34", "3 3 1 1", [(92, ["23 501", "26 502"])]), t03)
    # input gamout tape: made in pre.input, then merged by the run
    case("in_lo", PRE3M + gam("32 31 0 33", "3 3 2 1", [(26, [])]),
         gam("32 31 33 34", "3 3 2 1", [(92, [])]), t03)
    case("in_eq", PRE3M + gam("32 31 0 33", "3 3 2 1", [(26, []), (92, [])]),
         gam("32 31 33 34", "3 3 2 1", [(26, ["23 501", "26 502"])]), t03)
    case("in_hi", PRE3M + gam("32 31 0 33", "3 3 2 1", [(92, [])]),
         gam("32 31 33 34", "3 3 2 1", [(26, [])]), t03)
    case("in_two", PRE3M + gam("32 31 0 33", "3 3 2 1", [(26, [])]),
         gam("32 31 33 34", "3 3 2 1", [(1, ["23 501"]), (26, ["23 501", "26 504"]), (92, [])]), t03)
    # mat 1 on the input tape is taken as a tape header: end of file (exit 2)
    case("in_mat1", PRE03 + gam("32 31 0 33", "3 3 2 1", [(1, [])]),
         gam("32 31 33 34", "3 3 2 1", [(92, [])]), t03)
    # ENDF-6 photoatomic data
    case("pho_u", PREU, gam("20 21 0 22", "3 3 3 1", [(9200, [])]), tu)
    case("pho_u8", PREU, gam("20 21 0 22", "10 3 8 1", [(9200, [])]), tu)
    case("pho_h", PREH, gam("20 21 0 22", "1 1 2 1", [(100, [])], groups=GR, wt=WT), th)
    case("pho_h2", PREH, gam("20 21 0 -22", "7 2 0 1", [(100, ["23 522", "23 525", "26 516"])]), th)
    # error exits
    case("err_igg", PRE03, gam("32 31 0 33", "11 3 3 1", both), t03)
    case("err_iwt", PRE03, gam("32 31 0 33", "3 4 3 1", both), t03)
    case("err_ngg", PRE03, gam("32 31 0 33", "1 3 3 1", both, groups=" 400/\n 1.e4 2.e4/\n"), t03)
    case("err_mt", PREH, gam("20 21 0 22", "3 3 3 1", [(100, ["23 515"])]), th)
    case("err_mf", PRE03, gam("32 31 0 33", "3 3 3 1", [(1, ["24 501"])]), t03)
    print("cases: %d" % n)


if __name__ == "__main__":
    main()
