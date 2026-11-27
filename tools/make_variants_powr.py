#!/usr/bin/env python3
"""Make the POWR cases (no NJOY2016 test uses POWR).

pre.input makes the GENDF tapes with the Fortran: MODER, RECONR, BROADR
(sometimes UNRESR or THERMR) and GROUPR runs, one per group structure:

  tape24  68 groups (GAM-I, ign=6)   for lib=1 (fast)
  tape25  35 groups (ign=8)          for lib=2 (thermal)
  tape26  69 groups (ign=9)          for lib=3 (cpm)
  tape27  185 groups (ign=1, user)   collapsed to the above by icgrp

run.input is POWR alone (one or more runs).

  f_*      lib=1, EPRI-CELL fast library (GAMTAP input)
  t_*      lib=2, EPRI-CELL thermal library (LIBRAR input)
  c_*      lib=3, EPRI-CPM library (CLIB input)
  err_*    error exits

lib=3 takes the card-5 values of each material from the scratch array
scr, not from the variables it reads them into (an NJOY bug): in the
first POWR run of a deck they are all zero (nina=0, ires=0, sigma zero
0, ...). A material with nina=0 then stops in sfile2 (no reference
sigma zero), so the cases that run to the end use nina>0 on card 5 for
every material, at least one ires=1 (or the deallocation of glamt
fails), if4=0, and mode 1 or 2 (mode 0 uses position 0 in idnto and
iflgs).

With --cards, only the k_* cases are made: lib=3 decks for NJOY2016 with
tools/njoy-powr-cards.patch, to be compared with the Rust port built with
--features powr-cards (see cards_cases).

Usage: make_variants_powr.py [--cards] NJOY2016_DIR CORPUS_DIR
"""
import math
import os
import sys


def pendf(mats, temps, tol=".005", unresr=None):
    """moder + reconr + broadr for mats on tape20; the PENDF is tape22."""
    s = "moder\n 20 -21/\n"
    s += "reconr\n -21 -22/\n 'pendf'/\n"
    for m in mats:
        s += " %d 0/\n %s/\n" % (m, tol)
    s += " 0/\n"
    nt = len(temps.split())
    s += "broadr\n -21 -22 -23/\n %d %d 0 0 0./\n %s/\n %s/\n" % (mats[0], nt, tol, temps)
    for m in mats[1:]:
        s += " %d/\n" % m
    s += " 0/\nmoder\n -23 -22/\n"
    if unresr:
        s += "unresr\n -21 -22 -23/\n %s/\n" % unresr
        s += "moder\n -23 -22/\n"
    return s


def e185():
    """186 group bounds, 1.e-5 eV to 2.e7 eV, equal lethargy."""
    lo, hi = math.log(1.e-5), math.log(2.e7)
    v = [math.exp(lo + (hi - lo) * i / 185) for i in range(186)]
    lines = []
    for i in range(0, 186, 6):
        lines.append(" " + " ".join("%.6e" % x for x in v[i:i + 6]))
    return " 185/\n" + "\n".join(lines) + "/\n"


def gendf(pend, nout, mat, ign, temps, sigz, rx, more=""):
    """One groupr run: mat at temps and sigma zeros, reaction list rx for
    each temperature (more: further materials)."""
    nt = len(temps.split())
    nz = len(sigz.split())
    s = "groupr\n -21 %s 0 %d/\n %d %d 0 3 1 %d %d 0/\n 'gendf'/\n %s/\n %s/\n" % (
        pend, nout, mat, ign, nt, nz, temps, sigz)
    if ign == 1:
        s += e185()
    s += (rx + " 0/\n") * nt
    s += more
    return s + " 0/\n"


def main():
    args = sys.argv[1:]
    cards = "--cards" in args
    if cards:
        args.remove("--cards")
    njoy, corpus = args[0], args[1]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    n = 0

    def case(name, pre_in, run_in, tapes):
        nonlocal n
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        for u, f in tapes:
            os.symlink(os.path.join(res, f), os.path.join(d, "tapes", "tape%d" % u))
        with open(os.path.join(d, "pre.input"), "w") as fo:
            fo.write(pre_in + "stop\n")
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(run_in + "stop\n")
        n += 1

    if cards:
        cards_cases(case)
        print("%d cases" % n)
        return
    c511 = [(20, "t511")]
    RX = " 3/\n 6/\n"
    RXU = " 3/\n 3 452/\n 3 455/\n 5 455/\n 6/\n 6 18/\n"

    def gset(mat, temps, sigz, rx, more="", ngout=(24, 25, 26, 27)):
        s = ""
        for u, ign in zip(ngout, (6, 8, 9, 1)):
            s += gendf("-22", -u, mat, ign, temps, sigz, rx, more)
        return s

    # carbon: one temperature and sigma zero; two of each
    PC1 = pendf([1306], "300.") + gset(1306, "300.", "1.e10", RX)
    PC2 = pendf([1306], "300. 600.") + gset(1306, "300. 600.", "1.e10 10.", RX)
    # U-235 (ENDF/B-V): fission, nubar, delayed neutrons; three sigma zeros
    PU2 = pendf([1395], "300. 600.") + gset(1395, "300. 600.", "1.e10 100. 1.", RXU)
    # carbon and U-235 on one GENDF tape
    PCU = pendf([1306, 1395], "300. 600.")
    for u, ign in zip((24, 25, 26, 27), (6, 8, 9, 1)):
        extra = " 1395/\n" + (RXU + " 0/\n") * 2
        PCU += gendf("-22", -u, 1306, ign, "300. 600.", "1.e10 10.", RX, extra)

    def fast(units, mats, card2="1 1 0"):
        """powr lib=1: mats is a list of (card3, word, fsn)."""
        s = "powr\n %s/\n %s/\n" % (units, card2)
        for c3, w, f in mats:
            s += " %s/\n '%s'/\n '%s'/\n" % (c3, w, f)
        return s + " 0/\n"

    f24 = "-24 30"
    case("f_c1", PC1, fast(f24, [("1306 300.", "carbon", "none")]), c511)
    case("f_c2", PC2, fast(f24, [("1306 300.", "carbon", "none")]), c511)
    case("f_c2_ip0", PC2, fast(f24, [("1306 300.", "carbon", "none")], "1 0"), c511)
    case("f_c2_rt", PC2, fast(f24, [("1306 600.", "carbon hot", "none")]), c511)
    case("f_u5", PU2, fast(f24, [("1395 300.", "u235 endf/b-v", "u235 fission spectrum")]), c511)
    case("f_u5_iz", PU2, fast(f24, [("1395 300. 1 2 2", "u235", "chi u235")]), c511)
    case("f_u5_nf", PU2, fast(f24, [("1395 300. 0", "u235 no ff", "chi")]), c511)
    case("f_u5_ip0", PU2, fast(f24, [("1395 300.", "u235", "chi")], "1 0"), c511)
    case("f_two", PCU, fast(f24, [("1306 300.", "carbon", "x"), ("1395 300.", "u235", "u235 chi")]), c511)
    case("f_twice", PCU, fast(f24, [("1395 300.", "u235", "chi")]) +
         fast("-24 31", [("1306 300.", "carbon", "x"), ("1395 300.", "u235", "chi 2")]), c511)
    case("f_185", PCU, fast("-27 30", [("1306 300.", "carbon", "x"), ("1395 300.", "u235", "chi")]), c511)
    case("f_nneg", PC1, fast("-24 -30", [("1306 300.", "carbon", "x")]), c511)
    case("f_asc", PC1 + "moder\n -24 34/\n", fast("34 30", [("1306 300.", "carbon", "x")]), c511)
    # Be-9 (n,2n), Pu-238 (ENDF/B-IV) with unresolved self-shielding
    PBE = pendf([425], "300.") + gset(425, "300.", "1.e10", RX)
    case("f_be", PBE, fast(f24, [("425 300.", "be9", "x")]), [(20, "n-004_Be_009-ENDF8.0.endf")])
    PPU = pendf([1050], "300.", unresr="1050 1 3 1/\n 300./\n 1.e10 1.e3 10.\n 0") + \
        gset(1050, "300.", "1.e10 1.e3 10.", RX + " 6 18/\n")
    case("f_pu", PPU, fast(f24, [("1050 300.", "pu238 b-iv", "pu chi")]), [(20, "t404")])

    # thermal library
    def therm(units, mats, card2="2 1"):
        """powr lib=2: mats is a list of (card3, card4, card5)."""
        s = "powr\n %s/\n %s/\n" % (units, card2)
        for c3, c4, c5 in mats:
            s += " %s/\n %s/\n %s/\n" % (c3, c4, c5)
        return s + " 0/\n"

    t25 = "-25 30"
    z8 = ".1 .2 .3 .4 .5 .6 .7 0."
    case("t_c", PC2, therm(t25, [("1306 300 'carbon'", "1 2 0", z8), ("1306 600 'carb600'", "0 2 0", "")]), c511)
    case("t_c0", PC2, therm(t25, [("1306 300 'carbon'", "0 0 0", z8)], "2 0"), c511)
    case("t_u5", PU2, therm(t25, [("1395 300 'u235'", "1 2 0", "1. .5 .1 2.4 1.9 .3 .8 9.")]), c511)
    case("t_two", PCU, therm(t25, [("1306 300 'carbon'", "1 2 0", z8), ("1395 600 'u-235 hot'", "0 2 2", z8)]),
         c511)
    case("t_185", PCU, therm("-27 30", [("1306 300 'carbon'", "1 2 0", z8), ("1395 300 'u235'", "1 2 0", z8)]),
         c511)
    case("t_twice", PC2, therm(t25, [("1306 300 'carbon'", "1 2 0", z8)]) +
         therm("-25 31", [("1306 600 'carbon'", "1 2 0", "")]), c511)
    # H-1 free gas from THERMR (mti=221)
    th = ("moder\n 20 -21/\nreconr\n -21 -30/\n 'pendf'/\n 125 0/\n .005/\n 0/\n"
          "broadr\n -21 -30 -31/\n 125 1 0 0 0./\n .005/\n 296./\n 0/\n"
          "thermr\n 0 -31 -22/\n 0 125 16 1 1 0 0 1 221 0/\n 296./\n .005 1./\n")
    PH = th + gendf("-22", -25, 125, 8, "296.", "1.e10", " 3/\n 6/\n 6 221/\n")
    case("t_h", PH, therm(t25, [("125 296 'h free gas'", "1 221 0", ".5 1. .67 0. 0. 0. 0. 20.")]),
         [(20, "n-001_H_001-ENDF8.0-Beta6.endf")])

    # cpm library
    def cpm(units, card2, card3, mats, card5s, rest):
        s = "powr\n %s/\n %s/\n %s/\n %s/\n" % (units, card2, card3, mats)
        for c in card5s:
            s += " %s/\n" % c
        return s + rest

    lam13 = " .1 .2 .3 .4 .5 .6 .7 .8 .9 1. 1.1 1.2 1.3/\n"
    lam15 = " .1 .2 .3 .4 .5 .6 .7 .8 .9 1. 1.1 1.2 1.3 1.4 1.5/\n"
    burn = " 2 1/\n 92235/\n 54135 2.1e-5 .065/\n 53135 2.9e-5 .063/\n"
    c26 = "-26 30"
    case("c_burn", PC2, cpm(c26, "3 1 1", "1 20260101 1 0 1 1 0", "1306", ["1 0 0 1.e10 1 4.7"], burn + lam13),
         c511)
    case("c_burn15", PC2, cpm(c26, "3 0 0", "7 1999 1 0 1 1 0", "1306", ["2 0 0 1.e10 1 4.7"], burn + lam15),
         c511)
    case("c_mode2", PU2, cpm(c26, "3 1 1", "2 20260102 1 0 2 0 0", "1395", ["1 0 0 1.e10 1 11.6"], lam13), c511)
    case("c_mode2_185", PCU, cpm("-27 30", "3 1 0", "2 20260103 2 0 2 0 0", "1306 1395",
                                 ["1 0 0 1.e10 1 4.7", "2 0 0 1.e10 1 11.6"], lam15 + lam15), c511)
    case("c_two", PCU, cpm(c26, "3 0 1", "4 1 2 0 1 1 0", "1395 -1306",
                           ["1 2 3 1.e10 1 11.6", "2 1 1 1.e10 0 4.7"],
                           " 3 2/\n 92235 94239/\n 54135 2.1e-5 .065 .07/\n 53135 2.9e-5 0. 0./\n"
                           " 55135 0. .01 .02/\n" + lam13), c511)
    case("c_yld0", PC2, cpm(c26, "3 1 1", "1 2 1 0 2 1 0", "1306", ["1 0 0 1.e10 1 4.7"],
                            " 1 1/\n 92235/\n 54135 2.1e-5 0./\n" + lam13), c511)
    # a lib=1 run first leaves its last GENDF record in scr (after a
    # lib=2 run, the values there send pinit past its material search and
    # the Fortran prints its uninitialised group-bound array: not used)
    case("c_after_fast", PCU, fast(f24, [("1306 300.", "carbon", "x")]) +
         cpm(c26, "3 1 1", "1 2 1 0 1 0 0", "1306", ["1 0 0 1.e10 1 4.7"], lam13), c511)


    # error exits
    case("err_lib", PC1, "powr\n -24 30/\n 4/\n", c511)
    case("err_f_mat", PC1, fast(f24, [("1395 300.", "u", "x")]), c511)
    case("err_f_rt", PC2, fast(f24, [("1306 900.", "c", "x")]), c511)
    case("err_f_grp", PC1, fast("-25 30", [("1306 300.", "c", "x")]), c511)
    PC6 = pendf([1306], "300. 400. 500. 600. 700. 800. 900.") + \
        gendf("-22", -24, 1306, 6, "300. 400. 500. 600. 700. 800. 900.", "1.e10 10.", RX)
    case("err_f_nt", PC6, fast(f24, [("1306 300.", "c", "x")]), c511)
    PC8 = pendf([1306], "300.") + gendf("-22", -24, 1306, 6, "300.", "1.e10 1.e5 1.e4 1.e3 100. 10. 1. .1", RX)
    case("err_f_sz", PC8, fast(f24, [("1306 300.", "c", "x")]), c511)
    case("err_t_temp", PC1, therm(t25, [("1306 400 'carbon'", "1 2 0", z8)]), c511)
    case("err_t_ngn", PC1, therm(f24, [("1306 300 'carbon'", "1 2 0", z8)]), c511)
    PCNF = pendf([1306], "300.") + gendf("-22", -25, 1306, 8, "300.", "1.e10", " 3 2/\n 6 2/\n")
    case("err_t_flux", PCNF, therm(t25, [("1306 300 'carbon'", "1 2 0", z8)]), c511)
    case("err_c_sz", PC2, cpm(c26, "3 1 1", "1 2 1 0 1 0 0", "1306", ["0 0 0 1.e10 1 4.7"], lam13), c511)
    case("err_c_glam", PC2, cpm(c26, "3 1 1", "1 2 1 0 1 0 0", "1306", ["1 0 0 1.e10 0 4.7"], ""), c511)
    case("err_c_mat", PC2, cpm(c26, "3 1 1", "1 2 1 0 1 0 0", "1395", ["1 0 0 1.e10 1 4.7"], lam13), c511)
    case("err_c_if5", PC2, cpm(c26, "3 1 1", "1 2 1 0 1 2 0", "1306", ["1 0 0 1.e10 1 4.7"],
                               burn + " 1. 0 1/\n 2. 0 1/\n" + lam13), c511)
    case("err_c_twice", PC2, cpm(c26, "3 1 1", "1 2 1 0 1 0 0", "1306", ["1 0 0 1.e10 1 4.7"], lam13) * 2,
         c511)
    case("err_c_grp", PC2, cpm(f24, "3 1 1", "1 2 1 0 1 0 0", "1306", ["1 0 0 1.e10 1 4.7"], lam13), c511)
    print("%d cases" % n)


def cards_cases(case):
    """lib=3 cases for NJOY2016 with tools/njoy-powr-cards.patch (and the
    Rust port built with --features powr-cards): the card-5 values and
    lambdas are kept, so subfiles 2 to 5 and the CLIB output run.

    The GENDF tape (69 groups, 300 and 600 K) has U-235 (ENDF/B-V) with
    elastic, fission and a free-gas thermal matrix (THERMR, mt 221), and
    carbon. Constraints of the Fortran, kept by every deck:
      * every scattering-matrix row has at most 5 groups, or the patched
        format is needed (the patch has it; real data: only U-235 here)
      * ip1opt=1 for nina=0 materials (sfile6 skips every MF6 section and
        stops with "no required p1 matrices")
      * at least one material with ires=1 (indf is set in sfile3 only)
      * ntemp > 0 on card 5 (sf4out copies ntemp temperatures, cpmout
        reads the number it found)
      * ires=1 only with nina 1 or 2 (sf3out runs for nina>0 only, and
        sfile3 would search the tape forever for a nina=3 ident)
    """
    c511 = [(20, "t511")]
    pre = ("moder\n 20 -21/\nreconr\n -21 -22/\n 'pendf'/\n 1306 0/\n .005/\n 1395 0/\n .005/\n 0/\n"
           "broadr\n -21 -22 -23/\n 1306 2 0 0 0./\n .005/\n 300. 600./\n 1395/\n 0/\n"
           "thermr\n 0 -23 -24/\n 0 1395 16 2 1 0 0 1 221 0/\n 300. 600./\n .005 4./\n")
    ru = " 3/\n 3 452/\n 3 455/\n 5 455/\n 6 2/\n 6 18/\n 6 221/\n 0/\n"
    rc = " 3/\n 6/\n 0/\n"
    # THERMR writes U-235 only; carbon is grouped from the BROADR tape and
    # added to the GENDF tape of U-235 (input gout tape)
    pre += ("groupr\n -21 -24 0 -25/\n 1395 9 0 3 1 2 2 0/\n 'g'/\n 300. 600./\n 1.e10 10./\n"
            + ru * 2 + " 0/\n")
    pre += ("groupr\n -21 -23 -25 -26/\n 1306 9 0 3 1 2 2 0/\n 'g'/\n 300. 600./\n 1.e10 10./\n"
            + rc * 2 + " 0/\n")
    lam = " .1 .2 .3 .4 .5 .6 .7 .8 .9 1. 1.1 1.2 1.3/\n"
    g = 69

    def num(v):
        return " ".join("%g" % x for x in v)

    # nina=3 data for an ident not on the tape: card 11, 11a (sfile2),
    # card 12, 12a, 12b/12c (sfile4)
    def nina3(scale):
        s11 = " %s/\n %s/\n" % (num([2.4 + .01 * i for i in range(13)]), num([10. + i for i in range(13)]))
        s12 = " 239. 300. %s/\n" % num([scale * (1. + .1 * i) for i in range(g)])
        s12 += " %s/\n %s/\n %s/\n" % (num([.5 * scale] * g), num([2.9] * g), num([3. + .01 * i for i in range(g)]))
        for ia in range(1, g + 1):
            if ia % 7 == 0:
                s12 += " 0 0 0/\n"
            else:
                hi = min(ia + 2, g)
                s12 += " %d %d %d/\n %s/\n" % (ia, ia, hi, num([.3 * (k + 1) for k in range(hi - ia + 1)]))
        return s11, s12

    def nina2():
        return " 16. 300. %s/\n" % num([1.e-4 * (i + 1) for i in range(g)])

    def deck(card2, card3, mats, c5, extra_before_lam="", nlam=1, s2="", s4=""):
        s = "powr\n -26 30/\n %s/\n %s/\n %s/\n" % (card2, card3, mats)
        for c in c5:
            s += " %s/\n" % c
        return s + extra_before_lam + lam * nlam + s2 + s4

    U = "0 2 0 1.e10 0 11.6 221 0 1 0 50 0"
    C = "1 2 0 1.e10 1 4.7 0 0 0 0 5 0"
    P3 = "3 1 0 1.e10 0 10. 0 0 0 0 55 0"
    O2 = "2 1 0 1.e10 0 3.8 0 0 0 0 6 0"
    s11, s12 = nina3(1.)
    case("k_uc", pre, deck("3 1 1", "1 20260101 2 0 1 0 1", "1395 1306", [U, C]), c511)
    case("k_uc_ip0", pre, deck("3 0 1", "1 20260101 2 0 1 0 1", "1395 1306", [U, C]), c511)
    case("k_mode0", pre, deck("3 1 1", "1 20260101 2 0 0 0 1", "1395 1306", [U, C]), c511)
    case("k_mode2", pre, deck("3 1 1", "1 20260101 2 0 2 0 1", "1395 1306", [U, C]), c511)
    case("k_all", pre, deck("3 1 1", "1 20260101 4 0 0 0 1", "1395 1306 94239 8016", [U, C, P3, O2],
                            nlam=2, s2=s11, s4=s12 + nina2()), c511)
    # (the first material must be on the tape: pinit reads the group
    # bounds only there, else the Fortran uses its unset array)
    case("k_all1", pre, deck("3 1 1", "2 3 4 0 1 0 1", "1306 8016 1395 94239", [C, O2, U, P3],
                             nlam=2, s2=s11, s4=nina2() + s12), c511)
    case("k_burn", pre, deck("3 1 1", "1 20260101 2 0 1 1 1", "1395 1306", [U, C],
                             " 2 1/\n 92235/\n 54135 2.1e-5 .065/\n 53135 2.9e-5 .063/\n"), c511)
    case("k_u1t", pre, deck("3 1 1", "1 20260101 2 0 1 0 1", "1395 1306",
                            ["0 1 1 1.e10 0 11.6 221 0 1 0 50 0", "1 1 1 1.e10 1 4.7 0 0 0 0 5 0"]), c511)
    case("k_szref", pre, deck("3 1 1", "1 20260101 2 0 1 0 1", "1395 1306",
                              ["0 2 0 10. 0 11.6 221 0 1 0 50 0", "1 2 0 10. 1 4.7 0 0 0 0 5 0"]), c511)
    # U-235 as the resonance material (fission in sfile3), carbon with all
    # its matrices (mti=51 is carbon's first inelastic level)
    UR = "1 2 0 1.e10 1 11.6 0 0 0 0 50 1"
    CF = "0 2 0 1.e10 0 4.7 51 0 1 0 5 0"
    case("k_ures", pre, deck("3 1 1", "1 20260101 2 0 1 0 1", "1306 1395", [CF, UR]), c511)
    case("k_ures0", pre, deck("3 1 1", "1 20260101 2 0 0 0 1", "1306 1395", [CF, UR]), c511)
    case("k_ures2", pre, deck("3 0 1", "1 20260101 2 0 2 1 1", "1306 1395", [CF, UR],
                              " 1 1/\n 92235/\n 54135 2.1e-5 .065/\n"), c511)
    case("k_norf", pre, deck("3 1 1", "1 20260101 2 0 1 0 1", "1306 1395",
                             ["1 2 0 1.e10 1 4.7 0 0 0 1 5 0", "2 2 0 1.e10 1 11.6 0 0 0 1 50 0"],
                             nlam=1, s4=" 235. 300. %s/\n" % " ".join(["1."] * 69)), c511)


if __name__ == "__main__":
    main()
