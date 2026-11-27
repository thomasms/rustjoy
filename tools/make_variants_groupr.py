#!/usr/bin/env python3
"""Make extra GROUPR cases.

pre.input makes the PENDF tape (MODER, RECONR, BROADR and sometimes
UNRESR, THERMR or a first GROUPR run, by the Fortran); run.input is
GROUPR alone.

  ign_*    each built-in neutron group structure (C-nat, ENDF/B-V)
  igg_*    each gamma group structure (C-nat photon production)
  iwt_*    each weight function, the flux calculator (iwt<0) and a
           flux read back from ninwt (iwt=0)
  opt_*    Legendre orders, print option, two temperatures, sigma
           zeroes with unresolved data, auto lists, special mt numbers
  out_*    gout tapes: replace a material, add mts, a material below the
           first one on the tape, two materials, all materials
  u5_*     U-235: fission spectra, nubar, delayed neutrons, photons
  t76      the GROUPR call of test 76 (its deck quotes module names, so
           build_corpus.py does not find it)
  mf6_*    File 4, 5 and 6 matrices (laws 1, 2, 3, 6, 7, discrete
           capture gamma), particle and nuclide production
  th_*     thermal MF6 data from THERMR (H in ZrH, Al)
  cp_*     incident protons, deuterons and alphas
  g_*      incident photons (photonuclear data)
  err_*    error exits

Usage: make_variants_groupr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys


def pre(mats, files, temps="296.", tol=".005", unresr=None, extra=""):
    """moder + reconr + broadr for mats (list of mat numbers) on tape20."""
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
    return s + extra


def gr(units, mat, card2, title, temps, sigz, cards, more=""):
    """A groupr deck. card2 is 'ign igg iwt lord ntemp nsigz iprint ismooth'
    (matb is put first); cards are the group/weight cards and reactions."""
    s = "groupr\n %s/\n %d %s/\n '%s'/\n %s/\n %s/\n" % (units, mat, card2, title, temps, sigz)
    s += cards
    s += more
    return s


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
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

    c511 = [(20, "t511")]
    PC = pre([1306], "t511")
    U = "-21 -22 0 -24"

    # every built-in neutron group structure
    for ign in range(2, 37):
        case("ign_%d" % ign, PC, gr(U, 1306, "%d 0 3 0 1 1 1" % ign, "c ign %d" % ign, "296.", "1.e10",
                                    " 3 1 'total'/\n 3 2/\n 0/\n 0/\n"), c511)
    # every gamma group structure, photon production cross sections and matrices
    for igg in range(1, 11):
        g = " 6/\n 1.e4 1.e5 1.e6 3.e6 7.e6 2.e7/\n" if igg == 1 else ""
        case("igg_%d" % igg, PC, gr(U, 1306, "3 %d 3 2 1 1 1" % igg, "c igg %d" % igg, "296.", "1.e10",
                                    g + " 13/\n 16/\n 0/\n 0/\n"), c511)

    # weight functions (iwt 1 and 4 read cards)
    wt1 = " 0. 0. 0 0 1 4 4 2 1.e-5 1. 1.e3 .5 1.e6 2. 2.e7 1./\n"
    wt4 = " .1 .025 .8208e6 1.4e6/\n"
    for iwt in range(1, 13):
        w = wt1 if iwt == 1 else wt4 if iwt == 4 else ""
        case("iwt_%d" % iwt, PC, gr(U, 1306, "3 0 %d 1 1 1 1" % iwt, "c iwt %d" % iwt, "296.", "1.e10",
                                    w + " 3 1/\n 3 2/\n 6 2/\n 0/\n 0/\n"), c511)

    # self-shielding with unresolved data: Pu-238 (ENDF/B-IV), 3 sigma zeroes
    p404 = [(20, "t404")]
    PPU = pre([1050], "t404", temps="300.",
              unresr="1050 1 3 1/\n 300./\n 1.e10 1.e3 10./\n 0")
    case("opt_pu", PPU, gr(U, 1050, "3 0 3 3 1 3 1", "pu238 sigz", "300.", "1.e10 1.e3 10.",
                           " 3 1/\n 3 2/\n 3 18/\n 3 102/\n 6 2/\n 6 18/\n 0/\n 0/\n"), p404)
    case("opt_pu0", PPU, gr(U, 1050, "3 0 3 1 1 3 0", "pu238 sigz", "300.", "1.e10 1.e3 10.",
                            " 3 1/\n 3 2/\n 6 2/\n 0/\n 0/\n"), p404)
    # flux calculator (iwt<0), writing ninwt, then a run with iwt=0 reading it
    fc = " 50. 10. 3000 30 2 0 0 0 0 0/\n"
    case("iwt_m3", PPU, gr("-21 -22 0 -24", 1050, "3 0 -3 1 1 3 1", "pu flux", "300.", "1.e10 1.e3 10.",
                           fc + " 3 1/\n 3 2/\n 3 102/\n 0/\n 0/\n"), p404)
    case("iwt_m4", PPU, gr("-21 -22 0 -24", 1050, "3 0 -4 1 1 3 1", "pu flux", "300.", "1.e10 1.e3 10.",
                           " 50. 10. 3000 0 1 .1 2. .2 .05 .3/\n" + wt4 + " 3 1/\n 3 102/\n 0/\n 0/\n"), p404)
    case("iwt_0", PPU + gr("-21 -22 0 -24", 1050, "3 0 -3 1 1 3 1", "pu flux", "300.", "1.e10 1.e3 10.",
                           fc + " 3 1/\n 0/\n 0/\n"),
         gr("-21 -22 0 -25", 1050, "3 0 0 1 1 3 1", "pu ninwt", "300.", "1.e10 1.e3 10.",
            " 30/\n 3 1/\n 3 2/\n 3 102/\n 0/\n 0/\n"), p404)

    # options on C: Legendre orders, print, two temperatures, auto lists
    PC2 = pre([1306], "t511", temps="296. 600.")
    for lord in (0, 3, 5, 8):
        case("opt_lord%d" % lord, PC, gr(U, 1306, "3 3 3 %d 1 1 1" % lord, "c lord", "296.", "1.e10",
                                         " 6/\n 16/\n 0/\n 0/\n"), c511)
    case("opt_ip0", PC, gr(U, 1306, "3 3 3 1 1 1 0", "c iprint 0", "296.", "1.e10",
                           " 3/\n 6/\n 13/\n 16/\n 0/\n 0/\n"), c511)
    case("opt_t2", PC2, gr(U, 1306, "3 3 3 1 2 1 1", "c two temps", "296. 600.", "1.e10",
                           " 3 1/\n 3 2/\n 6 2/\n 16 102/\n 0/\n 3 1/\n 6 2/\n 0/\n 0/\n"), c511)
    case("opt_auto3", PC, gr(U, 1306, "5 0 3 0 1 1 1", "c auto", "296.", "1.e10", " 3/\n 0/\n 0/\n"), c511)
    case("opt_mts", PC, gr(U, 1306, "9 0 3 2 1 1 1", "c special", "296.", "1.e10",
                           " 3 251/\n 3 252/\n 3 253/\n 3 257/\n 3 258/\n 3 259/\n 3 -4/\n 3 51 'x'/\n"
                           " 3 -55/\n 0/\n 0/\n"), c511)
    case("opt_sz", PC, gr(U, 1306, "3 0 3 1 1 3 1", "c sigz", "296.", "1.e10 100. 1.",
                          " 3 1/\n 3 2/\n 6 2/\n 0/\n 0/\n"), c511)

    # gout tapes: replace a material on an input gout tape, add mts to a
    # material (matb<0 on card 10: on card 2 NJOY2016 does not make matb
    # positive, and findf fails), a material below the first one on the
    # tape, all materials of the ENDF tape (matb<0, no input gout tape:
    # NJOY2016 takes the tape label as the first material and fails)
    PCU = pre([1306, 1395], "t511")
    g2 = gr("-21 -22 0 -24", 1306, "3 0 3 1 1 1 1", "c", "296.", "1.e10",
            " 3 1/\n 3 4/\n 0/\n", " 1395/\n 3 1/\n 3 18/\n 0/\n 0/\n")
    case("out_rep", PCU + g2, gr("-21 -22 -24 -25", 1395, "3 0 3 1 1 1 1", "u rep", "296.", "1.e10",
                                 " 3 1/\n 3 2/\n 0/\n 0/\n"), c511)
    case("out_add", PCU + g2, gr("-21 -22 -24 -25", 1306, "3 0 3 1 1 1 1", "c rep", "296.", "1.e10",
                                 " 3 2/\n 0/\n", " -1395/\n 3 2/\n 3 1/\n 3 102/\n 0/\n 0/\n"), c511)
    case("out_add1", PC + gr("-21 -22 0 -24", 1306, "3 0 3 1 1 1 1", "c", "296.", "1.e10",
                             " 3 1/\n 3 4/\n 0/\n 0/\n"),
         gr("-21 -22 -24 -25", -1306, "3 0 3 1 1 1 1", "c add", "296.", "1.e10",
            " 3 2/\n 0/\n 0/\n"), c511)
    case("out_lo", PCU + gr("-21 -22 0 -24", 1395, "3 0 3 1 1 1 1", "u", "296.", "1.e10",
                            " 3 1/\n 0/\n 0/\n"),
         gr("-21 -22 -24 -25", 1306, "3 0 3 1 1 1 1", "c lo", "296.", "1.e10", " 3 1/\n 0/\n 0/\n"), c511)
    case("out_all", PCU, gr(U, -1306, "3 0 3 0 1 1 1", "all mats", "296.", "1.e10",
                            " 3 1/\n 3 102/\n 0/\n 3 1/\n 3 102/\n 0/\n"), c511)
    case("out_two", PCU, gr("-21 -22 0 24", 1306, "3 3 3 1 1 1 1", "two mats", "296.", "1.e10",
                            " 3 1/\n 16/\n 13 102/\n 13 51/\n 16 51/\n 0/\n",
                            " 1395/\n 3 1/\n 3 18/\n 3 452/\n 5 18/\n 0/\n 0/\n"), c511)

    # thermal MF6 data from THERMR: H in ZrH (incoherent inelastic and
    # incoherent elastic) and Al (coherent elastic)
    th = ("moder\n 20 -30\nmoder\n 21 -40\nreconr\n -30 -31 /\n 'pendf'/\n %d 0 0 /\n 0.001 /\n 0 /\n"
          "broadr\n -30 -31 -32\n %d 1 /\n 0.001 /\n %s /\n 0 /\n"
          "thermr\n -40 -32 -50 /\n %d %d 20 1 2 1 0 1 %d 0 /\n %s /\n 0.001 1. /\n")
    case("th_zrh", th % (125, 125, "296", 7, 125, 225, "296"),
         gr("-30 -50 0 -24", 125, "3 0 3 3 1 1 1", "h in zrh", "296.", "1.e10",
            " 3 225/\n 3 226/\n 6 225/\n 6 226/\n 0/\n 0/\n"),
         [(20, "n-001_H_001-ENDF8.0-Beta6.endf"), (21, "tsl-HinZrH-ENDF8.0.endf")])
    case("th_al", th % (1325, 1325, "293.6", 53, 1325, 243, "293.6"),
         gr("-30 -50 0 -24", 1325, "3 0 3 3 1 1 1", "al", "293.6", "1.e10",
            " 3 243/\n 3 244/\n 6 243/\n 6 244/\n 0/\n 0/\n"),
         [(20, "n-013_Al_027-ENDF8.0.endf"), (21, "tsl-013_Al_027-ENDF8.0.endf")])

    # U-235 (ENDF/B-V): fission spectra, nubar, delayed data, photons
    PU = pre([1395], "t511")
    case("u5_fis", PU, gr(U, 1395, "3 3 3 1 1 1 1", "u235", "296.", "1.e10",
                          " 3 452/\n 3 455/\n 3 456/\n 5 18/\n 5 455/\n 6 18/\n 6 16/\n 0/\n 0/\n"), c511)
    case("u5_fis0", PU, gr(U, 1395, "3 3 3 1 1 1 0 0", "u235 nosmooth", "296.", "1.e10",
                           " 5 18/\n 5 455/\n 6 18/\n 0/\n 0/\n"), c511)
    case("u5_auto", PU, gr(U, 1395, "3 3 3 1 1 1 1", "u235 auto", "296.", "1.e10",
                           " 6/\n 13/\n 16/\n 0/\n 0/\n"), c511)
    # ENDF/B-VIII U-235 with unresolved data (test 76 deck)
    t76p = (" moder\n 20  -21\n reconr\n -21 -22\n 'pendf tape for  U235'/\n 9228 3 0\n .0005 /\n"
            " 'U 235'/\n 'processed by the njoy nuclear data processing system'/\n"
            " 'see original file for details of evaluation'/\n 0/\n broadr\n -21 -22 -23\n"
            " 9228  1 0 0 0\n .0005 1.e6 /\n 300.\n 0/\n unresr\n -21 -23 -24\n 9228    1 1 1\n"
            " 300.\n 1.e10\n 0/\n")
    case("t76", t76p, gr("-21 -24 0 25", 9228, "9 0 5 5 1 1 1 1", "U 235", "300.", "1.e10",
                         " 3 /\n 3 452 'nu-t' /\n 3 455 'nu-d' /\n 3 456 'nu-p' /\n 5 18 'chi ' /\n"
                         " 5 455 'chid' /\n 6 /\n 0/\n 0/\n"), [(20, "n-092_U_235-ENDF8.0.endf")])

    # File 6 matrices and production
    def m6(name, f, mat, cards, tol=".005", card2="3 3 3 1 1 1 1"):
        case(name, pre([mat], f, tol=tol), gr(U, mat, card2, name, "296.", "1.e10", cards), [(20, f)])
    m6("mf6_be9", "n-004_Be_009-ENDF8.0.endf", 425, " 6/\n 16/\n 21/\n 25/\n 10/\n 0/\n 0/\n")
    m6("mf6_h1", "n-001_H_001-ENDF8.0-Beta6.endf", 125, " 3/\n 6/\n 16/\n 26/\n 0/\n 0/\n")
    m6("mf6_h2", "n-001_H_002-ENDF8.0.endf", 128, " 3/\n 6/\n 16/\n 10/\n 0/\n 0/\n")
    m6("mf6_h2b", "n-001_H_002-ENDF8.0.endf", 128,
       " -1/\n 9 102 1003 0/\n -1/\n 3 102 1003 0/\n -1/\n 6 16 1001 0/\n 3 102/\n 0/\n 0/\n")
    m6("mf6_al", "n-013_Al_027-ENDF8.0.endf", 1325, " 6/\n 10/\n 21/\n 22/\n 23/\n 24/\n 25/\n 26/\n 0/\n 0/\n",
       tol=".02", card2="3 3 3 1 1 1 1")
    m6("mf6_alg", "n-013_Al_027-ENDF8.0.endf", 1325, " 16/\n 0/\n 0/\n", tol=".02")
    m6("mf6_fe", "n-026_Fe_056-JEFF3.3.endf", 2631, " 6/\n 21/\n 25/\n 26/\n 16/\n 0/\n 0/\n", tol=".02")
    m6("mf6_b10", "n-005_B_010-ENDF8.0.endf", 525, " 6/\n 13/\n 16/\n 21/\n 25/\n 31/\n 35/\n 26 51/\n 26 2/\n 0/\n 0/\n")
    m6("mf6_co", "n-027_Co_058-ENDF8.0.endf", 2722, " 21/\n 25/\n 16 102/\n 0/\n 0/\n", tol=".02")
    m6("mf6_ni", "eni61", 2834, " 6/\n 16/\n 21/\n 25/\n 0/\n 0/\n", tol=".02")
    m6("mf6_pu1", "e6pu241c", 9443, " 5 18/\n 6 16/\n 6 18/\n 8/\n 16/\n 0/\n 0/\n", tol=".02")

    # incident charged particles
    def cpre(mat, tol=".005"):
        return "moder\n 20 -21/\nreconr\n -21 -22/\n 'cp pendf'/\n %d 0/\n %s/\n 0/\n" % (mat, tol)
    # (MF3 MT2 for a charged particle makes the Fortran getmf6 write past
    # its 2-word feed-function array; that case is not used here)
    case("cp_pd", cpre(128), gr(U, 128, "3 0 3 3 1 1 1", "p+d", "0.", "1.e10",
                                " 3 1/\n 3 102/\n 21/\n 22/\n 0/\n 0/\n"), [(20, "p-001_H_002-ENDF8.0.endf")])
    case("cp_pd2", cpre(128), gr(U, 128, "3 0 3 3 1 1 1", "p+d", "0.", "1.e10",
                                 " 21 2/\n 22 2/\n 6/\n 0/\n 0/\n"), [(20, "p-001_H_002-ENDF8.0.endf")])
    case("cp_aa", cpre(228), gr(U, 228, "3 0 3 2 1 1 1", "a+a", "0.", "1.e10",
                                " 25 2/\n 0/\n 0/\n"), [(20, "a-002_He_004-ENDF8.0.endf")])

    case("cp_dd", cpre(128), gr(U, 128, "3 0 3 3 1 1 1", "d+d", "0.", "1.e10",
                                " 3 1/\n 6/\n 22/\n 0/\n 0/\n"), [(20, "d-001_H_002-ENDF8.0.endf")])
    case("cp_dpt", cpre(7807, ".02"), gr(U, 7807, "3 0 3 2 1 1 1", "d+pt184", "0.", "1.e10",
                                         " 3 5/\n 6/\n 21/\n 22/\n 26/\n 0/\n 0/\n"),
         [(20, "d-078_Pt_184-TENDL2021.endf")])
    case("cp_pn14", cpre(725, ".02"), gr(U, 725, "3 0 3 2 1 1 1", "p+n14", "0.", "1.e10",
                                         " 3 5/\n 6/\n 21/\n 25/\n 0/\n 0/\n"), [(20, "epn14")])
    # incident photons (photonuclear data): initial groups are gamma groups
    case("g_co59", cpre(2725, ".02"), gr(U, 2725, "3 3 3 1 1 1 1", "g+co59", "0.", "1.e10",
                                         " 3/\n 6/\n 16/\n 21/\n 26/\n 0/\n 0/\n"),
         [(20, "g-027_Co_59-IAEA.endf")])
    case("g_u5", cpre(9228, ".02"), gr(U, 9228, "3 3 3 1 1 1 1", "g+u235", "0.", "1.e10",
                                       " 3/\n 6/\n 16/\n 0/\n 0/\n"), [(20, "g-092_U_235-ENDF8.0.endf")])
    m6("mf6_ca", "n-020_Ca_040-ENDF8.0.endf", 2025, " 6/\n 16/\n 21/\n 25/\n 0/\n 0/\n", tol=".02")

    # error exits
    case("err_mfd", PC, gr(U, 1306, "3 0 3 0 1 1 1", "e", "296.", "1.e10", " 7 1/\n 0/\n 0/\n"), c511)
    case("err_ign", PC, gr(U, 1306, "37 0 3 0 1 1 1", "e", "296.", "1.e10", " 3 1/\n 0/\n 0/\n"), c511)
    case("err_ignm", PC, gr(U, 1306, "-1 0 3 0 1 1 1", "e", "296.", "1.e10", " 3 1/\n 0/\n 0/\n"), c511)
    case("err_igg", PC, gr(U, 1306, "3 11 3 0 1 1 1", "e", "296.", "1.e10", " 3 1/\n 0/\n 0/\n"), c511)
    case("err_iwt", PC, gr(U, 1306, "3 0 13 0 1 1 1", "e", "296.", "1.e10", " 3 1/\n 0/\n 0/\n"), c511)
    case("err_sm", PC, gr(U, 1306, "3 0 3 0 1 1 1 2", "e", "296.", "1.e10", " 3 1/\n 0/\n 0/\n"), c511)
    case("err_g0", PC, gr(U, 1306, "3 0 3 0 1 1 1", "e", "296.", "1.e10", " 16 102/\n 0/\n 0/\n"), c511)
    case("err_ord", PC, gr(U, 1306, "1 0 3 0 1 1 1", "e", "296.", "1.e10",
                           " 3/\n 1.e-5 1.e3 1.e2 2.e7/\n 3 1/\n 0/\n 0/\n"), c511)
    case("err_temp", PC, gr(U, 1306, "3 0 3 0 1 1 1", "e", "900.", "1.e10", " 3 1/\n 0/\n 0/\n"), c511)
    case("err_mt", PC, gr(U, 1306, "3 0 3 0 1 1 1", "e", "296.", "1.e10", " 3 1/\n 3 600/\n 0/\n 0/\n"), c511)
    case("err_ff", PC, gr(U, 1306, "3 0 3 0 1 1 1", "e", "296.", "1.e10", " 3 1/\n 6 102/\n 0/\n 0/\n"), c511)
    case("err_9a", PC, gr(U, 1306, "3 0 3 0 1 1 1", "e", "296.", "1.e10", " -1/\n 4 1 1 0/\n 0/\n 0/\n"), c511)
    print("cases: %d" % n)


if __name__ == "__main__":
    main()
