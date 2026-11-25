#!/usr/bin/env python3
"""Make extra THERMR cases with other materials and input options.

Each case runs the Fortran RECONR and BROADR on a neutron ENDF file from
the NJOY2016 test resources (pre.input, loose tolerances to keep it fast)
and then THERMR on the PENDF tape (run.input).

S(alpha,beta) cases (iinc=2) use the thermal files of the test resources
(tape30) with the matching neutron file:

  _sab     E-mu-E' format (iform=0), iprint=1
  _sabf    E-E'-mu format (iform=1)
  _sabp    iprint=2 (full listing), nbin=4

Free-gas and built-in cases (nendf=0, iinc=1 or 0):

  fg_*     free gas, several temperatures, nbin, iform, iprint values
  coh_*    built-in coherent elastic lattices (icoh=1,2,3), with and
           without free-gas inelastic
  iel_*    built-in incoherent elastic (icoh=11,12,13)
  hot      a temperature above 3000 K (the energy grid is scaled)
  err_*    error exits (bad mt, nendf=0 with iinc=2, mode conversion)

Usage: make_variants_thermr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

SAB = [
    # name, neutron file, mat, thermal file, mat, natom, temperature, mtref
    ("h2o", "n-001_H_001-ENDF8.0-Beta6.endf", 125, "tsl-HinH2O-ENDF8.0-Beta6.endf", 1, 2, 300., 222),
    ("al", "n-013_Al_027-ENDF8.0.endf", 1325, "tsl-013_Al_027-ENDF8.0.endf", 53, 1, 400., 243),
    ("dlid", "n-001_H_002-ENDF8.0.endf", 128, "tsl-DinLiD-ENDF8.1-beta.endf", 15, 1, 400., 229),
    ("hzrh", "n-001_H_001-ENDF8.0-Beta6.endf", 125, "tsl-HinZrH-ENDF8.0.endf", 7, 1, 296., 225),
    ("zrzrh", "n-040_Zr_090-ENDF8.0.endf", 4025, "tsl-ZrinZrH-ENDF8.0.endf", 58, 1, 296., 227),
]


def pre(mat, temps):
    t = " ".join("%g" % x for x in temps)
    return ("reconr\n 20 21/\n 'thermr variant'/\n %d 0/\n .01/\n 0/\n"
            "broadr\n 20 21 22/\n %d %d 0 0 0./\n .01/\n %s/\n 0/\nstop\n" % (mat, mat, len(temps), t))


def thermr(nendf, matde, matdp, nbin, temps, iinc, icoh, iform, natom, mtref, iprint,
           tol=0.05, emax=1.2, units=None):
    u = units or (nendf, 22, 23)
    t = " ".join("%g" % x for x in temps)
    return ("thermr\n %d %d %d/\n %d %d %d %d %d %d %d %d %d %d/\n %s/\n %g %g/\nstop\n"
            % (u[0], u[1], u[2], matde, matdp, nbin, len(temps), iinc, icoh, iform, natom, mtref,
               iprint, t, tol, emax))


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    n = 0

    def case(name, tapes, pre_in, run_in):
        nonlocal n
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        for tp, f in tapes.items():
            os.symlink(os.path.join(res, f), os.path.join(d, "tapes", tp))
        with open(os.path.join(d, "pre.input"), "w") as fo:
            fo.write(pre_in)
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(run_in)
        n += 1

    for name, nf, mn, tf, mt, natom, temp, mtref in SAB:
        tapes = {"tape20": nf, "tape30": tf}
        p = pre(mn, [temp])
        case(name + "_sab", tapes, p, thermr(30, mt, mn, 8, [temp], 2, 1, 0, natom, mtref, 1))
        case(name + "_sabf", tapes, p, thermr(30, mt, mn, 8, [temp], 2, 1, 1, natom, mtref, 0))
        case(name + "_sabp", tapes, p, thermr(30, mt, mn, 4, [temp], 2, 1, 0, natom, mtref, 2, tol=0.1))

    h1 = "n-001_H_001-ENDF8.0-Beta6.endf"
    be = "n-004_Be_009-ENDF8.0.endf"
    fe = "n-026_Fe_056-JEFF3.3.endf"
    zr = "n-040_Zr_090-ENDF8.0.endf"
    # free gas
    case("fg_h1_2t", {"tape20": h1}, pre(125, [296., 600.]),
         thermr(0, 0, 125, 8, [296., 600.], 1, 0, 0, 1, 221, 1))
    case("fg_h1_1t", {"tape20": h1}, pre(125, [296., 600.]),
         thermr(0, 0, 125, 16, [600.], 1, 0, 0, 1, 221, 0, tol=0.02, emax=4.))
    case("fg_fe", {"tape20": fe}, pre(2631, [300., 900.]),
         thermr(0, 0, 2631, 8, [300., 900.], 1, 0, 0, 1, 221, 1, emax=2.))
    case("fg_iform1", {"tape20": h1}, pre(125, [296.]),
         thermr(0, 0, 125, 8, [296.], 1, 0, 1, 1, 221, 1))
    case("fg_iform1p", {"tape20": h1}, pre(125, [296.]),
         thermr(0, 0, 125, 4, [296.], 1, 0, 1, 1, 221, 2, tol=0.1, emax=0.5))
    case("fg_print2", {"tape20": be}, pre(425, [296.]),
         thermr(0, 0, 425, 12, [296.], 1, 0, 0, 1, 221, 2, tol=0.1))
    case("hot", {"tape20": h1}, pre(125, [3500.]),
         thermr(0, 0, 125, 8, [3500.], 1, 0, 0, 1, 221, 1))
    # built-in coherent elastic (graphite, Be, BeO) on the Be data
    case("coh_gr", {"tape20": be}, pre(425, [296., 600.]),
         thermr(0, 0, 425, 8, [296., 600.], 1, 1, 0, 1, 229, 1))
    case("coh_be", {"tape20": be}, pre(425, [400.]),
         thermr(0, 0, 425, 8, [400.], 1, 2, 0, 1, 231, 1))
    case("coh_beo", {"tape20": be}, pre(425, [296.]),
         thermr(0, 0, 425, 8, [296.], 0, 3, 0, 2, 233, 1))
    case("coh_only", {"tape20": be}, pre(425, [296.]),
         thermr(0, 0, 425, 8, [296.], 0, 2, 0, 1, 231, 0))
    # built-in incoherent elastic (H in ZrH, Zr in ZrH, polyethylene)
    case("iel_11", {"tape20": h1}, pre(125, [296.]),
         thermr(0, 0, 125, 8, [296.], 1, 11, 0, 1, 225, 1))
    case("iel_12", {"tape20": zr}, pre(4025, [400.]),
         thermr(0, 0, 4025, 8, [400.], 1, 12, 0, 1, 227, 2, tol=0.1))
    case("iel_13", {"tape20": zr}, pre(4025, [600.]),
         thermr(0, 0, 4025, 8, [600.], 1, 13, 0, 1, 227, 1))
    # error exits
    case("err_mt", {"tape20": h1}, pre(125, [296.]),
         thermr(0, 0, 125, 8, [296.], 1, 0, 0, 1, 260, 1))
    case("err_nendf", {"tape20": h1}, pre(125, [296.]),
         thermr(0, 0, 125, 8, [296.], 2, 0, 0, 1, 221, 1))
    case("err_mode", {"tape20": h1}, pre(125, [296.]),
         thermr(0, 0, 125, 8, [296.], 1, 0, 0, 1, 221, 1, units=(0, 22, -23)))
    case("err_temp", {"tape20": h1}, pre(125, [296.]),
         thermr(0, 0, 125, 8, [500.], 1, 0, 0, 1, 221, 1))
    print("variants:", n)


if __name__ == "__main__":
    main()
