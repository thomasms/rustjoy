#!/usr/bin/env python3
"""Make extra HEATR cases with other materials and input options.

Each case runs the Fortran RECONR on an ENDF file from the NJOY2016 test
resources (pre.input, a loose tolerance to keep it fast) and then HEATR on
the PENDF tape with one of these option sets (run.input):

  A  14 partial kermas incl. damage (302 303 304 318 401 402 403 407 442 443
     444 445 446 447), iprint=1 (no plot file: without iprint=2 the plot
     reads columns that are never set, see README)
  B  kinematic check (iprint=2) with 303 304 402 443 444 447, plot file
  C  user q values: mt102 with an energy-dependent q (qa >= 99e6 and a
     TAB1 card), mt16 with a constant q; 301 in the kerma list
  D  local photon deposition, given damage energy (35 eV), ntemp=1
  E  defaults (only the material number)

Every material gets two option sets, rotated over the list, so that each
set meets several kinds of data (File 6, File 12/13 photons, fission, old
ENDF formats).

Usage: make_variants_heatr.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

MATERIALS = [
    ("n-001_H_001-ENDF8.0-Beta6.endf", 125),
    ("n-001_H_002-ENDF8.0.endf", 128),
    ("n-004_Be_009-ENDF8.0.endf", 425),
    ("n-005_B_010-ENDF8.0.endf", 525),
    ("n-013_Al_027-ENDF8.0.endf", 1325),
    ("n-018_Ar_37-tendl2023.endf", 1828),
    ("n-020_Ca_040-ENDF8.0.endf", 2025),
    ("n-025_Mn_055-ENDF8.0.endf", 2525),
    ("n-026_Fe_000-IRDFF-II.endf", 2600),
    ("n-026_Fe_056-JEFF3.3.endf", 2631),
    ("n-026_Fe_056-TENDL19.endf", 2631),
    ("n-026_Fe_058-ENDF8.0-Beta4.endf", 2637),
    ("n-027_Co_058-ENDF8.0.endf", 2722),
    ("n-027_Co_058m1-ENDF8.0.endf", 2723),
    ("n-028_Ni_063-ENDF8.0.endf", 2840),
    ("n-030_Zn_067-ENDF8.0.endf", 3034),
    ("n-032_Ge_071-ENDF8.0.endf", 3228),
    ("n-036_Kr_083-ENDF8.0.endf", 3640),
    ("n-038_Sr_088-ENDF8.1.endf", 3837),
    ("n-040_Zr_090-ENDF8.0.endf", 4025),
    ("n-042_Mo_095-beta.endf", 4234),
    ("n-047_Ag_109-ENDF8.0.endf", 4731),
    ("n-050_Sn_119-ENDF8.0.endf", 5046),
    ("n-072_Hf_177-TENDL2021.endf", 7234),
    ("n-078_Pt_190-ENDF8.1.endf", 7825),
    ("n-082_Pb_208-TENDL2021.endf", 8237),
    ("n-094_Pu_241-ENDF8.0.endf", 9443),
    ("n-094_Pu_245-ENDF8.0-Beta6.endf", 9455),
    ("n-095_Am_242-ENDF8.0.endf", 9546),
    ("n-095_Am_242m1-ENDF8.0.endf", 9547),
    ("t404", 1043),
    ("t511", 1146),
    ("eni61", 2834),
    ("epn14", 725),
    ("e6pu241c", 9443),
    ("cl35rml", 1725),
    ("J33U235", 9228),
    ("J33U238", 9237),
    ("J33Pu239", 9437),
    # charged-particle and light-ion incident data (File 6 discrete lines,
    # other kinematics)
    ("p-001_H_001-ENDF8.0.endf", 125),
    ("p-001_H_002-ENDF8.0.endf", 128),
    ("p-001_H_003-ENDF8.0.endf", 131),
    ("d-001_H_002-ENDF8.0.endf", 128),
    ("d-002_He_003.endf", 225),
    ("a-002_He_004-ENDF8.0.endf", 228),
    ("d-078_Pt_184-TENDL2021.endf", 7807),
]

QBAR = " 0. 0. 0 0 1 3 3 2 1.e-5 5.0e6 1.e6 5.5e6 2.e7 6.0e6/"


def deck(opt, mat):
    if opt == "A":
        return ("heatr\n 20 21 22/\n %d 14 0 0 0 1/\n"
                " 302 303 304 318 401 402 403 407 442 443 444 445 446 447/\n" % mat)
    if opt == "B":
        return "heatr\n 20 21 22 31/\n %d 6 0 0 0 2/\n 303 304 402 443 444 447/\n" % mat
    if opt == "C":
        return ("heatr\n 20 21 22/\n %d 3 2 0 0 1/\n 301 444 445/\n 102 16/\n"
                " 1.0e8 -1.0e7/\n%s\n" % (mat, QBAR))
    if opt == "D":
        return "heatr\n 20 21 22/\n %d 2 0 1 1 0 35./\n 302 444/\n" % mat
    return "heatr\n 20 21 22/\n %d/\n" % mat


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    opts = "ABCDE"
    n = 0
    for k, (f, mat) in enumerate(MATERIALS):
        for o in (opts[k % 5], opts[(k + 2) % 5]):
            name = "%s_%s" % (f.split("-ENDF")[0].split(".endf")[0].replace(".", "_"), o)
            d = os.path.join(corpus, name)
            if os.path.exists(d):
                sys.exit("%s exists (nothing is deleted)" % d)
            os.makedirs(os.path.join(d, "tapes"))
            os.symlink(os.path.join(res, f), os.path.join(d, "tapes", "tape20"))
            with open(os.path.join(d, "pre.input"), "w") as fo:
                fo.write("reconr\n 20 21/\n 'heatr variant'/\n %d 0/\n .01/\n 0/\nstop\n" % mat)
            with open(os.path.join(d, "run.input"), "w") as fo:
                fo.write(deck(o, mat) + "stop\n")
            n += 1
    print("variants:", n)


if __name__ == "__main__":
    main()
