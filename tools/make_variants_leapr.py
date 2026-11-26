#!/usr/bin/env python3
"""Make extra LEAPR cases with other input options.

LEAPR needs no input tapes: each case is one leapr deck (run.input) and an
empty pre.input. The decks use small alpha and beta grids to keep them
fast. The frequency spectra are the ones of NJOY tests 09 (H in H2O) and
23 (Be in BeO); the s(kappa) table is a smooth made-up curve.

  h2o_*    H in H2O: free gas, diffusion, oscillators, print levels,
           several temperatures (with reuse of parameters), lat=1,
           isabt=1, ilog=1, a pure solid (incoherent elastic), a large
           phonon order, secondary scatterers (free and mixed), 320 alpha
           values, large alpha values, and a low-energy oscillator
  coh_*    coherent elastic for every lattice (iel 1..6), 1 to 3 temps
  cold_*   ortho and para hydrogen and deuterium (ncold 1..4), with
           isabt and ilog, and two temperatures with a wide alpha grid
  sk_*     s(kappa) options (nsk 1 and 2)
  err_*    cases that stop with an error

Usage: make_variants_leapr.py CORPUS_DIR
"""
import math
import os
import sys

H2O = """0 .0005 .001 .002 .0035 .005 .0075 .01 .013 .0165 .02 .0245
 .029 .034 .0395 .045 .0506 .0562 .0622 .0686 .075 .083 .091
 .099 .107 .115 .1197 .1214 .1218 .1195 .1125 .1065 .1005 .09542
 .09126 .0871 .0839 .0807 .07798 .07574 .0735 .07162 .06974
 .06804 .06652 .065 .0634 .0618 .06022 .05866 .0571 .05586
 .05462 .0535 .0525 .0515 .05042 .04934 .04822 .04706 .0459
 .04478 .04366 .04288 .04244 .042 0./"""

BE = """0.0 .3 .7 .9 1. 1.2 1.6 2.0 2.2 3.0 3.5 4.5 5.5 6.8 8.0 9.2 10.9 12.9
15.5 18.6 22.0 26.0 30.5 35.0 39.0 40.0 34.0 28.0 26.0 24.4 23.0
21.3 19.8 17.0 14.1 12.0 10.0 9.0 9.0 8.5 7.5 6.0 4.6 3.1 1.6
0.5 0. 0.0 4.0 15.0 38.0 52.0 70.0 105.0 165.0 230.0 200.0 170.0
145.0 136.0 134.0 112.0 96.0 89.0 84.0 75.0 87.0 81.0 66.0 59.0
68.0 105.0 95.0 97.0 135.0 163.0 130.0 111.0 92.0 67.0 45.0 19.0
7.0 0.0/"""

ALPHA = ".01 .02 .05 .1 .2 .35 .5 .75 1. 1.5 2. 3. 4. 6. 8. 10. 15. 20. 30. 40./"
BETA = ("0. .05 .1 .2 .3 .5 .75 1. 1.25 1.5 2. 2.5 3. 3.5 4. 5. 6. 7. 8. 10. "
        "12. 15. 18. 20. 25./")
ALPHA_C = ".005 .01 .05 .1 .3 .7 1. 2. 4. 8. 15. 30./"
BETA_C = "0. .05 .2 .5 1. 2. 4. 5.8 6. 8.7 10. 15. 20. 30. 40./"
COMMENTS = " ' leapr variant for testing'/\n ' '/\n ' a second comment line'/\n /\n"


def ska(n, dk):
    vals = []
    for k in range(1, n + 1):
        x = k * dk
        vals.append("%.5e" % (1.0 - 0.9 * math.exp(-x) + 0.3 * math.exp(-(x - 2.5) ** 2)))
    lines = [" ".join(vals[i:i + 5]) for i in range(0, len(vals), 5)]
    return "\n".join(lines) + "/"


def head(nout, title, ntempr, iprint, nphon, mat, za, isabt, ilog,
         awr, spr, npr, iel, ncold, nsk, nss="0", na=20, nb=25, lat=0,
         alpha=ALPHA, beta=BETA):
    return ("leapr\n %d/\n '%s'/\n %d %d %d/\n %d %d %d %d 1e-75/\n"
            " %s %s %d %d %d %d/\n %s/\n %d %d %d/\n %s\n %s\n"
            % (nout, title, ntempr, iprint, nphon, mat, za, isabt, ilog,
               awr, spr, npr, iel, ncold, nsk, nss, na, nb, lat, alpha, beta))


def h2o_temp(t, twt="0.055556", c="0.", tbeta="0.444444", nd=2,
             osc=(".205 .48", ".166667 .333333"), extra=""):
    s = " %s/\n .00255 67/\n %s\n %s %s %s/\n %d/\n" % (t, H2O, twt, c, tbeta, nd)
    if nd > 0:
        s += " %s/\n %s/\n" % osc
    return s + extra


def be_temp(t):
    return " %s/\n .0016518 84/\n %s\n 0. 0. 1./\n 0/\n" % (t, BE)


def main():
    corpus = sys.argv[1]
    n = 0

    def case(name, deck):
        nonlocal n
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        with open(os.path.join(d, "pre.input"), "w") as fo:
            fo.write("stop\n")
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(deck + "stop\n")
        n += 1

    h = lambda **k: head(24, "h in h2o variant", **{**dict(
        ntempr=1, iprint=1, nphon=60, mat=101, za=1001, isabt=0, ilog=0,
        awr="0.99917", spr="20.449", npr=2, iel=0, ncold=0, nsk=0), **k})

    # H in H2O
    case("h2o_p0", h(iprint=0) + h2o_temp(296) + COMMENTS)
    case("h2o_p1", h() + h2o_temp(296) + COMMENTS)
    case("h2o_p2", h(iprint=2, nphon=40) + h2o_temp(296) + COMMENTS)
    case("h2o_t3", h(ntempr=3) + h2o_temp(296) + " -350/\n" + h2o_temp(400) + COMMENTS)
    case("h2o_lat", h(ntempr=2, lat=1) + h2o_temp(296) + h2o_temp(600) + COMMENTS)
    case("h2o_diff", h(iprint=2, nphon=40) + h2o_temp(296, c="120.") + COMMENTS)
    case("h2o_diff1", h() + h2o_temp(296, c="12.") + COMMENTS)
    case("h2o_sym2", h(isabt=1) + h2o_temp(296) + COMMENTS)
    case("h2o_log", h(ilog=1) + h2o_temp(296) + COMMENTS)
    case("h2o_log2", h(ilog=1, isabt=1, ntempr=2) + h2o_temp(296) + " -400/\n" + COMMENTS)
    case("h2o_solid", h(iprint=2, nphon=40) + h2o_temp(296, twt="0.", tbeta="0.5") + COMMENTS)
    case("h2o_solid2", h(ntempr=2) + h2o_temp(296, twt="0.", tbeta="0.5")
         + h2o_temp(500, twt="0.", tbeta="0.5") + COMMENTS)
    case("h2o_nod", h(ntempr=1) + h2o_temp(296, twt="0.", tbeta="1.", nd=0) + COMMENTS)
    case("h2o_nod2", h(ntempr=2, iprint=2, nphon=30) + h2o_temp(296, twt="0.", tbeta="1.", nd=0)
         + " -800/\n" + COMMENTS)
    case("h2o_gas", h() + h2o_temp(296, twt="0.5", tbeta="0.5", nd=0) + COMMENTS)
    case("h2o_osc1", h(iprint=2, nphon=30) + h2o_temp(296, nd=1, tbeta="0.777777",
                                                     osc=(".48", ".166667")) + COMMENTS)
    case("h2o_nphon", h(iprint=0, nphon=260) + h2o_temp(296, nd=0, tbeta="0.944444") + COMMENTS)
    case("h2o_free2", h(nss=" 1 1. 15.858 3.7481 1") + h2o_temp(296) + COMMENTS)
    case("h2o_free2d", h(nss=" 1 2. 15.858 3.7481 1") + h2o_temp(296) + COMMENTS)
    case("h2o_mix", h(ntempr=2, nss=" 1 0. 15.858 3.7481 1") + h2o_temp(296) + " -400/\n"
         + h2o_temp(296, twt="0.", tbeta="1.", nd=0) + " -400/\n" + COMMENTS)
    case("h2o_mix_coh", h(iel=1, nss=" 1 0. 15.858 3.7481 1") + h2o_temp(296)
         + h2o_temp(296, twt="0.", tbeta="1.", nd=0) + COMMENTS)

    # coherent elastic, all lattices
    lats = {1: ("11.898", "4.7392", 1), 2: ("8.93478", "6.15", 1), 3: ("8.93478", "6.15", 1),
            4: ("26.7495", "1.495", 1), 5: ("205.2", "11.115", 1), 6: ("55.454", "11.22", 1)}
    for iel, (awr, spr, npr) in lats.items():
        hc = lambda **k: head(24, "coherent %d" % iel, **{**dict(
            ntempr=1, iprint=1, nphon=40, mat=200 + iel, za=4009, isabt=0, ilog=0,
            awr=awr, spr=spr, npr=npr, iel=iel, ncold=0, nsk=0), **k})
        case("coh_%d" % iel, hc() + be_temp(296) + COMMENTS)
        case("coh_%d_t2" % iel, hc(ntempr=2) + be_temp(296) + " -600/\n" + COMMENTS)
    case("coh_3_t3", head(24, "coherent beo 3t", 3, 0, 30, 203, 4009, 0, 0, "8.93478", "6.15", 2,
                           3, 0, 0) + be_temp(296) + be_temp(500) + " -900/\n" + COMMENTS)

    # cold hydrogen and deuterium
    kap = " 120 .05/\n" + ska(120, 0.05) + "\n"
    for ncold in (1, 2, 3, 4):
        awr, spr, za = ("0.99917", "20.478", 1001) if ncold <= 2 else ("1.9968", "3.395", 1002)
        for tag, isabt, ilog, ipr in (("", 0, 0, 1), ("_s", 1, 0, 2), ("_l", 0, 1, 1), ("_sl", 1, 1, 1)):
            deck = head(20, "cold %d" % ncold, 1, ipr, 30, 300 + ncold, za, isabt, ilog, awr, spr, 2,
                        0, ncold, 0, na=12, nb=15, alpha=ALPHA_C, beta=BETA_C)
            deck += (" 20/\n .00025 48/\n 0. .01563 .0625 .141 .25 .391 .5625 .766 1. 1.266"
                     " 1.5625 1.89\n 2.25 2.64 3.0625 3.52 4. 4.6 5.5 7.0 8.5 9.2 9.5 9.4\n"
                     " 9.2 8.9 8.5 8.0 7.5 7.05 6.7 6.4 6.2\n 6.1 6.2 6.45 6.7 6.95 7.1 6.55"
                     " 5.5 3. 0.\n 0. 0. 0. 0. 0./\n .025 40 .475/\n 0/\n") + kap
            case("cold_%d%s" % (ncold, tag), deck + COMMENTS)
    deck = head(20, "cold 2 lat", 2, 1, 30, 302, 1001, 0, 0, "0.99917", "20.478", 2,
                0, 2, 1, na=12, nb=15, lat=1, alpha=ALPHA_C, beta=BETA_C)
    temp = (" %s/\n .00025 48/\n 0. .01563 .0625 .141 .25 .391 .5625 .766 1. 1.266"
            " 1.5625 1.89\n 2.25 2.64 3.0625 3.52 4. 4.6 5.5 7.0 8.5 9.2 9.5 9.4\n"
            " 9.2 8.9 8.5 8.0 7.5 7.05 6.7 6.4 6.2\n 6.1 6.2 6.45 6.7 6.95 7.1 6.55"
            " 5.5 3. 0.\n 0. 0. 0. 0. 0./\n .025 0. .475/\n 0/\n") + kap + " .3/\n"
    case("cold_2_lat", deck + temp % "20" + temp % "25" + COMMENTS)
    # two temperatures (LIST records), every isabt/ilog pair
    alpha_w = "1e-7 1e-5 .001 .1 1. 10. 100. 500. 2000. 5000./"
    for tag, isabt, ilog in (("", 0, 0), ("_s", 1, 0), ("_l", 0, 1), ("_sl", 1, 1)):
        deck = head(20, "cold 1 2t", 2, 1, 30, 301, 1001, isabt, ilog, "0.99917", "20.478", 2,
                    0, 1, 0, na=10, nb=15, alpha=alpha_w, beta=BETA_C)
        case("cold_1_t2%s" % tag, deck + (temp.replace(" .3/\n", "") % "20") + " -25/\n" + COMMENTS)

    # more H in H2O cases for the LIST records and the rarer branches
    case("h2o_t2_l", h(ntempr=2, ilog=1) + h2o_temp(296) + " -500/\n" + COMMENTS)
    case("h2o_t2_s", h(ntempr=2, isabt=1) + h2o_temp(296) + " -500/\n" + COMMENTS)
    a320 = " ".join("%.6g" % (0.01 * 1.02 ** k) for k in range(320)) + "/"
    case("h2o_na320", h(ntempr=2, nphon=30, alpha=a320, na=320, nb=6, beta="0. .1 .5 1. 2. 5./")
         + h2o_temp(296) + " -500/\n" + COMMENTS)
    case("h2o_bigal", h(alpha=".01 .1 1. 10. 50. 100. 200. 400. 800./", na=9,
                        beta=".02 .1 .3 .6 1. 2. 4. 8. 15./", nb=9)
         + h2o_temp(296) + COMMENTS)
    case("h2o_osc0", h(iprint=2, nphon=30) + h2o_temp(296, twt="0.", tbeta="0.9", nd=1,
                                                     osc=(".0005", ".1")) + COMMENTS)

    # s(kappa) options
    kap2 = " 80 .1/\n" + ska(80, 0.1) + "\n .4/\n"
    case("sk_1", h(nsk=1) + h2o_temp(296, extra=kap2) + COMMENTS)
    case("sk_2", h(nsk=2) + h2o_temp(296, extra=kap2) + COMMENTS)
    case("sk_2p", h(nsk=2, iprint=2, nphon=30, ntempr=2) + h2o_temp(296, extra=kap2) + " -500/\n"
         + COMMENTS)
    case("sk_2lat", h(nsk=2, lat=1) + h2o_temp(296, extra=kap2) + COMMENTS)

    # errors
    case("err_iel", h(iel=7) + h2o_temp(296) + COMMENTS)
    case("err_eof", h() + h2o_temp(296))
    print("%d cases in %s" % (n, corpus))


if __name__ == "__main__":
    main()
