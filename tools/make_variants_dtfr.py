#!/usr/bin/env python3
"""Make extra DTFR cases from the GROUPR runs of the NJOY2016 test suite.

Only test 03 calls DTFR (on GAMINR output). Here pre.input is the whole
module chain of a test up to and including its GROUPR run (made by
build_corpus.py --module groupr into WORK_DIR), and run.input is DTFR
on that GENDF tape (and on the PENDF tape for the plots).

  gr_*     carbon in graphite (test 01): 30 groups, thermal groups with
           MT229/230, 12 photon groups; edits, plots, CLAW format
  pu_*     Pu-238 (test 02): 50 groups, 3 temperatures, 7 sigma zeros,
           fission chi; temperature and sigma-zero choice
  ni_*     Ni-61 (test 08): 5 Legendre orders, photon tables
  h_*      H-1 (test 30): 175 groups, 42 photon groups, a short table
  u5g_*    U-235 photon production only (test 07): photon tables only
  pho_*    U photoatomic data (ENDF-6) through GAMINR
  u5d_*    U-235 (JENDL-3.3) with delayed neutron data (MT455, MF5)
  err_*    error exits (also an edit MT that is not on the PENDF tape)

Usage: make_variants_dtfr.py NJOY2016_DIR WORK_DIR CORPUS_DIR
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def body(path):
    """A deck without its final 'stop' card."""
    with open(path) as f:
        lines = f.read().splitlines()
    while lines and lines[-1].strip().lower() in ("stop", ""):
        lines.pop()
    return "\n".join(lines) + "\n"


def dtfr(nin, nout, npend, nplot, opts, card3, mats, therm=None, names=None, edits=None,
         photon=(0, 0)):
    s = "dtfr\n %d %d %d %d/\n %s/\n %s/\n" % (nin, nout, npend, nplot, opts, card3)
    if therm:
        s += " %s/\n" % therm
    if names:
        s += " %s/\n" % " ".join("'%s'" % n for n in names)
    if edits:
        s += " %s/\n" % " ".join("%d %d %d" % e for e in edits)
    s += " %d %d/\n" % photon
    for m in mats:
        s += " %s/\n" % m
    return s + " /\n"


def main():
    njoy, work, corpus = sys.argv[1], sys.argv[2], sys.argv[3]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    src = os.path.join(work, "groupr-chains")
    if not os.path.exists(src):
        subprocess.check_call([sys.executable, os.path.join(HERE, "build_corpus.py"), njoy, src,
                               "--module", "groupr"], stdout=subprocess.DEVNULL)
    n = 0

    def case(name, base, run_in, pre_extra=None, tapes=None):
        nonlocal n
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        if base:
            b = os.path.join(src, base)
            for t in sorted(os.listdir(os.path.join(b, "tapes"))):
                p = os.path.join(b, "tapes", t)
                os.symlink(os.path.realpath(p), os.path.join(d, "tapes", t))
            pre = body(os.path.join(b, "pre.input")) + body(os.path.join(b, "run.input"))
        else:
            for u, f in tapes:
                os.symlink(os.path.join(res, f), os.path.join(d, "tapes", "tape%d" % u))
            pre = pre_extra
        with open(os.path.join(d, "pre.input"), "w") as fo:
            fo.write(pre + "stop\n")
        with open(os.path.join(d, "run.input"), "w") as fo:
            fo.write(run_in + "stop\n")
        n += 1

    # carbon in graphite: GENDF 24, PENDF 23, mat 1306, 296 K
    g_edits = [(1, 102, 1), (1, 107, 1), (2, 301, 1), (2, 444, 1)]
    g3 = "4 30 5 7 36 4 4"
    case("gr_tab", "01_1_a", dtfr(24, 61, 0, 0, "0 0 0", g3, ["'c' 1306 1 296."],
                                  "229 230 1", ["nabs", "heat"], g_edits, (2, 12)))
    case("gr_plot1", "01_1_a", dtfr(24, 61, 23, 62, "1 1 0", g3, ["'c' 1306 1 296."],
                                    "229 230 1", ["nabs", "heat"], g_edits, (2, 12)))
    case("gr_plot2", "01_1_a", dtfr(24, 61, 23, 62, "0 2 0", g3, ["'c' 1306 1 296."],
                                    "229 0 0", ["nabs", "heat"], g_edits, (1, 12)))
    case("gr_short", "01_1_a", dtfr(24, 61, 0, 0, "1 0 0", "3 30 3 5 12 0 0", ["'c' 1306 1 296."]))
    case("gr_claw", "01_1_a", dtfr(24, 61, 0, 0, "1 0 1", "4 30", ["'c' 1306 1 296."],
                                   photon=(1, 12)))
    case("gr_bin", "01_1", dtfr(-24, 61, -23, 62, "0 1 0", g3, ["'c' 1306 1 296."],
                                "229 230 1", ["nabs", "heat"], g_edits, (2, 12)))
    case("gr_plot2b", "01_1_a", dtfr(24, 61, 23, 62, "0 2 0", g3,
                                     ["'c' 1306 1 296.", "'c2' 1306 1 296."],
                                     "229 230 1", ["nabs", "heat"], g_edits, (2, 12)))
    case("gr_none", "01_1_a", dtfr(24, 61, 0, 0, "0 0 0", g3,
                                   ["'x' 9999 1 296.", "'c' 1306 1 300.", "'c' 1306 1 296."],
                                   "229 230 1", ["nabs", "heat"], g_edits))

    # Pu-238: GENDF 25, PENDF 24, mat 1050, 300/900/2100 K, 7 sigma zeros
    p_edits = [(1, 18, 1), (2, 102, 1), (3, 470, 1), (1, 300, 0), (4, 16, 1), (4, 17, 2)]
    p3 = "4 50 7 9 58 6 0"
    pn = ["fiss", "capt", "chi", "nxn"]
    case("pu_300", "02_1_a", dtfr(25, 61, 0, 0, "1 0 0", p3, ["'pu238' 1050 1 300."],
                                  names=pn, edits=p_edits))
    case("pu_900", "02_1_a", dtfr(25, 61, 0, 0, "0 0 0", p3, ["'pu238' 1050 3 900."],
                                  names=pn, edits=p_edits))
    case("pu_2100", "02_1_a", dtfr(25, 61, 0, 0, "0 0 0", p3, ["'pu238' 1050 9 2100.",
                                                               "'pu' 1050 1 450.", "'pu' 1050 7"],
                                   names=pn, edits=p_edits))
    case("pu_claw", "02_1_a", dtfr(25, 61, 0, 0, "1 0 1", "4 50", ["'pu238' 1050 2 300."]))
    case("pu_plot", "02_1_a", dtfr(25, 61, 24, 62, "0 1 0", p3, ["'pu238' 1050 1 900."],
                                   names=pn, edits=[e for e in p_edits if e[1] != 300]))

    # Ni-61: GENDF 22, PENDF 24, mat 2834, 300 K, 5 Legendre orders, 12 photon groups
    n_edits = [(1, 2, 1), (2, 4, 1), (2, 51, 1), (2, 91, 1), (3, 102, 1), (4, 103, 1),
               (5, 107, 1), (6, 16, 1)]
    n3 = "5 30 9 11 40 8 0"
    nn = ["elas", "inel", "ngam", "np", "na", "n2n"]
    case("ni_all", "08_1_a", dtfr(22, 61, 24, 62, "1 2 0", n3, ["'ni61' 2834 1 300."],
                                  names=nn, edits=n_edits, photon=(5, 12)))
    case("ni_plot1", "08_1_a", dtfr(22, 61, 24, 62, "0 1 0", n3, ["'ni61' 2834 1 300."],
                                    names=nn, edits=n_edits, photon=(2, 12)))
    m_edits = [(1, 2, 1)] + [(2, m, 1) for m in (4, 51, 52, 53, 54, 55, 56, 57, 58, 91)]
    case("ni_many", "08_1_a", dtfr(22, 61, 24, 62, "0 1 0", "5 30 4 6 41 11 0",
                                   ["'ni61' 2834 1 300."], names=["elas", "inel"],
                                   edits=m_edits))
    case("ni_phot", "08_1_a", dtfr(22, 61, 0, 0, "1 0 0", "0 30 3 4 4 0 0",
                                   ["'ni61' 2834 1 300."], photon=(3, 12)))
    case("ni_claw", "08_1_a", dtfr(22, 61, 0, 0, "0 0 1", "5 30", ["'ni61' 2834 1 300."],
                                   photon=(2, 12)))

    # H-1: GENDF 31, PENDF 23, mat 125, 300 K, 175 groups, 42 photon groups
    case("h_short", "30_1_a", dtfr(31, 61, 23, 62, "0 1 0", "7 175 3 4 20 0 0",
                                   ["'h1' 125 1 300."], photon=(2, 42)))
    case("h_up", "30_1_a", dtfr(31, 61, 0, 0, "1 0 0", "2 175 4 8 30 1 3",
                                ["'h1' 125 1 300."], "222 0 0", ["elas"], [(1, 2, 1)]))

    # U-235 photon production only: GENDF 25, PENDF 24, mat 1395, 300 K
    case("u5g", "07_1_a", dtfr(25, 61, 24, 62, "1 1 0", "0 30 3 4 4 0 0",
                               ["'u235' 1395 1 300."], photon=(1, 94)))
    case("u5g_2", "07_1_a", dtfr(25, 61, 24, 62, "0 2 0", "0 30 3 4 4 0 0",
                                 ["'u235' 1395 1 300."], photon=(1, 94)))

    # U photoatomic data (ENDF-6) through GAMINR
    pre = ("reconr\n 20 21/\n 'u photoatomic'/\n 9200 0 0/\n .001/\n 0/\n"
           "gaminr\n 20 21 0 22/\n 9200 3 3 4 1/\n 'u 12 groups'/\n -1/\n 0/\n")
    t = [(20, "photoat-092_U_000-ENDF8.0.endf")]
    case("pho_u", None, dtfr(22, 61, 21, 62, "1 1 0", "5 12 4 5 16 1 0",
                             ["'u' 9200 1 0."], names=["pheat"], edits=[(1, 621, 1)]), pre, t)
    case("pho_u2", None, dtfr(22, 61, 21, 62, "0 2 0", "5 12 4 5 16 1 0",
                              ["'u' 9200 1 0.", "'u2' 9200 1 0."], names=["pheat"],
                              edits=[(1, 621, 1)]), pre, t)

    # U-235 with delayed neutrons
    pre = ("reconr\n 22 42/\n 'u-235 for dtfr'/\n 9228 0 0/\n .005/\n 0/\n"
           "broadr\n 22 42 32/\n 9228 1 0 0 0/\n .005/\n 300./\n 0/\n"
           "groupr\n 22 32 0 92/\n 9228 3 0 3 1 1 1 0/\n 'u-235'/\n 300./\n 1.e10/\n"
           " 3/\n 3 452/\n 3 455/\n 5 455/\n 6/\n 0/\n 0/\n")
    t = [(22, "J33U235")]
    u_edits = [(1, 18, 1), (2, 102, 1), (3, 470, 1), (4, 471, 1), (5, 455, 1), (5, 452, 1)]
    u3 = "2 30 8 9 38 6 0"
    un = ["fiss", "capt", "chi", "chid", "nu"]
    case("u5d", None, dtfr(92, 61, 32, 62, "1 1 0", u3, ["'u235' 9228 1 300."], names=un,
                           edits=u_edits), pre, t)
    case("u5d_claw", None, dtfr(92, 61, 0, 0, "1 0 1", "2 30", ["'u235' 9228 1 300."]), pre, t)
    case("u5d_nochi", None, dtfr(92, 61, 0, 0, "0 0 0", "2 30 4 5 34 2 0",
                                 ["'u235' 9228 1 300."], names=["fiss"],
                                 edits=[(1, 18, 1), (1, 19, 1)]), pre, t)

    # errors
    case("err_ng", "01_1_a", dtfr(24, 61, 0, 0, "0 0 0", "4 31 5 7 37 0 0", ["'c' 1306 1 296."]))
    case("err_ping", "01_1_a", dtfr(24, 61, 0, 0, "0 0 0", "4 30 5 5 36 0 0", ["'c' 1306 1 296."]))
    # an edit MT that is not on the PENDF tape (error in findf)
    case("err_plot", "02_1_a", dtfr(25, 61, 24, 62, "0 1 0", p3, ["'pu238' 1050 1 900."],
                                    names=pn, edits=p_edits))
    # a PENDF tape (RECONR, 0 K) without the GENDF temperature
    case("err_ptemp", "02_1_a", dtfr(25, 61, 22, 62, "0 1 0", p3, ["'pu238' 1050 1 300."],
                                     names=pn, edits=[e for e in p_edits if e[1] != 300]))
    case("err_ngp", "01_1_a", dtfr(24, 61, 0, 0, "0 0 0", "4 30 3 4 33 0 0", ["'c' 1306 1 296."],
                                   photon=(1, 13)))
    print("%d cases in %s" % (n, corpus))


if __name__ == "__main__":
    main()
