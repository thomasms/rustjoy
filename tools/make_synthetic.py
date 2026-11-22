#!/usr/bin/env python3
"""Make synthetic RECONR cases for resonance formats that no evaluation in
the NJOY2016 test suite uses: Adler-Adler (LRF=4) and hybrid R-function
(LRF=6).

The MF2/MT151 section of a real ENDF material (Pu-241, ENDF/B-VIII.0) is
replaced by made-up parameters. The values are not physical; they only make
RECONR run the Adler-Adler and hybrid code (and their temperature branches)
so that the Fortran and the Rust results can be compared.

The hybrid data avoid the paths where the Fortran reads undefined values:
no competing channels (NIRE=NCRE=0) and no tabulated background R-function
or phase shift (LBK=LPS=0), because cshybr then uses an unset `eeff`.

Usage: make_synthetic.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

BASE = "n-094_Pu_241-ENDF8.0.endf"
MAT = 9443


def f11(x):
    """ENDF 11-column real, like NJOY's a11 for ordinary values."""
    if x == 0:
        return " 0.000000+0"
    s = "%.6e" % x  # d.dddddde+xx
    m, e = s.split("e")
    e = int(e)
    return ("%s%s%d" % (m, "+" if e >= 0 else "-", abs(e))).rjust(11)


def i11(n):
    return ("%d" % n).rjust(11)


class Sec:
    def __init__(self):
        self.lines = []

    def card(self, fields):
        self.lines.append("".join(fields))

    def cont(self, c1, c2, l1, l2, n1, n2):
        self.card([f11(c1), f11(c2), i11(l1), i11(l2), i11(n1), i11(n2)])

    def list(self, c1, c2, l1, l2, n2, data):
        self.cont(c1, c2, l1, l2, len(data), n2)
        for k in range(0, len(data), 6):
            row = [f11(v) for v in data[k:k + 6]]
            row += [" " * 11] * (6 - len(row))
            self.card(row)


def head(s, awr, el, eh, lrf):
    s.cont(94241.0, awr, 0, 0, 1, 0)
    s.cont(94241.0, 1.0, 0, 0, 1, 0)
    s.cont(el, eh, 1, lrf, 0, 0)


def adler(li):
    awr = 238.978
    s = Sec()
    head(s, awr, 1.0e-5, 100.0, 4)
    s.cont(0.5, 0.95, 0, 0, 2, 0)   # SPI, AP, NLS=2
    bkg = [2.0e-9, 1.0e-10, 0.0, 0.0, 1.0e-12, 0.0,      # total
           1.0e-9, 0.0, 1.0e-11, 0.0, 0.0, 1.0e-14,      # fission
           5.0e-10, 2.0e-10, 0.0, 1.0e-13, 0.0, 0.0]     # capture
    s.list(awr, 0.0, li, 0, 3, bkg)
    res = {  # (l, J): [(E, half width, GRT, GIT, GRF, GIF, GRC, GIC)]
        (0, 0.0): [(-0.45, 0.030, 2.4e-4, 1.0e-5, 0.8e-4, 2.0e-6, 0.7e-4, -1.0e-6),
                   (0.264, 0.020, 2.1e-4, 2.0e-5, 0.9e-4, 3.0e-6, 0.6e-4, 1.0e-6),
                   (4.28, 0.045, 3.0e-4, -1.0e-5, 1.2e-4, 1.0e-6, 0.9e-4, 2.0e-6),
                   (14.77, 0.060, 4.5e-4, 3.0e-5, 1.5e-4, -2.0e-6, 1.1e-4, 1.0e-6),
                   (52.3, 0.090, 6.0e-4, 1.0e-5, 2.0e-4, 1.0e-6, 1.6e-4, 0.0)],
        (0, 1.0): [(1.13, 0.025, 2.6e-4, 0.0, 0.7e-4, 1.0e-6, 0.8e-4, 1.0e-6),
                   (9.21, 0.050, 3.4e-4, 2.0e-5, 1.1e-4, 0.0, 1.0e-4, -1.0e-6),
                   (88.6, 0.110, 7.0e-4, 0.0, 2.2e-4, 3.0e-6, 1.9e-4, 2.0e-6)],
        (1, 1.5): [(22.5, 0.070, 1.5e-4, 1.0e-5, 0.5e-4, 0.0, 0.4e-4, 0.0),
                   (61.0, 0.100, 2.5e-4, -1.0e-5, 0.6e-4, 1.0e-6, 0.7e-4, 1.0e-6)],
    }
    for l in (0, 1):
        groups = [(j, v) for (ll, j), v in sorted(res.items()) if ll == l]
        s.cont(0.0, 0.0, l, 0, len(groups), 0)
        for j, rr in groups:
            data = []
            for (e, w, grt, git, grf, gif, grc, gic) in rr:
                data += [e, w, grt, git, e, w, grf, gif, e, w, grc, gic]
            s.list(j, 0.0, 0, 0, len(rr), data)
    return s


def hybrid():
    awr = 238.978
    s = Sec()
    head(s, awr, 1.0e-5, 100.0, 6)
    s.cont(0.5, 0.95, 0, 0, 3, 0)          # SPI, AP, LAD=0, NLS=3
    s.cont(0.0, 0.0, 1, 1, 0, 0)           # NGRE=1, NFRE=1, NIRE=NCRE=0
    s.cont(0.0, 0.0, 0, 0, 0, 0)
    s.list(0.0, 0.0, 0, 0, 0, [0.0, 0.0, 0.0, 0.0])
    groups = {  # l: [(J, AC, [(ER, GN, GG, GF)])]
        0: [(0.0, 0.95, [(-0.6, 2.0e-3, 0.040, 0.10), (0.30, 8.0e-5, 0.042, 0.12),
                         (3.9, 4.0e-4, 0.038, 0.31), (17.2, 1.5e-3, 0.041, 0.05)]),
            (1.0, 0.95, [(1.27, 2.0e-4, 0.040, 0.20), (8.6, 9.0e-4, 0.039, 0.02),
                         (44.0, 3.0e-3, 0.040, 0.45)])],
        1: [(0.5, 0.90, [(12.3, 5.0e-5, 0.040, 0.08), (71.0, 4.0e-4, 0.040, 0.15)]),
            (1.5, 0.90, [(27.4, 1.0e-4, 0.041, 0.01)])],
        4: [(3.5, 0.85, [(36.1, 2.0e-6, 0.040, 0.0), (93.0, 1.0e-5, 0.040, 0.03)])],
    }
    for l in (0, 1, 4):
        s.cont(awr, 0.0, l, 0, 1, 0)        # one spin sequence (NSS=1)
        s.cont(0.0, 0.0, 0, 0, len(groups[l]), 0)
        for j, ac, rr in groups[l]:
            data = []
            for (er, gn, gg, gf) in rr:
                data += [er, gn, gg, gf] + [0.0] * 8
            s.list(j, ac, 0, 0, len(rr), data)
    return s


def write_tape(src, dst, sec):
    out = []
    done = False
    for line in open(src):
        line = line.rstrip("\n")
        mat, mf, mt = line[66:70], line[70:72], line[72:75]
        if mf.strip() == "2" and mt.strip() == "151" and mat.strip() == str(MAT):
            if not done:
                for k, l in enumerate(sec.lines, 1):
                    out.append("%-66s%4d%2d%3d%5d" % (l, MAT, 2, 151, k % 100000))
                done = True
            continue
        out.append(line)
    assert done
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


CASES = {
    # name: (section, err, tempr)
    "syn_aa7": (lambda: adler(7), 0.001, 0.0),
    "syn_aa7_t": (lambda: adler(7), 0.001, 293.6),
    "syn_aa5": (lambda: adler(5), 0.002, 0.0),
    "syn_aa6_t": (lambda: adler(6), 0.002, 1200.0),
    "syn_hy": (hybrid, 0.001, 0.0),
    "syn_hy_t": (hybrid, 0.001, 293.6),   # Fortran stops: no Doppler for hybrid
}


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
    src = os.path.join(njoy, "tests", "resources", BASE)
    for name, (mk, err, temp) in CASES.items():
        d = os.path.join(corpus, name)
        if os.path.exists(d):
            sys.exit("%s exists (nothing is deleted)" % d)
        os.makedirs(os.path.join(d, "tapes"))
        write_tape(src, os.path.join(d, "tapes", "tape20"), mk())
        with open(os.path.join(d, "pre.input"), "w") as f:
            f.write("stop\n")
        with open(os.path.join(d, "run.input"), "w") as f:
            f.write("reconr\n 20 98/\n 'synthetic %s'/\n %d 1/\n %g %g/\n"
                    " 'synthetic resonance data, not physical'/\n 0/\nstop\n"
                    % (name, MAT, err, temp))
    print("synthetic cases:", len(CASES))


if __name__ == "__main__":
    main()
