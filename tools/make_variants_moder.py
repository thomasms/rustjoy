#!/usr/bin/env python3
"""Make extra MODER cases from the ENDF files of the NJOY2016 test resources.

For every resource file (incident neutrons, photons, charged particles,
thermal scattering, photo-atomic, atomic relaxation, decay and covariance
data) two cases are made:

  <file>_tb   text -> binary   (run: moder 20 -21)
  <file>_bt   binary -> text   (pre: moder 20 -21 by the Fortran;
                                run: moder -21 22)

and a few material-selection (loop mode, nin=1) cases:

  sel_tt      two text tapes, one material each, onto a text tape
  sel_tb      the same onto a binary tape
  sel_bb      binary input tapes (made by the Fortran moder in pre.input)
  sel_same    two materials from the same tape in one call
  sel_order   materials in descending order (the error exit)
  sel_blank   a blank tape id

Usage: make_variants_moder.py NJOY2016_DIR CORPUS_DIR
"""
import os
import sys

SKIP = {"thermalscattering.ace"}


def mat_of(path):
    with open(path, errors="replace") as f:
        f.readline()
        line = f.readline()
    return int(line[66:70])


def case(corpus, name, tapes, pre, run):
    d = os.path.join(corpus, name)
    if os.path.exists(d):
        sys.exit("%s exists (nothing is deleted)" % d)
    os.makedirs(os.path.join(d, "tapes"))
    for t, src in tapes.items():
        os.symlink(src, os.path.join(d, "tapes", t))
    with open(os.path.join(d, "pre.input"), "w") as f:
        f.write(pre + "stop\n")
    with open(os.path.join(d, "run.input"), "w") as f:
        f.write(run + "stop\n")


def main():
    njoy, corpus = sys.argv[1], sys.argv[2]
    res = os.path.join(os.path.abspath(njoy), "tests", "resources")
    n = 0
    for f in sorted(os.listdir(res)):
        if f in SKIP:
            continue
        src = os.path.join(res, f)
        base = f.replace(".endf", "").replace(".", "_")
        case(corpus, base + "_tb", {"tape20": src}, "", "moder\n 20 -21/\n")
        case(corpus, base + "_bt", {"tape20": src}, "moder\n 20 -21/\n", "moder\n -21 22/\n")
        n += 2

    # material selection (loop mode)
    a = os.path.join(res, "n-001_H_002-ENDF8.0.endf")      # mat 128
    b = os.path.join(res, "n-026_Fe_000-IRDFF-II.endf")    # mat 2600
    c = os.path.join(res, "n-050_Sn_119-ENDF8.0.endf")     # mat 5046
    ma, mb, mc = mat_of(a), mat_of(b), mat_of(c)
    two = {"tape20": a, "tape21": b, "tape24": c}
    case(corpus, "sel_tt", two, "",
         "moder\n 1 22/\n 'selected materials'/\n 20 %d/\n 21 %d/\n 0/\n" % (ma, mb))
    case(corpus, "sel_tb", two, "",
         "moder\n 1 -22/\n 'selected materials, binary'/\n 20 %d/\n 21 %d/\n 24 %d/\n 0/\n" % (ma, mb, mc))
    case(corpus, "sel_bb", two, "moder\n 20 -30/\nmoder\n 21 -31/\n",
         "moder\n 1 -22/\n 'from binary tapes'/\n -30 %d/\n -31 %d/\n 0/\n" % (ma, mb))
    # a tape with two materials, made by the Fortran moder
    case(corpus, "sel_same", two,
         "moder\n 1 23/\n 'two materials'/\n 20 %d/\n 21 %d/\n 0/\n" % (ma, mb),
         "moder\n 1 22/\n 'both from one tape'/\n 23 %d/\n 23 %d/\n 0/\n" % (ma, mb))
    case(corpus, "sel_order", two, "",
         "moder\n 1 22/\n 'wrong order'/\n 21 %d/\n 20 %d/\n 0/\n" % (mb, ma))
    case(corpus, "sel_blank", two, "",
         "moder\n 1 22/\n ' '/\n 20 %d/\n 0/\n" % ma)
    n += 6
    print("variants:", n)


if __name__ == "__main__":
    main()
