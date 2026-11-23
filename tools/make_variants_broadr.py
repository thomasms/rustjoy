#!/usr/bin/env python3
"""Make extra BROADR cases with other input options, from cases made by
build_corpus.py --module broadr.

For each chosen base case, these variants are written (same tapes and
pre.input; only the broadr deck changes, except for the restart variant):
  t3   three final temperatures, no bootstrap
  bs   three final temperatures with bootstrap
  th   thinning only (temp2 = temp1 = 0)
  tn   negative thnmax (fixed upper limit) and explicit errmax/errint
  hi   high temperatures (> 9999 K, the other terminal format)
  rt   293.6 K (thermal quantities, integrals, g-factors)
  rs   restart: pre.input also runs the Fortran broadr to 300 K onto a
       new tape, and the run continues from it to 900 K with istart=1

With --multi BASE MAT2 (after the bases), two more cases are made from
BASE, whose ENDF tape must also hold MAT2 (after BASE's material):
  mm   reconr and broadr process both materials in one call each
  m2   reconr makes both materials, broadr processes MAT2 only (it skips
       the first material on the PENDF tape)

Usage: make_variants_broadr.py CORPUS_DIR base ... [--multi BASE MAT2]
"""
import os
import re
import shutil
import sys


def ints_and_reals(lines):
    """Tokens of the list-directed cards, one list per card ('/' ends)."""
    cards = []
    for line in lines:
        body = line.split("/")[0]
        toks = [t for t in re.split(r"[\s,]+", body.strip()) if t]
        if toks:
            cards.append(toks)
    return cards


def parse(run):
    lines = open(run).read().splitlines()
    assert lines[0].strip().lower().startswith("broadr")
    body = lines[1:]
    # drop the trailing 'stop'
    body = [l for l in body if l.strip().lower() != "stop"]
    cards = ints_and_reals(body)
    nendf, nin, nout = (int(x) for x in cards[0][:3])
    mat = int(cards[1][0])
    return nendf, nin, nout, mat


def deck(units, mat, card2, card3, temps):
    nendf, nin, nout = units
    return ("broadr\n %d %d %d/\n %d %s/\n %s/\n %s/\n 0/\nstop\n"
            % (nendf, nin, nout, mat, card2, card3, " ".join(temps)))


def multi(corpus, b, mat2):
    src = os.path.join(corpus, b)
    nendf, nin, nout, mat = parse(os.path.join(src, "run.input"))
    u = (nendf, nin, nout)
    pre = open(os.path.join(src, "pre.input")).read()
    # add MAT2 to the reconr deck: insert its cards before the final '0/'
    i = pre.lower().rindex("reconr")
    tail = pre[i:]
    m = re.search(r"\n\s*0\s*/[^\n]*\n", tail)
    assert m, "no reconr terminator"
    extra = "\n %d 1/\n .005/\n 'second material'" % mat2 + "/"
    pre2 = pre[:i] + tail[:m.start()] + extra + tail[m.start():]
    cases = {
        "mm": "broadr\n %d %d %d/\n %d 1 0 0 0./\n .002/\n 300./\n %d/\n 0/\nstop\n" % (u + (mat, mat2)),
        "m2": deck(u, mat2, "1 0 0 0.", ".002", ["300."]),
    }
    for v, d in cases.items():
        dst = os.path.join(corpus, "%s_%s" % (b, v))
        if os.path.exists(dst):
            sys.exit("%s exists (nothing is deleted)" % dst)
        shutil.copytree(src, dst, symlinks=True)
        with open(os.path.join(dst, "pre.input"), "w") as f:
            f.write(pre2)
        with open(os.path.join(dst, "run.input"), "w") as f:
            f.write(d)
    return len(cases)


def main():
    corpus = sys.argv[1]
    args = sys.argv[2:]
    extra = []
    if "--multi" in args:
        i = args.index("--multi")
        extra = args[i + 1:i + 3]
        args = args[:i]
    n = 0
    for b in args:
        src = os.path.join(corpus, b)
        nendf, nin, nout, mat = parse(os.path.join(src, "run.input"))
        u = (nendf, nin, nout)
        decks = {
            "t3": deck(u, mat, "3 0 0 0.", ".002", ["300.", "600.", "1200."]),
            "bs": deck(u, mat, "3 0 1 0.", ".002", ["300.", "900.", "2400."]),
            "th": deck(u, mat, "1 0 0 0.", ".01", ["0."]),
            "tn": deck(u, mat, "1 0 0 0.", ".002 -2.e4 .02 1.e-6", ["400."]),
            "hi": deck(u, mat, "2 0 0 0.", ".005", ["15000.", "40000."]),
            "rt": deck(u, mat, "1 0 0 0.", ".001", ["293.6"]),
        }
        for v, d in decks.items():
            dst = os.path.join(corpus, "%s_%s" % (b, v))
            if os.path.exists(dst):
                sys.exit("%s exists (nothing is deleted)" % dst)
            shutil.copytree(src, dst, symlinks=True)
            with open(os.path.join(dst, "run.input"), "w") as f:
                f.write(d)
            n += 1
        # restart: the Fortran broadens to 300 K onto unit 80 (same mode as
        # nin), then the run restarts from it
        dst = os.path.join(corpus, "%s_rs" % b)
        if os.path.exists(dst):
            sys.exit("%s exists (nothing is deleted)" % dst)
        shutil.copytree(src, dst, symlinks=True)
        t80 = -80 if nin < 0 else 80
        t81 = -81 if nin < 0 else 81
        pre = open(os.path.join(src, "pre.input")).read()
        assert pre.rstrip().endswith("stop")
        pre = pre.rstrip()[:-4] + deck((nendf, nin, t80), mat, "1 0 0 0.", ".002", ["300."])
        with open(os.path.join(dst, "pre.input"), "w") as f:
            f.write(pre)
        with open(os.path.join(dst, "run.input"), "w") as f:
            f.write(deck((nendf, t80, t81), mat, "1 1 0 300.", ".002", ["900."]))
        n += 1
    if extra:
        n += multi(corpus, extra[0], int(extra[1]))
    print("variants:", n)


if __name__ == "__main__":
    main()
