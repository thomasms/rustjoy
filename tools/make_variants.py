#!/usr/bin/env python3
"""Make extra RECONR cases with other input options, from the ASCII corpus
cases (which read the ENDF tape directly).

For each chosen base case and material, three variants are written:
  v1: err=0.01 (looser tolerance)
  v2: err=0.005, tempr=293.6 K (psi-chi broadening code paths, or the
      'not coded for temperature' error exits for MLBW/RM/hybrid)
  v3: err=0.001 with explicit errmax/errint, 2 descriptive cards and 5 user
      grid points

Usage: make_variants.py CORPUS_DIR [base ...]
"""
import os
import re
import shutil
import sys


def parse(run):
    """Return (nin, mat) of the first material of a reconr deck."""
    lines = open(run).read().splitlines()
    # line 0 'reconr', line 1 units, line 2 label, line 3 'mat ...'
    body = [l for l in lines if l.strip()]
    units = [int(x) for x in re.split(r"[\s,/]+", body[1].strip()) if re.fullmatch(r"[+-]?\d+", x)]
    i = 3
    mat = None
    while i < len(body):
        toks = [t for t in re.split(r"[\s,/]+", body[i].strip()) if t]
        if toks and re.fullmatch(r"[+-]?\d+", toks[0]):
            mat = int(toks[0])
            break
        i += 1
    return units[0], mat


VARIANTS = {
    "v1": (" {mat} 0 0/\n 0.01/\n", ),
    "v2": (" {mat} 0 0/\n 0.005 293.6/\n", ),
    "v3": (" {mat} 2 5/\n 0.001 0. 0.004 1.e-6/\n 'variant case card one'/\n 'variant case card two'/\n"
           " 0.5 2.53e-2 7.77 1234.5678 2.e6/\n", ),
}


def main():
    corpus = sys.argv[1]
    bases = sys.argv[2:]
    n = 0
    for b in bases:
        src = os.path.join(corpus, b)
        nin, mat = parse(os.path.join(src, "run.input"))
        for v, (cards,) in VARIANTS.items():
            dst = os.path.join(corpus, "%s_%s" % (b, v))
            if os.path.exists(dst):
                shutil.rmtree(dst)
            shutil.copytree(src, dst, symlinks=True)
            with open(os.path.join(dst, "run.input"), "w") as f:
                f.write("reconr\n %d 98/\n 'variant of %s'/\n" % (nin, b))
                f.write(cards.format(mat=mat))
                f.write(" 0/\nstop\n")
            n += 1
    print("variants:", n)


if __name__ == "__main__":
    main()
