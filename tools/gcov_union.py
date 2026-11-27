#!/usr/bin/env python3
"""Union of the line coverage of one Fortran source from two gcov reports:
one of the unchanged source, one of the source with a patch applied.

A line counts as run when either report ran it. Lines of the patched report
are mapped back to the unchanged source through the hunks of the patch;
lines that the patch adds are left out, and the executable lines are those
of the unchanged report.

Usage: gcov_union.py SOURCE.gcov PATCHED.gcov PATCH_FILE
"""
import re
import sys


def load(p):
    d = {}
    for line in open(p, errors="replace"):
        parts = line.split(":", 2)
        if len(parts) < 3:
            continue
        c = parts[0].strip()
        try:
            n = int(parts[1])
        except ValueError:
            continue
        if n == 0 or c == "-":
            continue
        d[n] = c != "#####"
    return d


def mapping(patch):
    """new line number -> old line number (None for added lines), as a
    function."""
    hunks = []
    lines = open(patch).read().splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r"@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@", lines[i])
        if not m:
            i += 1
            continue
        o, nn = int(m.group(1)), int(m.group(3))
        hm = {}
        i += 1
        while i < len(lines) and not lines[i].startswith("@@"):
            t = lines[i][:1]
            if t == " ":
                hm[nn] = o
                o += 1
                nn += 1
            elif t == "-":
                o += 1
            elif t == "+":
                hm[nn] = None
                nn += 1
            i += 1
        hunks.append((min(hm) if hm else nn, max(hm) if hm else nn - 1, hm, nn - o))

    def f(n):
        delta = 0
        for lo, hi, hm, d in hunks:
            if n < lo:
                break
            if n <= hi:
                return hm.get(n)
            delta = d
        return n - delta
    return f


def main():
    src, pat, patch = sys.argv[1:4]
    a = load(src)
    b = load(pat)
    f = mapping(patch)
    u = dict(a)
    for n, c in b.items():
        o = f(n)
        if o is not None and o in u:
            u[o] = u[o] or c
    t = len(u)
    c = sum(1 for v in u.values() if v)
    print("%s + %s: %d of %d executable lines run (%.1f%%)" % (src, pat, c, t, 100.0 * c / t))


if __name__ == "__main__":
    main()
