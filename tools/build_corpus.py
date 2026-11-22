#!/usr/bin/env python3
"""Build the RECONR comparison corpus from the NJOY2016 test suite.

For every test deck that calls reconr, and for every reconr call in it, a case
directory is created with:
  tapes/      the input tapes the test copies in (from tests/resources)
  pre.input   the moder modules that run before this reconr call (or nothing)
  run.input   the reconr block alone, followed by 'stop'

A second, ASCII variant of each case (suffix _a) reads the original ASCII
ENDF tape directly and writes an ASCII PENDF tape, so that both the
formatted and the binary code paths are exercised.

Usage: build_corpus.py NJOY_SRC_DIR OUT_DIR
"""
import os
import re
import shutil
import sys

MODULES = {"reconr", "broadr", "unresr", "heatr", "thermr", "groupr", "gaminr",
           "errorr", "covr", "moder", "dtfr", "ccccr", "matxsr", "resxsr", "acer",
           "powr", "wimsr", "plotr", "viewr", "mixr", "purr", "leapr", "gaspr",
           "stop", "--"}


def module_of(line):
    t = line.strip()
    if not t:
        return None
    tok = re.split(r"[\s/,]", t, maxsplit=1)[0].lower()
    if tok.startswith("--"):
        return "--"
    return tok if tok in MODULES else None


def split_modules(text):
    """Return a list of (name, [lines]) for the deck."""
    mods = []
    for line in text.splitlines():
        m = module_of(line)
        if m is not None:
            mods.append([m, [line]])
        elif mods:
            mods[-1][1].append(line)
    return mods


def first_ints(lines, n):
    """The first n integers of a card sequence (list-directed, '/' ends)."""
    vals = []
    for line in lines:
        body = line.split("/")[0]
        for tok in re.split(r"[\s,]+", body.strip()):
            if tok:
                vals.append(int(tok))
                if len(vals) == n:
                    return vals
        if "/" in line:
            break
    return vals


def tapes_of(cmake):
    m = {}
    for res, tape in re.findall(r'configure_file\(\s*"\$\{RESOURCES\}/([^"]+)"\s*"\$\{CMAKE_CURRENT_BINARY_DIR\}/(tape\d+)"', cmake):
        m[tape] = res
    return m


def main():
    src, out = sys.argv[1], sys.argv[2]
    tests = os.path.join(src, "tests")
    os.makedirs(out, exist_ok=True)
    ncase = 0
    for t in sorted(os.listdir(tests)):
        d = os.path.join(tests, t)
        inp = os.path.join(d, "input")
        if not os.path.isfile(inp):
            continue
        text = open(inp, errors="replace").read()
        mods = split_modules(text)
        tapes = tapes_of(open(os.path.join(d, "CMakeLists.txt")).read())
        k = 0
        for idx, (name, lines) in enumerate(mods):
            if name != "reconr":
                continue
            pre = [m for m in mods[:idx] if m[0] in ("moder", "--")]
            # only keep cases whose earlier modules are all moder (or comments)
            if any(m[0] not in ("moder", "--") for m in mods[:idx]):
                # tapes from other modules (e.g. a previous reconr) are not
                # needed by reconr itself unless it reads them; skip if the
                # reconr input tape is produced by a non-moder module.
                nin = first_ints(lines[1:] if len(lines) > 1 else [], 1)
                card_lines = lines[1:]
                if lines[0].strip().lower() != "reconr":
                    card_lines = [lines[0].split(None, 1)[1] if len(lines[0].split(None, 1)) > 1 else ""] + lines[1:]
                nin = first_ints(card_lines, 1)
                producers = [m for m in mods[:idx] if m[0] not in ("moder", "--")]
                made = set()
                for m in producers:
                    cl = m[1][1:]
                    try:
                        u = first_ints(cl, 2)
                        made.update(abs(x) for x in u[1:])
                    except Exception:
                        pass
                if nin and abs(nin[0]) in made:
                    continue
            k += 1
            name_c = "%s_%d" % (t, k)
            case = os.path.join(out, name_c)
            if os.path.exists(case):
                sys.exit("%s exists: give a new corpus directory (nothing is deleted)" % case)
            os.makedirs(os.path.join(case, "tapes"))
            for tape, res in tapes.items():
                os.symlink(os.path.join(os.path.abspath(tests), "resources", res), os.path.join(case, "tapes", tape))
            with open(os.path.join(case, "pre.input"), "w") as f:
                for m in pre:
                    f.write("\n".join(m[1]) + "\n")
                f.write("stop\n")
            with open(os.path.join(case, "run.input"), "w") as f:
                f.write("\n".join(lines) + "\nstop\n")
            ncase += 1

            # ASCII variant: replace a moder-produced binary input by its
            # ASCII source, and write an ASCII PENDF tape on unit 98.
            card_lines = lines[1:]
            head = lines[0].strip()
            rest = head.split(None, 1)
            inline = rest[1] if len(rest) > 1 else ""
            all_cards = ([inline] if inline.strip() else []) + card_lines
            try:
                nin, npend = first_ints(all_cards, 2)
            except ValueError:
                continue
            src_tape = None
            if nin > 0:
                src_tape = nin
            else:
                for m in pre:
                    if m[0] != "moder":
                        continue
                    u = first_ints(m[1][1:] if m[1][0].strip().lower() == "moder" else m[1], 2)
                    if len(u) == 2 and u[1] == nin and u[0] > 0:
                        src_tape = u[0]
            if src_tape is None or ("tape%d" % src_tape) not in tapes:
                continue
            # rewrite the first card
            new_lines = []
            replaced = False
            body = "\n".join(lines)
            # find first two integers in the text after 'reconr' and replace
            mod_pos = body.lower().index("reconr") + len("reconr")
            tail = body[mod_pos:]
            mm = re.search(r"([+-]?\d+)([\s,]+)([+-]?\d+)", tail)
            if not mm:
                continue
            tail2 = tail[:mm.start()] + "%d%s%d" % (src_tape, mm.group(2), 98) + tail[mm.end():]
            casea = case + "_a"
            if os.path.exists(casea):
                sys.exit("%s exists: give a new corpus directory (nothing is deleted)" % casea)
            shutil.copytree(case, casea, symlinks=True)
            with open(os.path.join(casea, "pre.input"), "w") as f:
                f.write("stop\n")
            with open(os.path.join(casea, "run.input"), "w") as f:
                f.write(body[:mod_pos] + tail2 + "\nstop\n")
            ncase += 1
    print("cases:", ncase)


if __name__ == "__main__":
    main()
