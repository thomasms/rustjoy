#!/usr/bin/env python3
"""Run every corpus case with the Fortran NJOY2016 and the Rust port, and
compare the results byte for byte.

What is compared, per case:
  * every tape file that the reconr run created or changed (byte for byte)
  * the listing file 'output'
  * the terminal output (stdout) and the exit status

The only bytes that are masked before comparison are:
  * run-time values that differ between any two runs (CPU seconds at the end
    of listing lines, the date and time in the banner)
  * bytes 68..135 of the tape-label record of a binary (unformatted) PENDF
    tape. The Fortran fills these from uninitialised stack memory
    (equivalence of character(4)::tlabel(17) with real(kr)::rlabel(17)), so
    they differ even between two runs of the Fortran program itself.

Usage: compare.py FORTRAN_NJOY RUST_BIN CORPUS_DIR WORK_DIR [case ...]
"""
import filecmp
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import time

TIME_RE = re.compile(rb"[ 0-9.-]{8}s$", re.M)
TIME_FIELD = re.compile(rb"^ *-?\d+\.\d$")
DATE_RE = re.compile(rb"(date: )\d\d/\d\d/\d\d")
CLOCK_RE = re.compile(rb"(time: )\d\d:\d\d:\d\d")
TERM_DT_RE = re.compile(rb"\d\d/\d\d/\d\d \d\d:\d\d:\d\d")


def mask_text(b):
    b = TIME_RE.sub(lambda m: b"#time#s" if TIME_FIELD.match(m.group(0)[:8]) else m.group(0), b)
    b = DATE_RE.sub(rb"\1##/##/##", b)
    b = CLOCK_RE.sub(rb"\1##:##:##", b)
    b = TERM_DT_RE.sub(b"##/##/## ##:##:##", b)
    return b


def is_unformatted(b):
    if len(b) < 8:
        return False
    n = struct.unpack("<i", b[:4])[0]
    if n <= 0 or n + 8 > len(b):
        return False
    return b[4 + n:8 + n] == b[:4]


def mask_tape(b):
    """Mask the bytes of a binary tape that the Fortran leaves undefined.

    1. bytes 68..135 of the tape-label record (rlabel/tlabel equivalence in
       reconr);
    2. bytes 68..135 of every 136-byte Hollerith card block in MF1/MT451
       records written by hdatio (z/tz equivalence in recout).
    """
    if not is_unformatted(b):
        return b
    b = bytearray(b)
    off = 0
    first = True
    in_hol = False  # inside an MF1/MT451 Hollerith record sequence
    while off + 4 <= len(b):
        n = struct.unpack("<i", b[off:off + 4])[0]
        d0 = off + 4
        mat, mf, mt, nb, nw = struct.unpack("<5i", b[d0:d0 + 20])
        data = d0 + 20
        if first:
            if mf == 0 and mt == 0 and nw == 17 and n == 20 + 136:
                for i in range(data + 68, data + 136):
                    b[i] = 0
            first = False
        elif mf == 1 and mt == 451:
            blocks_at = None
            if not in_hol and nw > 6 and (nw - 6) % 17 == 0 and n == 20 + 8 * nw:
                # hdatio first record: 6 header words then cards
                # (only the record that follows the CONT records with
                # c1/c2 = temperature/err has cards; detect by size)
                blocks_at = data + 48
                ncard = (nw - 6) // 17
                in_hol = nb > 0
            elif in_hol and nw % 17 == 0:
                blocks_at = data
                ncard = nw // 17
                in_hol = nb > 0
            if blocks_at is not None:
                for k in range(ncard):
                    s0 = blocks_at + 136 * k
                    for i in range(s0 + 68, s0 + 136):
                        b[i] = 0
        else:
            in_hol = False
        off += 8 + n
    return bytes(b)


def snapshot(d):
    s = {}
    for f in os.listdir(d):
        p = os.path.join(d, f)
        if os.path.isfile(p) or os.path.islink(p):
            st = os.lstat(p)
            s[f] = (st.st_size, st.st_mtime_ns)
    return s


def run(binary, d, inp, timeout):
    t0 = time.time()
    with open(os.path.join(d, inp), "rb") as fi:
        p = subprocess.run([binary], stdin=fi, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, cwd=d, timeout=timeout)
    return p.returncode, p.stdout, p.stderr, time.time() - t0


def first_diff(a, b):
    n = min(len(a), len(b))
    for i in range(n):
        if a[i] != b[i]:
            return i
    return n if len(a) != len(b) else -1


def line_context(a, b, off):
    la = a[:off].count(b"\n") + 1
    al = a.splitlines()
    bl = b.splitlines()
    i = la - 1
    return {
        "line": la,
        "fortran": al[i].decode("latin1") if i < len(al) else None,
        "rust": bl[i].decode("latin1") if i < len(bl) else None,
    }


def compare_case(fortran, rust, case_dir, work, timeout):
    name = os.path.basename(case_dir)
    base = os.path.join(work, name)
    if os.path.exists(base):
        sys.exit("%s exists: give a new work directory (nothing is deleted)" % base)
    pre = os.path.join(base, "pre")
    shutil.copytree(os.path.join(case_dir, "tapes"), pre)
    shutil.copy(os.path.join(case_dir, "pre.input"), pre)
    shutil.copy(os.path.join(case_dir, "run.input"), pre)
    rc, _, _, _ = run(fortran, pre, "pre.input", timeout)
    if rc != 0:
        return {"case": name, "status": "setup-failed"}
    # keep the moder listing, but not under the name the RECONR run writes
    if os.path.exists(os.path.join(pre, "output")):
        os.rename(os.path.join(pre, "output"), os.path.join(pre, "pre_output.txt"))
    before = snapshot(pre)
    res = {"case": name}
    dirs = {}
    for who, binary in (("fortran", fortran), ("rust", rust)):
        d = os.path.join(base, who)
        os.makedirs(d)
        for f in os.listdir(pre):
            os.symlink(os.path.join(pre, f), os.path.join(d, f))
        snap0 = snapshot(d)
        try:
            rc, out, err, dt = run(binary, d, "run.input", timeout)
        except subprocess.TimeoutExpired:
            rc, out, err, dt = "timeout", b"", b"", timeout
        res[who] = {"rc": rc, "seconds": round(dt, 2), "stderr": err.decode("latin1")[-400:]}
        dirs[who] = (d, out, snap0)
    fd, fout, f0 = dirs["fortran"]
    rd, rout, r0 = dirs["rust"]
    files = {}
    changed_f = {f for f, v in snapshot(fd).items() if f0.get(f) != v}
    changed_r = {f for f, v in snapshot(rd).items() if r0.get(f) != v}
    allf = sorted(changed_f | changed_r)
    ok = res["fortran"]["rc"] == res["rust"]["rc"]
    for f in allf:
        pf, pr = os.path.join(fd, f), os.path.join(rd, f)
        if not os.path.exists(pf) or not os.path.exists(pr):
            files[f] = {"status": "missing-in-" + ("rust" if os.path.exists(pf) else "fortran")}
            ok = False
            continue
        a, b = open(pf, "rb").read(), open(pr, "rb").read()
        if f == "output":
            a, b = mask_text(a), mask_text(b)
        else:
            a, b = mask_tape(a), mask_tape(b)
        if a == b:
            files[f] = {"status": "identical", "bytes": len(a)}
        else:
            ok = False
            off = first_diff(a, b)
            info = {"status": "differs", "bytes_fortran": len(a), "bytes_rust": len(b), "first_diff_offset": off}
            if not is_unformatted(a):
                info.update(line_context(a, b, off))
            files[f] = info
    a, b = mask_text(fout), mask_text(rout)
    files["<stdout>"] = {"status": "identical" if a == b else "differs"}
    if a != b:
        ok = False
        off = first_diff(a, b)
        files["<stdout>"].update(line_context(a, b, off))
    res["files"] = files
    res["status"] = "identical" if ok else "differs"
    return res


def main():
    fortran, rust, corpus, work = [os.path.abspath(x) for x in sys.argv[1:5]]
    only = sys.argv[5:]
    timeout = int(os.environ.get("CASE_TIMEOUT", "1800"))
    os.makedirs(work, exist_ok=True)
    results = []
    for c in sorted(os.listdir(corpus)):
        if only and c not in only:
            continue
        r = compare_case(fortran, rust, os.path.join(corpus, c), work, timeout)
        results.append(r)
        print("%-10s %-10s F:%-4s R:%-4s %s" % (
            r["case"], r["status"],
            r.get("fortran", {}).get("rc"), r.get("rust", {}).get("rc"),
            " ".join("%s=%s" % (k, v["status"]) for k, v in r.get("files", {}).items() if v["status"] != "identical")),
            flush=True)
    json.dump(results, open(os.path.join(work, "results.json"), "w"), indent=1)
    n = len(results)
    ni = sum(1 for r in results if r["status"] == "identical")
    print("identical: %d / %d" % (ni, n))


if __name__ == "__main__":
    main()
