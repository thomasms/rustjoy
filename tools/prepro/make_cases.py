#!/usr/bin/env python3
"""Make test cases for prustpro: one directory per case with the ENDF input
and the program's .INP file.

usage: make_cases.py PROGRAM ENDF_DIR OUT_DIR

ENDF_DIR holds ENDF files (for example NJOY2016 tests/resources). Each
program has a list of input decks (variants of the PREPRO test inputs);
every deck is combined with every ENDF file it can use. A case directory
is OUT_DIR/PROGRAM/<deck>-<file>.
"""
import os
import shutil
import sys

def endf_files(d):
    out = []
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        if not os.path.isfile(p) or f.endswith('.ace'):
            continue
        with open(p, 'rb') as h:
            head = h.read(400)
        # an ENDF tape: second line ends with MF=1 MT=451
        lines = head.split(b'\n')
        if len(lines) > 2 and lines[1][70:75].strip() == b'1451':
            out.append(f)
    return out

def mat_range(s):
    return s

LINEAR = {
    # the PREPRO test input: its first line does not fit the format
    # (2I11,11A1,2I11), so LINEAR takes the ERR= branch and the defaults
    'std': """          0 1.0000E-30          1          1          1          1
ENDFB.IN
ENDFB.OUT
          1       9999
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.0000E-00 1.0000E-03
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
""",
    'law': """          0          1 1.0000E-30          0          0
ENDFB.IN
ENDFB.OUT
     0 0  0  999999999
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.0000E-00 1.0000E-03
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
 0.0000E-00 1.0000E-05
 1.0000E+00 1.0000E-05
 2.0000E+00 1.0000E-04
 2.0000E+07 1.0000E-04
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
""",
    'law3': """          0          0 1.0000E-30          1          1
ENDFB.IN
ENDFB.OUT
     0 0  0  999999999
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 1.0000E-05 1.0000E-02
 1.0000E+02 5.0000E-03
 1.0000E+05
 2.0000E+07 2.0000E-03

""",
    'za': """          1          0 1.0000E-30          0          0
ENDFB.IN
ENDFB.OUT
      1 0  0999999 0  0
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.0000E-00 1.0000E-04
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
""",
    'sel': """          0          1 1.0000E-30          1          0
ENDFB.IN
ENDFB.OUT
     1 3  1  9999 3  4
     1 3102  9999 3999
     1 9  0  999910  0
     123 0  9999270
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.0000E-00 3.0000E-03
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
""",
    'none': None,  # no LINEAR.INP: all defaults
    'eof': """          0          0 1.0000E-30          0          1
ENDFB.IN
ENDFB.OUT
     0 0  0  999999999
 0.0000E-00 1.0000E-03
 1.0000E+01 2.0000E-03""",
}

DICTIN = {
    'std': 'ENDFB.IN\nENDFB.OUT\n',
    'none': None,
    'blank': '\n\n',
}

def legend_deck(err, npt, tabs, legs, fix, var, table=''):
    return '%11s%11d%11d%11d%11d%11d\nENDFB.IN\nENDFB.OUT\n%s\n' % (err, npt, tabs, legs, fix, var, table)

LEGEND = {
    'std': legend_deck(' 1.0000E-02', 20000, 2, 1, 2, 0, '     0 0  0  999999999'),
    'def': None,
    'coef': legend_deck(' 1.0000E-03', 0, 1, 2, 1, 1),
    'copy': legend_deck(' 1.0000E-03', 5, 0, 0, 0, 0),
    'fix1': legend_deck(' 1.0000E-04', 3000, 2, 1, 1, 1,
                        '     1 4  2  9999 4  2 0.0        2.0+7      1.0-4      0.02\n     1 4  0  9999 4999 1.0+5      0.0        0.002      0.05'),
    'fix2': legend_deck(' 5.0000E-03', 100, 1, 1, 2, 0, '  2000 4  2  3000 4  2'),
}

DECKS = {'legend': ('LEGEND.INP', LEGEND), 'linear': ('LINEAR.INP', LINEAR), 'dictin': ('DICTIN.INP', DICTIN),
         'activate': ('ACTIVATE.INP', DICTIN)}

def link(case, name, src):
    dst = os.path.join(case, name)
    if not os.path.lexists(dst):
        os.symlink(os.path.abspath(src), dst)

def write(case, name, text):
    with open(os.path.join(case, name), 'w') as h:
        h.write(text)

def merger_cases(endf_dir, out):
    """MERGER: several tapes, request lists by MAT or ZA."""
    files = endf_files(endf_dir)
    groups = [files[i:i + 3] for i in range(0, len(files), 3)]
    reqs = {
        'all': ('    0', ''),
        'mat': ('    0', '   100 0  0  2000 0  0          7\n  2600 3  1  9999 3999\n  9200 1  0  9499 3 20         33\n'),
        'za': ('    1', '  1000 0  0 30000 0  0\n 90000 3  1 99000 3102\n'),
        'mf3': ('    0', '     0 3  0     0 3999\n'),
    }
    n = 0
    for gi, g in enumerate(groups):
        for rname, (how, rq) in reqs.items():
            case = os.path.join(out, 'merger', '%s-%02d' % (rname, gi))
            os.makedirs(case, exist_ok=True)
            names = []
            for k, f in enumerate(g):
                nm = 'TAPE%d.ENDF' % k
                link(case, nm, os.path.join(endf_dir, f))
                names.append(nm)
            # the same tape twice (duplicates), and a label line
            label = 'prustpro merger test %-46s' % rname
            text = 'MERGED.OUT\n%s%s\n%s\n%s\nEND\n%s\n' % (
                label[:66], '6000' + how[-2:], '\n'.join(names), names[0], rq)
            write(case, 'MERGER.INP', text)
            n += 1
    # special cases
    sp = {
        'blankname': 'X.OUT\n%-66s   0 0\n\n\nEND\n\n' % 'label',
        'missing': '\n%-66s   0 0\nNOFILE\nTAPE0.ENDF\nEND\n\n' % 'label',
        'one': '\n%-66s   0 0\nTAPE0.ENDF\nend\n     0 3  0     0 3999\n\n' % '',
        'none': '\n%-66s   0 0\nEND\n' % 'label',
        'onlyone': 'OUT\n%-66s   0 0\nTAPE0.ENDF\nEND\n\n' % 'label',
        'unsorted': 'OUT\n%-66s   0 0\nREV.ENDF\nTAPE0.ENDF\nEND\n\n' % 'label',
        'empty': 'OUT\n%-66s   0 0\nEMPTY\nTAPE0.ENDF\nEND\n\n' % 'label',
        'nomatch': 'OUT\n%-66s   0 0\nTAPE0.ENDF\nTAPE1.ENDF\nEND\n  9990 3  1  9991 3  1\n\n' % 'label',
    }
    for name, text in sp.items():
        case = os.path.join(out, 'merger', 'sp-' + name)
        os.makedirs(case, exist_ok=True)
        link(case, 'TAPE0.ENDF', os.path.join(endf_dir, files[0]))
        link(case, 'TAPE1.ENDF', os.path.join(endf_dir, files[5]))
        link(case, 'ENDFB.IN1', os.path.join(endf_dir, files[1]))
        link(case, 'ENDFB.IN2', os.path.join(endf_dir, files[2]))
        write(case, 'EMPTY', '')
        # sections in reverse MT order
        lines = open(os.path.join(endf_dir, files[3]), 'rb').read().split(b'\n')
        write(case, 'MERGER.INP', text)
        body = [l for l in lines[1:] if len(l) >= 75]
        secs, cur = [], []
        for l in body:
            cur.append(l)
            if l[72:75].strip() in (b'0', b''):
                secs.append(cur)
                cur = []
        with open(os.path.join(case, 'REV.ENDF'), 'wb') as h:
            h.write(lines[0] + b'\n')
            for s in reversed(secs[:-3]):
                h.write(b'\n'.join(s) + b'\n')
            for s in secs[-3:]:
                h.write(b'\n'.join(s) + b'\n')
        n += 1
    print(n, 'cases')

def relabel_cases(src_dir, out):
    """RELABEL: FORTRAN source files (for example the PREPRO sources)."""
    n = 0
    for root, _, fs in os.walk(src_dir):
        for f in sorted(fs):
            if not (f.endswith('.f') or f.endswith('.h')):
                continue
            case = os.path.join(out, 'relabel', os.path.relpath(os.path.join(root, f), src_dir).replace('/', '_'))
            os.makedirs(case, exist_ok=True)
            link(case, 'SRC.F', os.path.join(root, f))
            write(case, 'RELABEL.INP', 'SRC.F\nSRC.OUT\n')
            n += 1
    for name, text in {'noinp': None, 'blank': '\n\n', 'missing': 'NOFILE\nX\n', 'one': 'SRC.F\n'}.items():
        case = os.path.join(out, 'relabel', 'sp-' + name)
        os.makedirs(case, exist_ok=True)
        link(case, 'RELABEL.IN', os.path.join(src_dir, 'merger.f'))
        link(case, 'SRC.F', os.path.join(src_dir, 'dictin.f'))
        if text is not None:
            write(case, 'RELABEL.INP', text)
        n += 1
    print(n, 'cases')

SYNTH = """      PROGRAM TEST
     1 (A,B)
C---- comment that is copied
C***** DOUBLE
      IMPLICIT REAL*8 (A-H,O-Z)
C***** DOUBLE
C***** SINGLE
C     IMPLICIT REAL*4 (A-H,O-Z)
C***** SINGLE
C*****CHARACTER
      CHARACTER*4 A
C*****CHARACTER
C***INTEGER
CXX   INTEGER A
C***INTEGER
C*** STOP
      STOP
C*** STOP
C***** EXIT
C     CALL EXIT
C***** EXIT
C***** CIVIC
C***** NAMES
      OPEN(10,FILE='X')
C***** NAMES
C***** CIVIC
C***** REWIND
C***** SAVE
      SAVE
C***** SAVE
C***** CLOSE
      CLOSE(10)
C***** CLOSE
C***** ACTION
C***** ACTION
C***** ENDFILE
C     ENDFILE 10
C***** ENDFILE
C***** REWIND
C***** UNKNOWN
      X = 1.0          longer line with trailing text in columns 73-80 ABCDEFGH
C*********************
      END
"""

def convert_cases(src_dir, out):
    """CONVERT: keyword lists on synthetic and real sources."""
    keys = ['DOUBLE INTEGER PROGRAM', 'SINGLE CHARACTER STOP NOID', 'DOUBLE SINGLE INTEGER CHARACTER',
            'EXIT CIVIC NAMES SAVE ACTION CLOSE REWIND ENDFILE', '', 'NOID PROGRAM',
            'DOUBLE BADWORD', ' DOUBLE  STOP   EXIT     INTEGER    NAMES  SAVE   ACTION   CLOSE  CIVIC REWIND',
            'PROGRAMX DOUBLE', 'DOUBLE' + ' ' * 60 + 'INTEGER']
    srcs = {'synth': None}
    for f in sorted(os.listdir(src_dir)):
        if f in ('recent.f', 'sixpak.f', 'sigma1.f', 'scratcha.f', 'merger.f'):
            srcs[f] = os.path.join(src_dir, f)
    n = 0
    for k, key in enumerate(keys):
        for s, path in srcs.items():
            case = os.path.join(out, 'convert', 'k%02d-%s' % (k, s))
            os.makedirs(case, exist_ok=True)
            if path is None:
                write(case, 'SRC.F', SYNTH)
            else:
                link(case, 'SRC.F', path)
            write(case, 'CONVERT.INP', key + '\nSRC.F\nSRC.OUT\n')
            n += 1
    for name, text in {'noinp': None, 'nonames': 'DOUBLE\n', 'missing': 'DOUBLE\nNOFILE\n\n', 'blank': 'DOUBLE\n\n\n'}.items():
        case = os.path.join(out, 'convert', 'sp-' + name)
        os.makedirs(case, exist_ok=True)
        write(case, 'CONVERT.IN', SYNTH)
        if text is not None:
            write(case, 'CONVERT.INP', text)
        n += 1
    print(n, 'cases')

def combine(paths, dst):
    """One tape from several: label of the first, every material, TEND."""
    with open(dst, 'wb') as h:
        for k, p in enumerate(paths):
            lines = open(p, 'rb').read().split(b'\n')
            if k == 0:
                h.write(lines[0] + b'\n')
            for l in lines[1:]:
                if l[66:70].strip() == b'-1':
                    break
                if l:
                    h.write(l + b'\n')
        h.write(b' ' * 66 + b'  -1 0  0\n')

def mixer_cases(lin_dir, out):
    """MIXER: mixtures from linearized tapes (made by LINEAR)."""
    L = lambda n: os.path.join(lin_dir, n + '.lin')
    tapes = {
        'hb': [L('n-001_H_002-ENDF8.0.endf'), L('n-005_B_010-ENDF8.0.endf')],
        'fesral': [L('n-026_Fe_056-JEFF3.3.endf'), L('n-038_Sr_088-ENDF8.1.endf'), L('n-013_Al_027-ENDF8.0.endf')],
        'uzr': [L('n-040_Zr_090-ENDF8.0.endf'), L('n-092_U_235-ENDF8.0.endf')],
        'gam': [L('photoat-001_H_000-ENDF8.0.endf'), L('photoat-092_U_000-ENDF8.0.endf')],
    }
    def comp(rows):
        return ''.join('%11d       %4d%11s\n' % r for r in rows)
    decks = {
        'hb-g': ('hb', '      10010  1001 3  1          0\n', [(1002, 1, '0.1'), (5010, 1, '0.9')]),
        'hb-a': ('hb', '      10010  1001 3  2          1\n', [(1002, 2, '2.0'), (5010, 2, '1.0')]),
        'hb-102': ('hb', '      10010  1001 3102          0\n', [(1002, 102, '1.0'), (5010, 102, '1.0'), (1002, 102, '0.5')]),
        'fsa-1': ('fesral', '      99000  9900 3  1          0\n', [(26056, 1, '0.5'), (38088, 1, '.3'), (13027, 1, '2.E-1')]),
        'fsa-102': ('fesral', '          0     0 0  0          1\n', [(26056, 102, '1.0E+00'), (38088, 102, '1.0'), (13027, 102, '1.0')]),
        'fsa-miss': ('fesral', '      99000  9900 3  1          0\n', [(26056, 1, '0.5'), (40090, 1, '0.5')]),
        'fsa-mt': ('fesral', '      99000  9900 3  4          0\n', [(26056, 4, '0.5'), (38088, 999, '0.5')]),
        'fsa-bad': ('fesral', '      99000  9900 3  1          0\n', [(26056, 1, '0.5'), (38088, 1, '-0.5')]),
        'fsa-over': ('fesral', '      99000  9900 3  1          0\n', [(26056, 1, '0.1')] * 11),
        'fsa-ten': ('fesral', '      99000  9900 3  1          0\n', [(26056, 1, '0.1')] * 5 + [(13027, 1, '0.1')] * 5),
        'fsa-none': ('fesral', '      99000  9900 3  1          0\n', []),
        'uzr-1': ('uzr', '      92000  9200 3  1          0\n', [(92235, 1, '0.9'), (40090, 1, '0.1')]),
        'uzr-18': ('uzr', '      92000  9200 3 18          0\n', [(92235, 18, '1.0'), (40090, 18, '1.0')]),
        'gam-501': ('gam', '      92000  9200 23501          0\n', [(1000, 501, '0.1'), (92000, 501, '0.9')]),
        'gam-502': ('gam', '      92000  9200 23502          1\n', [(1000, 502, '1.0'), (92000, 502, '1.0')]),
    }
    n = 0
    for name, (tape, outline, rows) in decks.items():
        case = os.path.join(out, 'mixer', name)
        os.makedirs(case, exist_ok=True)
        combine(tapes[tape], os.path.join(case, 'MIX.IN'))
        text = ('prustpro MIXER test %s\nsecond title line\nMIX.IN\nMIX.OUT\n' % name) + outline + comp(rows) + '\n'
        write(case, 'MIXER.INP', text)
        n += 1
    for name, text in {'noinp': None, 'short': 'title\n', 'eof': 'a\nb\nMIX.IN\nMIX.OUT\n',
                       'missing': 'a\nb\nNOFILE\nX\n          1     1 3  1          0\n       1002          1        1.0\n\n'}.items():
        case = os.path.join(out, 'mixer', 'sp-' + name)
        os.makedirs(case, exist_ok=True)
        combine(tapes['hb'], os.path.join(case, 'MIX.IN'))
        if text is not None:
            write(case, 'MIXER.INP', text)
        n += 1
    print(n, 'cases')

def main():
    prog, endf_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]
    if prog == 'mixer':
        return mixer_cases(endf_dir, out)
    if prog == 'convert':
        return convert_cases(endf_dir, out)
    if prog == 'relabel':
        return relabel_cases(endf_dir, out)
    if prog == 'merger':
        return merger_cases(endf_dir, out)
    inp_name, decks = DECKS[prog]
    files = endf_files(endf_dir)
    n = 0
    for deck, text in decks.items():
        for f in files:
            case = os.path.join(out, prog, '%s-%s' % (deck, f))
            os.makedirs(case, exist_ok=True)
            dst = os.path.join(case, 'ENDFB.IN')
            if not os.path.exists(dst):
                os.symlink(os.path.abspath(os.path.join(endf_dir, f)), dst)
            if text is not None:
                with open(os.path.join(case, inp_name), 'w') as h:
                    h.write(text)
            n += 1
    print(n, 'cases')

main()
