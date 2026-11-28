#!/usr/bin/env python3
"""Inputs for tools/prepro/numcheck.f and `prustpro _numcheck`: random and
edge-case doubles (as hex bits) and random 11-character fields."""
import random, struct, sys

def hx(x):
    return struct.pack('>d', x).hex().upper()

r = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
out = sys.stdout
edges = [0.0, -0.0, 1.0, -1.0, 10.0, 0.1, 1e-3, 1e9, 1e-9, 1e10, 1e-99, 1e99,
         9.99999999995, 0.999999999, 999999999.9, 1.0e-3 * (1 - 1e-16),
         2.0e7, 1.0e-5, 5.0e-5, 1.234567890123e5]
for e in edges:
    out.write('D %s\n' % hx(e))
    out.write('D %s\n' % hx(-e))
for _ in range(n):
    k = r.random()
    if k < 0.3:
        x = 10.0 ** r.uniform(-102, 102)
    elif k < 0.6:
        x = 10.0 ** r.uniform(-12, 12)
    elif k < 0.8:
        # near decimal boundaries and rounding points
        e = r.randint(-12, 12)
        d = r.randint(1, 10)
        x = (r.randint(10 ** (d - 1), 10 ** d) + r.choice([0, 0.5, 0.4999999, 0.5000001, -0.5])) * 10.0 ** (e - d)
    else:
        x = struct.unpack('>d', struct.pack('>Q', r.getrandbits(64)))[0]
        if x != x or abs(x) == float('inf'):
            x = 1.5
    if r.random() < 0.5:
        x = -x
    out.write('D %s\n' % hx(x))
chars = '0123456789.+-eEdD  '
for _ in range(n):
    k = r.random()
    if k < 0.5:
        s = ''.join(r.choice(chars) for _ in range(11))
    else:
        m = r.randint(1, 9)
        s = ('%' + '.%dE' % m) % (r.uniform(-1, 1) * 10.0 ** r.uniform(-15, 15))
        s = s.replace('E', r.choice(['E', '', 'e', 'D'])).replace('+0', '+').replace('-0', '-')
        s = s[:11].rjust(11)
    out.write('S %s\n' % s)
