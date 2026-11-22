#!/usr/bin/env python3
"""Write the input of the numeric-kernel check (tools/numcheck.f90 and
examples/numcheck.rs) to standard output: pairs of (real(8) x, int32 n),
little endian. Edge cases first, then 2,000,000 seeded random values."""
import random, struct, sys, math
random.seed(12345)
out = bytearray()
def put(x, n):
    out.extend(struct.pack('<di', x, n))
edges = [0.0, -0.0, 0.1, 0.0999999999, 0.1000000001, 1.0, 0.999999999, 0.9999999995, 9999999.99, 9999999.995,
         1e7, 9.9999995, 9.999995, 9.99995, 9.999999995, 9.99999995, 1e-5, 2e7, 1e-38, 1e38, 1e100, 1e-100,
         2.5e-5, 0.5, 1.5, 0.125, 1234567.5, 123456.75]
for x in edges:
    for s in (1, -1):
        put(s * x, 3)
# values near every power of ten and the a11 thresholds, nudged by ulps
for e in range(-40, 41):
    for base in (1.0, 9.9999995, 9.999995, 9.99995, 9.999999995, 9.99999995, 5.0, 2.0):
        v = base * 10.0 ** e
        for k in range(-3, 4):
            w = v
            for _ in range(abs(k)):
                w = math.nextafter(w, math.inf if k > 0 else -math.inf)
            put(w, random.randint(-20, 20))
            put(-w, random.randint(-20, 20))
# random: uniform in log magnitude, typical ENDF ranges and beyond
for i in range(2_000_000):
    r = random.random()
    if r < 0.6:
        x = 10.0 ** random.uniform(-12, 12)
    elif r < 0.8:
        x = 10.0 ** random.uniform(-100, 100)
    elif r < 0.9:
        # exact decimals with few digits (ties in formatting)
        x = random.randint(1, 10**8) / 10 ** random.randint(0, 9)
    else:
        x = random.uniform(0, 1e7)
    if random.random() < 0.3:
        x = -x
    put(x, random.randint(-30, 30))
sys.stdout.buffer.write(out)
