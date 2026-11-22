#!/usr/bin/env python3
import gzip
import re
import sys
from itertools import zip_longest

def reads(path, mate):
    with gzip.open(path, 'rb') as f:
        while True:
            h = f.readline()
            if not h:
                return
            s, p, q = f.readline(), f.readline(), f.readline()
            m = re.fullmatch(rb'@([^\s]+)#([12]):([YN]):([0-9]+):([^\s]+)\s*', h)
            if not m or m[2] != str(mate).encode() or not p.startswith(b'+') or not s.strip() or len(s.rstrip()) != len(q.rstrip()):
                raise ValueError('Invalid record in ' + path)
            # Move CASAVA metadata into its standard second whitespace-delimited field.
            header = b'@' + m[1] + b' ' + b':'.join(m.group(i) for i in (2,3,4,5)) + b'\n'
            yield m[1], header + s + b'+\n' + q

r1, r2, out1, out2 = sys.argv[1:]
n = 0
with gzip.open(out1, 'wb', compresslevel=1) as a, gzip.open(out2, 'wb', compresslevel=1) as b:
    for x,y in zip_longest(reads(r1,1), reads(r2,2)):
        if x is None or y is None or x[0] != y[0]:
            raise ValueError('Mate mismatch at pair ' + str(n+1))
        a.write(x[1]); b.write(y[1]); n += 1
if n == 0:
    raise ValueError('Empty input')
print('pairs\t' + str(n))
