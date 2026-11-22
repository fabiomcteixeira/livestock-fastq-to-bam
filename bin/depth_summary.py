#!/usr/bin/env python3
"""Stream samtools depth -aa output: include zero-depth bases in denominator."""
import sys

print('contig\tpositions\tcovered_positions\tmean_depth\tpct_covered')
previous = None
n = covered = total = 0
def report(chrom, n, covered, total):
    print('{}\t{}\t{}\t{:.6f}\t{:.4f}'.format(chrom, n, covered, total/n, 100*covered/n))
for line in sys.stdin:
    chrom, pos, depth = line.split()
    if previous is not None and chrom != previous:
        report(previous, n, covered, total)
        n = covered = total = 0
    previous = chrom
    d = int(depth)
    n += 1
    total += d
    covered += d > 0
if previous is not None:
    report(previous, n, covered, total)
