#!/usr/bin/env python3
import csv
import sys
import re
from pathlib import Path

source, target = sys.argv[1:]
with open(source) as f:
    reader = csv.DictReader(f, delimiter='\t')
    required = {'animal_id','run','R1','R2','header_format','library_id'}
    if not required.issubset(reader.fieldnames or []):
        raise SystemExit('Missing manifest columns')
    fields = reader.fieldnames
    rows = list(reader)
if not rows:
    raise SystemExit('Empty manifest')
runs, paths, animals = set(), set(), {}
for r in rows:
    for key in ('animal_id','run'):
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',r[key]):
            raise SystemExit('Invalid identifier: ' + r[key])
    if r['library_id'] and not re.fullmatch(r'[A-Za-z0-9_.-]+',r['library_id']):
        raise SystemExit('Invalid library_id')
    if r['run'] in runs:
        raise SystemExit('Duplicate run: ' + r['run'])
    runs.add(r['run'])
    identity = r['animal_id']
    if r['animal_id'] in animals and animals[r['animal_id']] != identity:
        raise SystemExit('Inconsistent sample metadata')
    animals[r['animal_id']] = identity
    if r['header_format'] not in ('standard','hash_mate'):
        raise SystemExit('Unknown header format')
    for key in ('R1','R2'):
        p = Path(r[key])
        if not p.is_absolute() or not p.is_file() or not p.stat().st_size or "'" in str(p) or '\n' in str(p):
            raise SystemExit('Invalid FASTQ path: ' + str(p))
        resolved = str(p.resolve())
        if resolved in paths:
            raise SystemExit('FASTQ reused: ' + resolved)
        paths.add(resolved)
        r[key] = resolved
out = Path(target)
out.parent.mkdir(parents=True, exist_ok=True)
with out.open('w') as f:
    w = csv.DictWriter(f, fields, delimiter='\t'); w.writeheader(); w.writerows(rows)
print('{} animals; {} runs; {} FASTQs'.format(len(animals),len(runs),len(paths)))
