#!/usr/bin/env python3
"""Small report only: does not delete files or decide which reads to map."""
import csv
import sys
import zipfile
from pathlib import Path

folder = Path(sys.argv[1] if len(sys.argv)>1 else 'results/01_fastqc_raw')
files = sorted(folder.glob('*/*_fastqc.zip'))
if not files:
    raise SystemExit('No FastQC ZIP files found')
with open('fastqc_raw_summary.tsv','w',newline='') as out:
    w = csv.writer(out,delimiter='\t')
    w.writerow(['run','file','adapter_content','per_base_quality','overrepresented_sequences','all_modules'])
    for p in files:
        with zipfile.ZipFile(p) as z:
            name = next(n for n in z.namelist() if n.endswith('/summary.txt'))
            lines=[line.split('\t') for line in z.read(name).decode().splitlines()]
            statuses={line[1]:line[0] for line in lines}
            w.writerow([p.parent.name,lines[0][2],statuses.get('Adapter Content','NOT_REPORTED'),
                        statuses.get('Per base sequence quality','NOT_REPORTED'),
                        statuses.get('Overrepresented sequences','NOT_REPORTED'),
                        '; '.join(k+'='+v for k,v in statuses.items())])
print('{} FastQC reports summarized in fastqc_raw_summary.tsv'.format(len(files)))
if len(files)!=38:
    print('NOTICE: expected 38 reports for the complete current manifest.')
