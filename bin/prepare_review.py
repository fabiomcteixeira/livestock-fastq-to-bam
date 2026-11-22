#!/usr/bin/env python3
"""Freeze adapter and review settings, preserving stable cached paths/content."""
import sys,re,hashlib,csv,json
from pathlib import Path

def fasta(path,allow_empty=False):
    records=[];name=None;seq=[]
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: records.append((name,''.join(seq)))
            name=line[1:];seq=[]
        else:
            if name is None: raise ValueError('FASTA sequence without header: '+str(path))
            seq.append(line.upper())
    if name is not None: records.append((name,''.join(seq)))
    if not records and not allow_empty: raise ValueError('Empty adapter file')
    for name,seq in records:
        if not name or not seq or not re.fullmatch('[ACGTN]+',seq): raise ValueError('Invalid adapter '+name)
    return records

def main(base,folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    records=fasta(base)+fasta('additional_adapters.fa',True)
    names=[n for n,s in records]
    if len(names)!=len(set(names)): raise ValueError('Duplicate adapter names')
    content=''.join('>'+n+'\n'+s+'\n' for n,s in records)
    digest=hashlib.sha256(content.encode()).hexdigest()
    target=folder/('adapters_'+digest+'.fa')
    if not target.exists():target.write_text(content)
    with open('trim_overrides.tsv') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    seen=set()
    for r in rows:
        if r['run'] in seen: raise ValueError('Duplicate override')
        seen.add(r['run'])
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',r['run']):raise ValueError('Invalid run')
        if not re.fullmatch(r'(?:[A-Z]+:[0-9:]+)(?: +[A-Z]+:[0-9:]+)*',r['options']):raise ValueError('Invalid trimming options')
    print(str(target.resolve()))
if __name__=='__main__':main(*sys.argv[1:])
