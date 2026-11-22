#!/usr/bin/env python3
"""QC decisions and tables; Python >=3.6, standard library only."""
import sys, json, csv, re, zipfile
from pathlib import Path

def qc(animal, run, stage, paths):
    mates=[]
    for path in sorted(paths):
        with zipfile.ZipFile(path) as z:
            names=[n for n in z.namelist() if n.endswith('/fastqc_data.txt')]
            if len(names)!=1: raise ValueError('Expected one fastqc_data.txt: '+path)
            text=z.read(names[0]).decode('utf-8')
        status={}
        count=None
        for line in text.splitlines():
            if line.startswith('>>') and not line.startswith('>>END_MODULE'):
                name, flag=line[2:].split('\t')
                status[name]=flag.upper()
            if line.startswith('Total Sequences\t'): count=int(line.split('\t')[1])
        mates.append(dict(file=Path(path).name, reads=count,
            adapters=status.get('Adapter Content','MISSING'),
            quality=status.get('Per base sequence quality','MISSING')))
    valid=len(mates)==2 and all(m['reads'] is not None and m['adapters'] in ('PASS','WARN','FAIL') and m['quality'] in ('PASS','WARN','FAIL') for m in mates)
    if valid: valid=mates[0]['reads']==mates[1]['reads']
    return dict(animal=animal, run=run, kind=stage+'_qc', mates=mates,valid=valid,
                total_reads=sum(m['reads'] or 0 for m in mates),
                adapter_pass=valid and all(m['adapters']!='FAIL' for m in mates),
                quality_pass=valid and all(m['quality']!='FAIL' for m in mates))

def trim(animal,run,log):
    text=Path(log).read_text()
    patterns=['Input Read Pairs','Both Surviving','Forward Only Surviving','Reverse Only Surviving','Dropped']
    values=[]
    for name in patterns:
        m=re.search(re.escape(name)+r':\s*(\d+)',text)
        if not m: raise ValueError('Missing trimming count: '+name)
        values.append(int(m.group(1)))
    raw,paired,u1,u2,dropped=values
    if raw!=paired+u1+u2+dropped: raise ValueError('Trimming counts inconsistent')
    return dict(animal=animal,run=run,kind='trim',raw_pairs=raw,paired_after_trim=paired,
                unpaired_R1=u1,unpaired_R2=u2,reads_after_trim=2*paired+u1+u2)

def bam(animal,run,kind,count,depth):
    with open(depth) as f: rows=list(csv.DictReader(f,delimiter='\t'))
    n=sum(int(r['positions']) for r in rows)
    mean=sum(int(r['positions'])*float(r['mean_depth']) for r in rows)/n if n else None
    return dict(animal=animal,run=run,kind=kind,mapped_primary_reads=int(Path(count).read_text()),
                reference_positions=n,mean_depth=mean)

def write_table(name,fields,rows):
    with open(name,'w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def report(manifest,paths):
    with open(manifest) as f: entries=list(csv.DictReader(f,delimiter='\t'))
    records={}
    for path in paths:
        d=json.loads(Path(path).read_text());key=(d['run'],d['kind'])
        if key in records: raise ValueError('Duplicate record '+str(key))
        records[key]=d
    stats=[]; review=[]; animals={}
    for r in entries:
        a,run=r['animal_id'],r['run']; animals.setdefault(a,[]).append(run)
        raw=records.get((run,'raw_qc'),{});t=records.get((run,'trim'),{})
        post=records.get((run,'trim_qc'),{});b=records.get((run,'run'),{})
        state='RAW_QC_ONLY'
        if not raw.get('valid'): state='HOLD_RAW_INVALID'
        elif not raw.get('adapter_pass') and not t: state='HOLD_RAW_ADAPTERS'
        if t: state='TRIMMED'
        if post:
            state='TRIM_QC_PASS' if post.get('adapter_pass') and post.get('quality_pass') and post.get('total_reads',0)>0 else 'HOLD_TRIM_QC'
        if b: state='MAPPED'
        row=dict(animal_id=a,run=run,status=state,raw_reads=raw.get('total_reads','NA'),
            raw_pairs=t.get('raw_pairs',raw.get('total_reads',0)//2 if raw.get('valid') else 'NA'),
            reads_after_trim=t.get('reads_after_trim','NA'),paired_after_trim=t.get('paired_after_trim','NA'),
            unpaired_R1=t.get('unpaired_R1','NA'),unpaired_R2=t.get('unpaired_R2','NA'),
            reads_sent_to_mapping=2*t['paired_after_trim'] if b and t else 'NA',
            mapped_primary_reads=b.get('mapped_primary_reads','NA'),mean_depth=b.get('mean_depth','NA'))
        stats.append(row)
        for q in [raw,post]:
            for m in q.get('mates',[]):
                review.append(dict(animal_id=a,run=run,stage=q['kind'],file=m['file'],adapter_content=m['adapters'],per_base_quality=m['quality'],run_status=state))
    fields=['animal_id','run','status','raw_reads','raw_pairs','reads_after_trim','paired_after_trim','unpaired_R1','unpaired_R2','reads_sent_to_mapping','mapped_primary_reads','mean_depth']
    write_table('run_statistics.tsv',fields,stats)
    write_table('qc_review.tsv',['animal_id','run','stage','file','adapter_content','per_base_quality','run_status'],review)
    final=[]
    for a,runs in sorted(animals.items()):
        done=[r for r in runs if (r,'run') in records]
        m=records.get((a,'merged'),{});d=records.get((a,'final'),{})
        final.append(dict(animal_id=a,expected_runs=len(runs),mapped_runs=len(done),pending_runs=','.join(sorted(set(runs)-set(done))),
            status='FINAL' if d else ('MERGED' if m else 'WAITING'),
            merged_mapped_primary_reads=m.get('mapped_primary_reads','NA'),merged_mean_depth=m.get('mean_depth','NA'),
            final_mapped_primary_reads=d.get('mapped_primary_reads','NA'),final_mean_depth=d.get('mean_depth','NA')))
    write_table('animal_statistics.tsv',['animal_id','expected_runs','mapped_runs','pending_runs','status','merged_mapped_primary_reads','merged_mean_depth','final_mapped_primary_reads','final_mean_depth'],final)

if __name__=='__main__':
    mode=sys.argv[1]; args=sys.argv[2:]
    if mode=='qc': result=qc(*args[:3],args[3:])
    elif mode=='trim': result=trim(*args)
    elif mode=='bam': result=bam(*args)
    elif mode=='report': report(args[0],args[1:]);sys.exit(0)
    else: raise ValueError('Unknown command')
    print(json.dumps(result,sort_keys=True))
