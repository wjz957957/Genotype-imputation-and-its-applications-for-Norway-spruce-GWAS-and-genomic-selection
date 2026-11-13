#!/usr/bin/env python3
"""Import successfully completed outer_01 pilots into official all-fold r2=0.8 namespace.
Read-only old results. Does not import their large stage BED/LD directories.
"""
import argparse
from pathlib import Path
import csv
import os

COUNTS=(1,2,3,4,5,10,20,30,40,50,75,100,150,200,300,400,500,600,700)

def lines(p):
    p=Path(p)
    if not p.is_file() or p.stat().st_size==0:raise RuntimeError(f'Missing/empty file: {p}')
    return p.read_text().splitlines()

def ids(p):
    x=lines(p)
    if any(len(s.split())!=1 for s in x) or len(set(x))!=len(x):
        raise RuntimeError(f'Invalid SNP list {p}')
    return x

def atomic_write(path,content):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(''.join(line+'\n' for line in content))
    tmp.replace(path)

def main():
    p=argparse.ArgumentParser()
    for x in ('source','dest','ranked','chip'):p.add_argument('--'+x,required=True)
    p.add_argument('--chip-count',required=True,type=int)
    p.add_argument('--training-size',required=True,type=int)
    a=p.parse_args()
    src=Path(a.source);dest=Path(a.dest)
    if not (src/'complete.flag').is_file():raise RuntimeError('Original pilot is not complete')
    new_rank=lines(src/'filtered_ranked_extra_snps.tsv')
    if new_rank[0].split()!=['SNP','P']:raise RuntimeError('Bad filtered ranking header')
    selected=[]
    for line in new_rank[1:]:
        cols=line.split()
        if len(cols)!=2:raise RuntimeError('Bad filtered ranking row')
        selected.append(cols[0])
    if len(set(selected))!=len(selected) or len(selected)<700:
        raise RuntimeError('Pilot ranking lacks >=700 unique SNPs')
    if len(ids(a.chip))!=a.chip_count:raise RuntimeError('CHIP count mismatch')
    if set(selected)&set(ids(a.chip)):raise RuntimeError('Pilot QTN overlaps selected CHIP')
    actual_ranking=lines(a.ranked)
    if actual_ranking[0].split()[:2]!=['SNP','P']:raise RuntimeError('Bad original ranking')
    all_ranked={x.split()[0] for x in actual_ranking[1:]}
    if not set(selected).issubset(all_ranked):raise RuntimeError('Pilot SNP not in current fold ranking')
    top=(src/'top700_extra.snps') if (src/'top700_extra.snps').is_file() else (src/'top1000_extra.snps')
    if ids(top)[:700]!=selected[:700]:raise RuntimeError('Pilot top SNPs differ from ranking')
    summary={r[0]:r[1] for r in csv.reader(lines(src/'filter_summary.tsv')[1:],delimiter='\t')}
    for key,val in (('CHIP_count',a.chip_count),('Inner_train_size',a.training_size)):
        if int(summary[key])!=val:raise RuntimeError(f'Pilot metadata mismatch: {key}')
    if float(summary['LD_threshold'])!=0.8:raise RuntimeError('Incorrect pilot threshold')
    if dest.joinpath('complete.flag').exists():
        if ids(dest/'top700_extra.snps')!=selected[:700]:raise RuntimeError('Existing output differs')
        print('[SKIP] Already imported:',dest)
        return
    # Publish output files first; completion last.
    atomic_write(dest/'filtered_ranked_extra_snps.tsv',new_rank)
    atomic_write(dest/'top700_extra.snps',selected[:700])
    for name in ('filter_summary.tsv','candidate_ld_annotated.tsv'):
        atomic_write(dest/name,lines(src/name))
    atomic_write(dest/'marker_sets'/'extra_0.snps',[])
    for n in COUNTS:atomic_write(dest/'marker_sets'/f'extra_{n}.snps',selected[:n])
    atomic_write(dest/'candidate_count_availability.tsv',['N_QTN\tAvailable']+[f'{n}\t1' for n in COUNTS])
    atomic_write(dest/'import_provenance.tsv',['Field\tValue',f'Source\t{src}',f'Raw_ranked\t{a.ranked}',f'CHIP_selected_count\t{a.chip_count}',f'Inner_train_size\t{a.training_size}'])
    atomic_write(dest/'complete.flag',['DONE verified pilot reuse: cross CHIP LD r2=0.8'])
    print('[IMPORTED]',src,'->',dest,'QTN=',len(selected))

if __name__=='__main__':main()
