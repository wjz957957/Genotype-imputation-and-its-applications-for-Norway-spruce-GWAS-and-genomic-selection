#!/usr/bin/env python3
"""Create a 2-column CHROM/POS targets file for bcftools view -T from exact keys."""
import argparse, os, sys
p=argparse.ArgumentParser()
p.add_argument('--keys', required=True)
p.add_argument('--targets', required=True)
a=p.parse_args()
tmp=a.targets+'.tmp'
rows=0; unique=0; last_chrom=None; last_pos=0
with open(a.keys) as src, open(tmp,'w') as dst:
    for n,line in enumerate(src,1):
        f=line.rstrip('\n').split('\t')
        if len(f)!=4:
            raise SystemExit(f'{a.keys}:{n}: expected 4 tab-separated columns')
        chrom,pos,ref,alt=f
        pos=int(pos)
        if chrom!=last_chrom:
            if last_chrom is not None and chrom < last_chrom:
                raise SystemExit(f'{a.keys}:{n}: chromosome order decreases')
            last_chrom=chrom; last_pos=0
        if pos < last_pos:
            raise SystemExit(f'{a.keys}:{n}: positions are not sorted')
        rows+=1
        if pos != last_pos:
            dst.write(f'{chrom}\t{pos}\n')
            unique+=1
        last_pos=pos
os.replace(tmp,a.targets)
print(f'keys={rows} unique_target_positions={unique} targets={a.targets}', flush=True)
