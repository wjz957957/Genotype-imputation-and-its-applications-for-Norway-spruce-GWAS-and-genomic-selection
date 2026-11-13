#!/usr/bin/env python3
"""Validate a posfile and write BCF region targets and exact SNP keys, streaming."""
import argparse,gzip,json,os,re
p=argparse.ArgumentParser()
for a in ('--chunk','--input','--keys','--regions','--stats'): p.add_argument(a,required=True)
a=p.parse_args()
reader=gzip.open if a.input.endswith('.gz') else open
n=ts=tv=other=0; last=-1; prev=None; duplicates=0; max_pos=0
with reader(a.input,'rt') as inp,open(a.keys+'.tmp','w') as keys,open(a.regions+'.tmp','w') as regions:
 for line_number,line in enumerate(inp,1):
  if not line.strip(): continue
  parts=line.split();
  if len(parts)<4: raise RuntimeError(f'{a.input}:{line_number} fewer than four columns')
  chrom,pos,ref,alt=parts[:4]
  if chrom!=a.chunk: raise RuntimeError(f'{a.input}:{line_number} unexpected contig {chrom} != {a.chunk}')
  pos=int(pos)
  if pos<1 or pos<last: raise RuntimeError(f'{a.input}:{line_number} unsorted/invalid POS {pos}; previous {last}')
  if pos==last: duplicates+=1
  last=pos;max_pos=max(max_pos,pos)
  ref,alt=ref.upper(),alt.upper()
  if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or ref==alt or ',' in alt:
   raise RuntimeError(f'{a.input}:{line_number} non-biallelic SNP {chrom}:{pos}:{ref}:{alt}')
  n+=1
  if (ref,alt) in {('A','G'),('G','A'),('C','T'),('T','C')}:ts+=1
  else:tv+=1
  keys.write(f'{chrom}\t{pos}\t{ref}\t{alt}\n')
  if pos!=prev:
   regions.write(f'{chrom}\t{pos}\t{pos}\n')
   prev=pos
os.replace(a.keys+'.tmp',a.keys);os.replace(a.regions+'.tmp',a.regions)
with open(a.stats+'.tmp','w') as out:json.dump(dict(chunk=a.chunk,candidate_rows=n,ts=ts,tv=tv,same_position_adjacent=duplicates,max_local_pos=max_pos),out)
os.replace(a.stats+'.tmp',a.stats)
print(a.chunk,'candidate',n,'same-position adjacent',duplicates,flush=True)
