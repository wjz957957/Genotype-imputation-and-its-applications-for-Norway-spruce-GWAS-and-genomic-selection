#!/usr/bin/env python3
"""Pair STITCH position files and original BCF files by chunk ID."""
import argparse,csv,glob,os,re,sys
p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--pos-root',required=True);p.add_argument('--bcf-root',required=True);a=p.parse_args()
pattern=re.compile(r'^PA_chr(\d{2})_(\d+)$')
pos={}; bcf={}
for f in glob.glob(os.path.join(a.pos_root,'PA_chr*','pos.tabulated.txt*')):
 if not (f.endswith('.txt') or f.endswith('.txt.gz')):continue
 k=os.path.basename(os.path.dirname(f))
 if pattern.fullmatch(k):
  if k in pos:raise SystemExit('Duplicate position file '+k)
  pos[k]=f
for f in glob.glob(os.path.join(a.bcf_root,'PA_chr*.bcf')):
 k=os.path.basename(f)[:-4]
 if pattern.fullmatch(k):
  if k in bcf:raise SystemExit('Duplicate BCF '+k)
  bcf[k]=f
if not pos or not bcf:raise SystemExit(f'Input discovery failed pos={len(pos)} bcf={len(bcf)}')
if pos.keys()!=bcf.keys():
 print('Position only:',sorted(pos.keys()-bcf.keys()),file=sys.stderr)
 print('BCF only:',sorted(bcf.keys()-pos.keys()),file=sys.stderr)
 raise SystemExit('Incomplete pairs; refusing to run')
for k,f in bcf.items():
 if not os.access(f,os.R_OK):raise SystemExit('Unreadable BCF '+f)
for k,f in pos.items():
 if not os.access(f,os.R_OK):raise SystemExit('Unreadable position file '+f)
os.makedirs(os.path.join(a.root,'config'),exist_ok=True)
with open(os.path.join(a.root,'config','manifest.tsv'),'w',newline='') as out:
 w=csv.writer(out,delimiter='\t');w.writerow(['chunk','chrom','part','position_file','bcf_file'])
 for k in sorted(pos,key=lambda s:tuple(map(int,pattern.fullmatch(s).groups()))):
  chrom,part=pattern.fullmatch(k).groups();w.writerow([k,'chr'+chrom,part,pos[k],bcf[k]])
print('Paired files:',len(pos),'chromosomes:',len({pattern.fullmatch(k).group(1) for k in pos}))
