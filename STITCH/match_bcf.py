#!/usr/bin/env python3
"""Subset indexed BCF via exact CHROM/POS targets, then match REF/ALT and accumulate GT-derived QC.
No additional Python dependencies; GT output is piped, never materialized.
"""
import argparse,collections,csv,json,math,os,subprocess,sys
p=argparse.ArgumentParser()
for x in ('--chunk','--bcf','--regions','--keys','--output'):p.add_argument(x,required=True)
a=p.parse_args()
# Re-query all selected sites; includes REF ALT to verify exact allele identity.
cmd1=['bcftools','view','-R',a.regions,'-Ou',a.bcf]
cmd2=['bcftools','query','-f','%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n','-']
p1=subprocess.Popen(cmd1,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
p2=subprocess.Popen(cmd2,stdin=p1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,bufsize=1024*1024)
p1.stdout.close()
counts=collections.Counter(); freq=collections.Counter(); maf=collections.Counter(); sample=[]
def next_key(handle):
 line=handle.readline()
 if not line:return None
 c,pos,ref,alt=line.rstrip('\n').split('\t');return (c,int(pos),ref,alt)
def genobin(p):return min(19,int(p*20))  # 20 bins, last includes 1.0
with open(a.keys) as kf:
 key=next_key(kf);previous_site=None
 for line in p2.stdout:
  f=line.rstrip('\n').split('\t')
  if len(f)<4:raise RuntimeError('Malformed query line')
  c,pos,ref,alt=f[:4]; pos=int(pos)
  counts['bcf_records_at_target_coordinates']+=1
  if ',' in alt:counts['multiallelic_bcf_at_target']+=1
  # Ordered merge over contig/position; exact REF ALT match, not just coordinate.
  while key is not None and (key[0],key[1])<(c,pos):
   counts['unmatched_candidate_keys']+=1;key=next_key(kf)
  if key is None:counts['unmatched_bcf_records']+=1;continue
  if (key[0],key[1])>(c,pos):
   counts['unmatched_bcf_records']+=1;continue
  if (c,pos,ref.upper(),alt.upper())!=key:
   counts['allele_mismatch_records']+=1
   # Need to keep key for other BCF records at same position
   continue
  counts['matched_records']+=1; key=next_key(kf)
  gt=f[4:]
  if not sample:sample=[collections.Counter() for _ in gt];counts['sample_count']=len(gt)
  if len(gt)!=len(sample):raise RuntimeError('Inconsistent sample count')
  alt_alleles=called_alleles=0
  for j,g in enumerate(gt):
   s=sample[j];s['eligible_snp_sites']+=1
   alleles=g.replace('|','/').split('/')
   if len(alleles)!=2: s['non_diploid']+=1;counts['non_diploid_gt']+=1;continue
   if '.' in alleles:s['missing_gt']+=1;counts['missing_gt']+=1;continue
   if alleles[0] not in ('0','1') or alleles[1] not in ('0','1'):
    s['other_gt']+=1;counts['other_gt']+=1;continue
   dosage=int(alleles[0])+int(alleles[1]); alt_alleles+=dosage;called_alleles+=2
   s['called_diploid_gt']+=1;counts['called_diploid_gt']+=1
   if dosage==1:s['het_gt']+=1;counts['het_gt']+=1
   else:s['hom_gt']+=1;counts['hom_gt']+=1
  if called_alleles:
   af=alt_alleles/called_alleles
   freq[genobin(af)]+=1;maf[genobin(min(af,1-af))]+=1
   counts['sites_with_called_alleles']+=1
  else:counts['sites_all_gt_missing_or_invalid']+=1
 while key is not None:
  counts['unmatched_candidate_keys']+=1;key=next_key(kf)
p2.stdout.close()
err2=p2.stderr.read()
err1=p1.stderr.read()
rc2=p2.wait();rc1=p1.wait()
if rc1 or rc2:raise RuntimeError(f'bcftools failed view={rc1} query={rc2} stderr={(err1.decode(errors="replace")+str(err2))[-6000:]}')
if counts['unmatched_candidate_keys'] or counts['allele_mismatch_records'] or counts['unmatched_bcf_records']:
 print('WARNING: candidate / BCF match diagnostics:',dict(counts),file=sys.stderr)
result=dict(chunk=a.chunk,counts=dict(counts),af_bins={str(i):freq[i] for i in range(20)},maf_bins={str(i):maf[i] for i in range(20)},sample_counts=[dict(s) for s in sample],
 methodology='Exact CHROM,POS,REF,ALT match; diploid GT=0/0,0/1,1/1; ./.,0/. missing; non-diploid and other allele GT separate; 0.05 AF bins')
with open(a.output+'.tmp','w') as out:json.dump(result,out)
os.replace(a.output+'.tmp',a.output)
print(a.chunk,'matched',counts['matched_records'],'not_found',counts['unmatched_candidate_keys'],'allele_mismatch',counts['allele_mismatch_records'],flush=True)
