#!/usr/bin/env python3
"""Stream BCF once using bcftools view -T, then exact-match candidate CHROM/POS/REF/ALT keys.
GT-derived QC is accumulated only for exact candidate SNP matches.
"""
import argparse, collections, json, os, subprocess, sys, tempfile, time

p=argparse.ArgumentParser()
for x in ('--chunk','--bcf','--targets','--keys','--output'):
    p.add_argument(x, required=True)
p.add_argument('--progress-every', type=int, default=1_000_000)
a=p.parse_args()

# Load one candidate key at a time; keys are sorted by CHROM/POS.
def key_reader(path):
    with open(path) as fh:
        for n,line in enumerate(fh,1):
            f=line.rstrip('\n').split('\t')
            if len(f)!=4:
                raise RuntimeError(f'{path}:{n}: expected 4 columns')
            c,pos,r,alt=f
            yield (c,int(pos),r.upper(),alt.upper())

keys=iter(key_reader(a.keys))
def next_key():
    try: return next(keys)
    except StopIteration: return None
key=next_key()
counts=collections.Counter(); af_bins=collections.Counter(); maf_bins=collections.Counter(); samples=[]
start=time.monotonic(); queried=0; last_progress=start

def bin20(x):
    return min(19, int(x*20))

cmd1=['bcftools','view','-T',a.targets,'-Ou',a.bcf]
cmd2=['bcftools','query','-f','%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n','-']
print('COMMAND:', ' '.join(cmd1), flush=True)
# Avoid stderr pipe deadlocks; capture each command's diagnostics in temporary files.
with tempfile.TemporaryFile(mode='w+b') as err1, tempfile.TemporaryFile(mode='w+b') as err2:
    p1=subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=err1)
    p2=subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE,
                        stderr=err2, text=True, bufsize=1024*1024)
    p1.stdout.close()
    try:
        for line in p2.stdout:
            f=line.rstrip('\n').split('\t')
            if len(f)<4:
                raise RuntimeError('Malformed bcftools query output')
            c,pos,ref,alt=f[:4]; pos=int(pos); queried+=1
            counts['bcf_records_at_candidate_positions']+=1
            site=(c,pos,ref.upper(),alt.upper())
            # Advance candidate keys that precede this BCF record.
            while key is not None and (key[0],key[1]) < (c,pos):
                counts['unmatched_candidate_keys']+=1
                key=next_key()
            if key is None:
                counts['unmatched_bcf_records']+=1
            elif (key[0],key[1]) > (c,pos):
                counts['unmatched_bcf_records']+=1
            elif site != key:
                counts['allele_mismatch_records']+=1
                if ',' in alt:
                    counts['multiallelic_bcf_records_at_candidate_positions']+=1
                # Keep key: another BCF record at same coordinate may match it.
            else:
                counts['matched_records']+=1
                gt=f[4:]
                if not samples:
                    samples=[collections.Counter() for _ in gt]
                    counts['sample_count']=len(gt)
                if len(gt)!=len(samples):
                    raise RuntimeError('Inconsistent sample count across records')
                alt_alleles=called_alleles=0
                for j,g in enumerate(gt):
                    s=samples[j]; s['candidate_sites_observed']+=1
                    alleles=g.replace('|','/').split('/')
                    if len(alleles)!=2:
                        s['non_diploid_gt']+=1; counts['non_diploid_gt']+=1
                        continue
                    if '.' in alleles:
                        s['missing_gt']+=1; counts['missing_gt']+=1
                        continue
                    if alleles[0] not in ('0','1') or alleles[1] not in ('0','1'):
                        s['other_gt']+=1; counts['other_gt']+=1
                        continue
                    dosage=int(alleles[0])+int(alleles[1])
                    alt_alleles+=dosage; called_alleles+=2
                    s['called_diploid_gt']+=1; counts['called_diploid_gt']+=1
                    if dosage==1:
                        s['het_gt']+=1; counts['het_gt']+=1
                    else:
                        s['hom_gt']+=1; counts['hom_gt']+=1
                if called_alleles:
                    af=alt_alleles/called_alleles
                    af_bins[bin20(af)]+=1
                    maf_bins[bin20(min(af,1-af))]+=1
                    counts['sites_with_called_alleles']+=1
                else:
                    counts['sites_all_gt_missing_or_invalid']+=1
                key=next_key()
            if a.progress_every>0 and queried % a.progress_every==0:
                elapsed=time.monotonic()-start
                print(f'PROGRESS queried_records={queried} matched={counts["matched_records"]} '
                      f'unmatched_keys={counts["unmatched_candidate_keys"]} elapsed_seconds={elapsed:.1f}',
                      flush=True)
        p2.stdout.close()
        rc2=p2.wait(); rc1=p1.wait()
    except BaseException:
        p2.kill(); p1.kill()
        p2.wait(); p1.wait()
        raise
    err1.seek(0); msg1=err1.read().decode(errors='replace')
    err2.seek(0); msg2=err2.read().decode(errors='replace')
    if rc1 or rc2:
        raise RuntimeError(f'bcftools failed view={rc1} query={rc2}\n{msg1[-4000:]}\n{msg2[-4000:]}')

while key is not None:
    counts['unmatched_candidate_keys']+=1
    key=next_key()

result=dict(
    chunk=a.chunk,
    counts=dict(counts),
    af_bins={str(i):af_bins[i] for i in range(20)},
    maf_bins={str(i):maf_bins[i] for i in range(20)},
    sample_counts=[dict(s) for s in samples],
    methodology='Single-pass bcftools view -T streaming targets; exact CHROM,POS,REF,ALT match; diploid GT only for AF/heterozygosity; fully/partially missing GT counted missing; non-diploid and other GT separately counted; AF/MAF 0.05 bins',
    elapsed_seconds=round(time.monotonic()-start,2),
)
os.makedirs(os.path.dirname(os.path.abspath(a.output)),exist_ok=True)
tmp=a.output+'.tmp'
with open(tmp,'w') as out:
    json.dump(result,out,separators=(',',':'))
os.replace(tmp,a.output)
print(f'COMPLETE {a.chunk} matched={counts["matched_records"]} '
      f'unmatched_candidate_keys={counts["unmatched_candidate_keys"]} '
      f'allele_mismatch={counts["allele_mismatch_records"]} '
      f'elapsed_seconds={result["elapsed_seconds"]}',flush=True)
