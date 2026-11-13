#!/usr/bin/env python3
import argparse,csv,json,os,collections
p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args()
r=a.root
with open(f'{r}/config/chromosome_lengths.tsv') as f: lengths={x['chromosome']:int(x['length_bp']) for x in csv.DictReader(f,delimiter='\t')}
with open(f'{r}/config/manifest.tsv') as f: manifest=list(csv.DictReader(f,delimiter='\t'))
by_chr=collections.defaultdict(collections.Counter);whole=collections.Counter();af=collections.defaultdict(collections.Counter);maf=collections.defaultdict(collections.Counter)
summary=[];problems=[]
for x in manifest:
 chunk=x['chunk'];chrom=x['chrom']
 with open(f'{r}/position_stats/{chunk}.json') as f:ps=json.load(f)
 with open(f'{r}/per_file/{chunk}.json') as f:bs=json.load(f)
 cnt=bs['counts'];v={'candidate_rows':ps['candidate_rows'],'candidate_ts':ps['ts'],'candidate_tv':ps['tv'],'candidate_same_pos_adjacent':ps['same_position_adjacent'],'matched_records':cnt.get('matched_records',0),'unmatched_candidate_keys':cnt.get('unmatched_candidate_keys',0),'allele_mismatch_records':cnt.get('allele_mismatch_records',0),'called_diploid_gt':cnt.get('called_diploid_gt',0),'het_gt':cnt.get('het_gt',0),'missing_gt':cnt.get('missing_gt',0),'other_gt':cnt.get('other_gt',0),'non_diploid_gt':cnt.get('non_diploid_gt',0),'sites_with_called_alleles':cnt.get('sites_with_called_alleles',0)}
 summary.append(dict(chunk=chunk,chrom=chrom,**v))
 for k,value in v.items():by_chr[chrom][k]+=value;whole[k]+=value
 for i,n in bs['af_bins'].items():af[chrom][int(i)]+=n;af['genome'][int(i)]+=n
 for i,n in bs['maf_bins'].items():maf[chrom][int(i)]+=n;maf['genome'][int(i)]+=n
 if v['candidate_rows']!=v['matched_records'] or v['unmatched_candidate_keys'] or v['allele_mismatch_records'] or v['candidate_same_pos_adjacent']:
  problems.append(chunk)
def row(label,v):
 ts=v['candidate_ts'];tv=v['candidate_tv'];called=v['called_diploid_gt'];nmiss=v['missing_gt'];bad=v['other_gt']+v['non_diploid_gt'];eligible=called+nmiss+bad
 return dict(chromosome=label,length_bp=lengths.get(label,''),candidate_rows=v['candidate_rows'],candidate_ts=ts,candidate_tv=tv,candidate_titv=ts/tv if tv else 'NA',candidate_density_per_Mb=v['candidate_rows']/(lengths[label]/1e6) if label in lengths else 'NA',matched_records=v['matched_records'],unmatched_candidate_keys=v['unmatched_candidate_keys'],allele_mismatch_records=v['allele_mismatch_records'],heterozygosity=v['het_gt']/called if called else 'NA',missingness=nmiss/eligible if eligible else 'NA',called_diploid_gt=called,het_gt=v['het_gt'],missing_gt=nmiss,other_gt=v['other_gt'],non_diploid_gt=v['non_diploid_gt'],sites_with_called_alleles=v['sites_with_called_alleles'])
def write(path,rows):
 if not rows:return
 with open(path,'w',newline='') as out:
  w=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write(f'{r}/summary/chunk_summary.tsv',summary)
write(f'{r}/summary/chromosome_summary.tsv',[row(c,by_chr[c]) for c in sorted(lengths)])
write(f'{r}/summary/overall_summary.tsv',[row('genome',whole)])
write(f'{r}/summary/frequency_spectrum.tsv',[dict(scope=scope,type=kind,lower=i*.05,upper=min((i+1)*.05,1),count=n) for kind,table in [('AF',af),('MAF',maf)] for scope,cnt in sorted(table.items()) for i,n in sorted(cnt.items())])
with open(f'{r}/summary/validation.txt','w') as out:
 out.write('Candidates expected: 2741098686\n')
 out.write(f"Observed rows: {whole['candidate_rows']}\n")
 out.write(f"Matched records: {whole['matched_records']}\n")
 out.write(f"Chunks requiring review: {','.join(problems) or 'none'}\n")
 out.write('Exact uniqueness across separate chunk contigs relies on unique contig labels and no overlapping genomic representations.\n')
if problems:raise SystemExit(f'VALIDATION FAILED: {len(problems)} chunks have incomplete/ambiguous matching or adjacent repeated positions; inspect summary/validation.txt')
if whole['candidate_rows']!=2741098686:raise SystemExit('VALIDATION FAILED: candidate count differs from manuscript')
print('VALIDATED all candidate SNP records match exact BCF variants',whole['candidate_rows'])
