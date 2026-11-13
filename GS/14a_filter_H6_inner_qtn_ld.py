#!/usr/bin/env python3
"""H6 inner-fold GWAS-ranked imputed QTN / selected-CHIP LD filter.

Uses existing inner_train SNP rankings and inner_train genotype-derived PLINK LD.
No phenotype fitting, no GWAS/clumping, no access to inner_valid/outer_test genotypes.
PLINK's sparse --ld-window-r2 threshold reports ONLY strong-LD pairs; unreported
pairs are below the threshold, NOT claimed to have an actual r2 of zero.
"""
import argparse
import csv
import gzip
import math
import os
from pathlib import Path

COUNTS = (1,2,3,4,5,10,20,30,40,50,75,100,150,200,300,400,500,750,1000)


def required(path):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        raise RuntimeError(f"File missing or empty: {p}")
    return p


def ids(path):
    p = required(path)
    result=[]
    with p.open() as fp:
        for ln, line in enumerate(fp,1):
            f=line.split()
            if len(f)!=1:
                raise RuntimeError(f"Expected one ID per line: {p}:{ln}")
            result.append(f[0])
    if len(result)!=len(set(result)):
        raise RuntimeError(f"Duplicate ID in {p}")
    return result


def ranked(path, n):
    p=required(path)
    out=[]; seen=set(); prev=-1.0
    with p.open() as fp:
        header=fp.readline().split()
        if header[:2]!=["SNP","P"] or len(header) not in (2,3):
            raise RuntimeError(f"Unexpected ranking header {header} in {p}")
        for ln,line in enumerate(fp,2):
            if len(out)==n:
                break
            fields=line.split()
            if len(fields)!=len(header):
                raise RuntimeError(f"Bad ranking row {p}:{ln}")
            snp,p_str=fields[:2]
            v=float(p_str)
            if not math.isfinite(v) or v<0 or v>1 or v<prev:
                raise RuntimeError(f"Invalid or unordered P value at {p}:{ln}")
            if snp in seen:
                raise RuntimeError(f"Duplicate ranked SNP {snp} at {p}:{ln}")
            if len(header)==3 and (snp.split(':',1)[0]!=fields[2]):
                raise RuntimeError(f"Chromosome mismatch at {p}:{ln}")
            out.append((snp,p_str)); prev=v; seen.add(snp)
    if not out:
        raise RuntimeError(f"No SNPs found in {p}")
    return out


def write(path, lines):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    with tmp.open('w',newline='') as fp:
        for line in lines:
            fp.write(line+'\n')
    os.replace(tmp,path)


def prepare(a):
    chip=ids(a.chip)
    r=ranked(a.ranked,a.pool)
    if len(set(s for s,_ in r) & set(chip)):
        raise RuntimeError('Ranked imputed candidate overlaps selected chip SNP')
    train=ids(a.train)
    valid=ids(a.valid)
    outer_test=ids(a.outer_test)
    if not 500 <= len(train) <= 600:
        raise RuntimeError(f'Unexpected inner_train sample size {len(train)}')
    if set(train)&set(valid) or set(train)&set(outer_test) or set(valid)&set(outer_test):
        raise RuntimeError('Nested-CV sample overlap')
    d=Path(a.stage);d.mkdir(parents=True,exist_ok=True)
    write(d/'selected_chip.snps',chip)
    write(d/'requested_qtn.snps',[s for s,p in r])
    write(d/'union.snps',chip+[s for s,p in r])
    write(d/'candidate_rank.tsv',['Original_Rank\tSNP\tP']+[f'{i}\t{s}\t{p}' for i,(s,p) in enumerate(r,1)])
    write(d/'train.keep',[f'{s}\t{s}' for s in train])
    print(f'[PREP] selected CHIP: {len(chip)}; original candidates in pool: {len(r)}; samples: {len(train)}')


def read_bim(path):
    result={}
    with required(path).open() as fp:
        for ln,line in enumerate(fp,1):
            x=line.split()
            if len(x)!=6: raise RuntimeError(f'Invalid BIM {path}:{ln}')
            if x[1] in result: raise RuntimeError(f'Duplicate BIM SNP {x[1]}')
            result[x[1]]=x[0]
    return result


def report(a):
    d=Path(a.stage)
    chip=ids(d/'selected_chip.snps')
    train=ids(a.train)
    ranked_rows=ranked(a.ranked,a.pool)
    staged=ids(d/'requested_qtn.snps')
    if staged != [s for s,p in ranked_rows]:
        raise RuntimeError('Candidate rank changed since extraction')
    bim=read_bim(a.bim)
    selected=set(chip)
    candidate=set(staged)
    missing=(selected|candidate)-set(bim)
    unexpected=set(bim)-(selected|candidate)
    if missing or unexpected:
        raise RuntimeError(f'PLINK marker/QC mismatch: missing={len(missing)} {sorted(missing)[:5]}, extra={len(unexpected)}')
    fam=required(a.fam)
    fam_ids=[]
    with fam.open() as fp:
        for line in fp:
            f=line.split()
            if len(f)<2 or f[0]!=f[1]:
                raise RuntimeError('Invalid FID/IID in FAM')
            fam_ids.append(f[1])
    if len(fam_ids)!=len(train) or set(fam_ids)!=set(train) or len(set(fam_ids))!=len(fam_ids):
        raise RuntimeError('PLINK genotype samples do not match inner_train exactly')
    chip_by_chr={}
    for s in chip: chip_by_chr[bim[s]]=chip_by_chr.get(bim[s],0)+1
    max_r2={};best_chip={};n_pairs=0
    # Each chromosome is processed independently to bound H6 work size.
    # Files are streamed; no whole-genome LD matrix is retained.
    for ld_path in a.ld:
        ld=required(ld_path)
        opening=gzip.open if ld.suffix=='.gz' else open
        with opening(ld,'rt') as fp:
            h=fp.readline().split()
            if not set(('SNP_A','SNP_B','R2')).issubset(h):
                raise RuntimeError(f'Bad PLINK LD header: {h}')
            ai,bi,ri=h.index('SNP_A'),h.index('SNP_B'),h.index('R2')
            for ln,line in enumerate(fp,2):
                t=line.split()
                if len(t)<=max(ai,bi,ri):raise RuntimeError(f'Bad LD row {ln}')
                snpa,snpb=t[ai],t[bi]
                if snpa in selected and snpb in candidate:
                    ch,qtn=snpa,snpb
                elif snpb in selected and snpa in candidate:
                    ch,qtn=snpb,snpa
                else:
                    continue
                if bim[qtn]!=bim[ch]:
                    raise RuntimeError(f'Interchromosomal LD unexpectedly present {qtn} {ch}')
                r=float(t[ri])
                if not math.isfinite(r) or r<0 or r>1.00001:
                    raise RuntimeError(f'Invalid LD r2 at row {ln}: {t[ri]}')
                r=min(1.0,r)
                if r+1e-8 < a.threshold:
                    raise RuntimeError(f'PLINK LD row has r2 {r} below report threshold {a.threshold}')
                n_pairs+=1
                if qtn not in max_r2 or r>max_r2[qtn] or (r==max_r2[qtn] and ch<best_chip[qtn]):
                    max_r2[qtn]=r;best_chip[qtn]=ch
    annotated=['Original_Rank\tSNP\tP\tChr\tLD_status\tMax_reported_high_r2\tTop_linked_CHIP\tKeep\tFiltered_Rank']
    retained=[];excluded=0;no_chip=0;old20_removed=0
    for i,(sid,p) in enumerate(ranked_rows,1):
        if chip_by_chr.get(bim[sid],0)==0:
            status='NO_SELECTED_CHIP_ON_CHR'; r2s='NA'; link='NA'; keep=True;no_chip+=1
        elif sid in max_r2:
            status='FILTERED_HIGH_LD';r2s=f'{max_r2[sid]:.7f}';link=best_chip[sid];keep=False
        else:
            status='BELOW_THRESHOLD';r2s=f'<{a.threshold:g}';link='NA';keep=True
        if keep:
            retained.append((sid,p,i,status,r2s,link))
        else:
            excluded+=1
            if i<=20:old20_removed+=1
        annotated.append(f'{i}\t{sid}\t{p}\t{bim[sid]}\t{status}\t{r2s}\t{link}\t{int(keep)}\t{len(retained) if keep else "NA"}')
    write(d/'candidate_ld_annotated.tsv',annotated)
    write(d/'filtered_ranked_extra_snps.tsv',['SNP\tP']+[f'{s}\t{p}' for s,p,*_ in retained])
    write(d/'excluded_high_ld.snps',[s for s,p in ranked_rows if s in max_r2])
    write(d/'filter_summary.tsv',['Item\tValue',f'Pool_size\t{len(ranked_rows)}',f'CHIP_count\t{len(chip)}',f'Inner_train_size\t{len(train)}',f'LD_threshold\t{a.threshold}',f'Filtered_high_LD\t{excluded}',f'Retained\t{len(retained)}',f'No_same_chr_chip\t{no_chip}',f'Original_top20_filtered\t{old20_removed}',f'High_LD_pairs_reported\t{n_pairs}',f'Can_make_top1000\t{int(len(retained)>=1000)}'])
    print(f'[RESULT] pool={len(ranked_rows)} retained={len(retained)} filtered={excluded}; no_same_chr_chip={no_chip}; top1000={len(retained)>=1000}')
    # One per-stage completion sentinel only after every report is committed.
    write(d/'stage.complete',[f'pool={len(ranked_rows)}\tretained={len(retained)}\tthreshold={a.threshold}'])
    print(f'[RESULT] {d / "filter_summary.tsv"}')
    print(f'[RESULT] {d / "candidate_ld_annotated.tsv"}')


def publish(a):
    stage=Path(a.stage)
    if not required(stage/'stage.complete'):
        raise RuntimeError('Stage not completed')
    f=required(stage/'filtered_ranked_extra_snps.tsv')
    rows=[]
    with f.open() as fp:
        hdr=fp.readline().split()
        if hdr!=['SNP','P']:raise RuntimeError('Bad filtered ranking header')
        for line in fp:
            v=line.split()
            if len(v)!=2:raise RuntimeError('Bad filtered ranking row')
            rows.append(v)
    if len(rows)<1000:
        raise RuntimeError(f'Cannot publish: only {len(rows)} survivors (need 1000)')
    out=Path(a.out)
    out.mkdir(parents=True,exist_ok=True)
    if (out/'complete.flag').exists():
        raise RuntimeError(f'Already published: {out}. Refusing to overwrite complete result')
    # Only publish top 1000, matching the pre-existing Step 9a matrix limit.
    write(out/'filtered_ranked_extra_snps.tsv',['SNP\tP']+[f'{s}\t{p}' for s,p in rows])
    write(out/'top1000_extra.snps',[s for s,p in rows[:1000]])
    write(out/'candidate_count_availability.tsv',['N_QTN\tAvailable'] + [f'{n}\t1' for n in COUNTS])
    write(out/'marker_sets'/'extra_0.snps',[])
    for n in COUNTS:write(out/'marker_sets'/f'extra_{n}.snps',[s for s,p in rows[:n]])
    # Copy stable report metadata into an independent, non-overwriting final folder.
    for nm in ('filter_summary.tsv','candidate_ld_annotated.tsv'):
        write(out/nm,required(stage/nm).read_text().splitlines())
    write(out/'pilot_metadata.tsv',['Item\tValue',f'Stage\t{stage.resolve()}',f'Full_filtered_ranked_count\t{len(rows)}',f'Usable_QTN_counts\t{len(COUNTS)}'])
    write(out/'complete.flag',['DONE cross-CHIP LD filter diagnostic; inner_train only'])
    print(f'[PUBLISHED] {out}; full_filtered={len(rows)}; all {len(COUNTS)} positive QTN sets available')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='mode',required=True)
    x=sub.add_parser('prepare')
    for k in ('ranked','chip','train','valid','outer-test','stage'):
        x.add_argument('--'+k,required=True)
    x.add_argument('--pool',type=int,required=True)
    y=sub.add_parser('report')
    for k in ('ranked','stage','bim','fam','train'):
        y.add_argument('--'+k,required=True)
    y.add_argument('--ld', required=True, action='append',
                   help='One gzipped PLINK sparse LD file per chromosome; repeat flag')
    y.add_argument('--pool',type=int,required=True)
    y.add_argument('--threshold',type=float,default=0.8)
    z=sub.add_parser('publish')
    for k in ('stage','out'):
        z.add_argument('--'+k,required=True)
    a=p.parse_args()
    if a.mode in ('prepare','report') and a.pool<1000:p.error('pool must be >= 1000')
    if a.mode=='report' and not 0<a.threshold<=1:p.error('threshold must be in (0,1]')
    {'prepare':prepare,'report':report,'publish':publish}[a.mode](a)

if __name__=='__main__':main()
