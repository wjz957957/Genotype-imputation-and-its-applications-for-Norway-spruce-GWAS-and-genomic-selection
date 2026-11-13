#!/usr/bin/env python3
"""Outer-train selected-CHIP vs imputed-QTN cross-panel LD filtering.

Uses existing outer_train GWAS/clumping rank and outer_train genotype LD.
No access to outer_test genotype during filtering; no GWAS reanalysis.
PLINK's sparse --ld-window-r2 threshold reports ONLY strong-LD pairs; unreported
pairs are below the threshold, NOT claimed to have an actual r2 of zero.
"""
import argparse
import gzip
import math
import os
from pathlib import Path

COUNTS = (1,2,3,4,5,10,20,30,40,50,75,100,150,200,300,400,500,600,700)


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


def ranked(path, n=None):
    p=required(path)
    out=[]; seen=set(); prev=-1.0
    with p.open() as fp:
        header=fp.readline().split()
        if header[:2]!=["SNP","P"] or len(header) not in (2,3):
            raise RuntimeError(f"Unexpected ranking header {header} in {p}")
        for ln,line in enumerate(fp,2):
            if n is not None and len(out)==n:
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



def rank_nonchip(a):
    """Remove exact overlaps with the frozen selected chip, preserving P order."""
    chip = set(ids(a.chip))
    rows = ranked(a.ranked)
    out_rows = []
    n_overlap = 0
    for snp, p in rows:
        if snp in chip:
            n_overlap += 1
        else:
            out_rows.append((snp, p))
    if not out_rows:
        raise RuntimeError("No non-chip candidates remain after exact-overlap removal")
    write(a.out, ["SNP\tP"] + [f"{s}\t{p}" for s, p in out_rows])
    write(a.summary, ["Item\tValue", f"Input_clump_indices\t{len(rows)}",
                      f"Exact_chip_overlaps_removed\t{n_overlap}",
                      f"Nonchip_candidates_retained\t{len(out_rows)}"])
    print(f"PASS rank: input={len(rows)} exact_chip_overlaps={n_overlap} nonchip={len(out_rows)}")


def prepare(a):
    chip=ids(a.chip)
    r=ranked(a.ranked,a.pool)
    if len(set(s for s,_ in r) & set(chip)):
        raise RuntimeError('Ranked imputed candidate overlaps selected chip SNP')
    train=ids(a.train)
    test=ids(a.test)
    if not 700 <= len(train) <= 900:
        raise RuntimeError(f'Unexpected outer_train sample size {len(train)}')
    if len(train)+len(test)!=904 or set(train)&set(test):
        raise RuntimeError('Outer_train/test sample overlap or sample count mismatch')
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
    # Each chromosome is processed independently to bound work size.
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
    write(d/'filter_summary.tsv',['Item\tValue',f'Pool_size\t{len(ranked_rows)}',f'CHIP_count\t{len(chip)}',f'Outer_train_size\t{len(train)}',f'LD_threshold\t{a.threshold}',f'Filtered_high_LD\t{excluded}',f'Retained\t{len(retained)}',f'No_same_chr_chip\t{no_chip}',f'Original_top20_filtered\t{old20_removed}',f'High_LD_pairs_reported\t{n_pairs}',f'Can_make_selected_count\t{int(len(retained)>=a.target)}'])
    print(f'[RESULT] pool={len(ranked_rows)} retained={len(retained)} filtered={excluded}; no_same_chr_chip={no_chip}; selected_{a.target}={len(retained)>=a.target}')
    # One per-stage completion sentinel only after every report is committed.
    write(d/'stage.complete',[f'pool={len(ranked_rows)}\tretained={len(retained)}\tthreshold={a.threshold}'])
    print(f'[RESULT] {d / "filter_summary.tsv"}')
    print(f'[RESULT] {d / "candidate_ld_annotated.tsv"}')


def publish(a):
    stage=Path(a.stage)
    required(stage/'stage.complete')
    rows=[]
    with required(stage/'filtered_ranked_extra_snps.tsv').open() as fp:
        hdr=fp.readline().split()
        if hdr!=['SNP','P']:
            raise RuntimeError(f'Invalid filtered ranking header: {hdr}')
        for line in fp:
            f=line.split()
            if len(f)!=2:raise RuntimeError('Invalid filtered ranking row')
            rows.append(f)
    if a.target not in COUNTS and a.target != 0:
        raise RuntimeError(f'Invalid target QTN count: {a.target}')
    if len(rows)<a.target:
        raise RuntimeError(f'Only {len(rows)} survivors, selected {a.target}')
    out=Path(a.out)
    out.mkdir(parents=True,exist_ok=True)
    if (out/'complete.flag').exists():
        raise RuntimeError(f'Already completed: {out}. Refusing to overwrite.')
    write(out/'filtered_ranked_extra_snps.tsv', ['SNP\tP']+[f'{s}\t{p}' for s,p in rows])
    write(out/'candidate_ld_annotated.tsv', required(stage/'candidate_ld_annotated.tsv').read_text().splitlines())
    write(out/'filter_summary.tsv', required(stage/'filter_summary.tsv').read_text().splitlines()+[
        f'Best_Extra\t{a.target}', f'Full_filtered_ranked_count\t{len(rows)}'])
    write(out/'marker_sets'/'extra_0.snps', [])
    for n in COUNTS:
        if n<=len(rows):
            write(out/'marker_sets'/f'extra_{n}.snps', [s for s,p in rows[:n]])
    # The selected QTN list must exist even when target=0 (empty list).
    selected_list=out/'marker_sets'/f'extra_{a.target}.snps'
    if not selected_list.is_file() or sum(1 for _ in selected_list.open())!=a.target:
        raise RuntimeError('Selected QTN file missing or incorrect length')
    write(out/'best_extra.txt', [str(a.target)])
    write(out/'selected_top_QTN.snps', [s for s,p in rows[:a.target]])
    write(out/'candidate_count_availability.tsv', ['N_QTN\tAvailable\tReason'] + [
        f'{n}\t{int(n<=len(rows))}\t{("AVAILABLE" if n<=len(rows) else "INSUFFICIENT_SURVIVORS")}'
        for n in (0,)+COUNTS])
    write(out/'complete.flag', ['DONE outer_train only, r2=0.8'])
    print(f'[PUBLISHED] {out}; selected QTN={a.target}, available={len(rows)}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='mode',required=True)
    r=sub.add_parser('rank')
    for k in ('ranked','chip','out','summary'):
        r.add_argument('--'+k,required=True)
    x=sub.add_parser('prepare')
    for k in ('ranked','chip','train','test','stage'):
        x.add_argument('--'+k,required=True)
    x.add_argument('--pool',type=int,required=True)
    y=sub.add_parser('report')
    for k in ('ranked','stage','bim','fam','train'):
        y.add_argument('--'+k,required=True)
    y.add_argument('--ld', required=True, action='append',
                   help='One gzipped PLINK sparse LD file per chromosome; repeat flag')
    y.add_argument('--pool',type=int,required=True)
    y.add_argument('--threshold',type=float,default=0.8)
    y.add_argument('--target',type=int,required=True)
    z=sub.add_parser('publish')
    for k in ('stage','out'):
        z.add_argument('--'+k,required=True)
    z.add_argument('--target',type=int,required=True)
    a=p.parse_args()
    if a.mode in ('prepare','report') and a.pool<1:p.error('pool must be >=1')
    if a.mode=='report' and not 0<a.threshold<=1:p.error('threshold must be in (0,1]')
    {'rank':rank_nonchip,'prepare':prepare,'report':report,'publish':publish}[a.mode](a)

if __name__=='__main__':main()
