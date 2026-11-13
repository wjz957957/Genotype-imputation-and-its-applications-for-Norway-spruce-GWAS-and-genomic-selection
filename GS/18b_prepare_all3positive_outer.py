#!/usr/bin/env python3
"""Build outer-genotype inputs for all-3-inner-positive count selection.

Never reruns GWAS, PLINK or LD. Reuses verified Step17 CHIP RAWs and
selects QTN columns by SNP ID from already extracted Step17 extra RAWs.
Only supports new count <= old count; any other situation fails closed.
"""
import argparse
import csv
import os
import shutil
from pathlib import Path

OLD_COUNTS = {0,1,2,3,4,5,10,20,30,40,50,75,100,150,200,300,400,500,600,700}
EXPECTED = {(f'{i:02d}',t) for i in range(1,11) for t in ('BB','H6')}
OPT_REL = Path('09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_allfolds/optimization')
OLD_GENO_REL = Path('10_outer_test/07_final_genotypes_crosschip_r2_0p8')
NEW_GENO_REL = Path('10_outer_test/08_final_genotypes_crosschip_r2_0p8_all3positive')
QTN_REL = Path('10_outer_test/06_crosschip_LD_r2_0p8_allfolds')
CV_REL = Path('02_CV')


def records(path):
    with path.open(newline='',encoding='utf-8-sig') as fp:
        return list(csv.DictReader(fp,delimiter='\t'))


def nonblank_lines(path):
    if not path.is_file():
        raise RuntimeError(f'Missing file: {path}')
    return [x.strip() for x in path.read_text().splitlines() if x.strip()]


def get_header(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing/empty PLINK RAW: {path}')
    with path.open() as fp:
        h=fp.readline().split()
    if h[:6] != ['FID','IID','PAT','MAT','SEX','PHENOTYPE'] or len(h)<7:
        raise RuntimeError(f'Invalid PLINK RAW header: {path}')
    snps=[x.rsplit('_',1)[0] for x in h[6:]]
    if len(snps)!=len(set(snps)):
        raise RuntimeError(f'Duplicate SNPs in {path}')
    return h,snps


def validate_sample_size(path, n):
    # Input Step17 is already fully verified; only count rows here.
    with path.open() as fp:
        count=sum(1 for _ in fp)-1
    if count!=n:
        raise RuntimeError(f'Wrong sample count {path}: {count} != {n}')


def subset_extra_raw(source, dest, selected, expected_ids, old_count):
    """Stream only required marker columns, respecting explicit selected order."""
    with source.open() as src, dest.open('w',newline='') as dst:
        h=src.readline().split()
        if h[:6] != ['FID','IID','PAT','MAT','SEX','PHENOTYPE']:
            raise RuntimeError(f'Bad RAW header: {source}')
        snps=[x.rsplit('_',1)[0] for x in h[6:]]
        if len(snps)!=old_count or len(set(snps))!=len(snps):
            raise RuntimeError(f'Invalid old extra SNP count: {source}')
        pos={s:i+6 for i,s in enumerate(snps)}
        missing=[s for s in selected if s not in pos]
        if missing:
            raise RuntimeError(f'New selection is not subset of old QTN RAW {source}: {missing[:5]}')
        indices=list(range(6))+[pos[s] for s in selected]
        dst.write('\t'.join(h[i] for i in indices)+'\n')
        for lineno,(line, expected) in enumerate(zip(src, expected_ids),2):
            f=line.split()
            if len(f)!=len(h) or f[0]!=expected or f[1]!=expected:
                raise RuntimeError(f'Wrong OLD QTN RAW row/order: {source}:{lineno} expected ID={expected}')
            dst.write('\t'.join(f[i] for i in indices)+'\n')
        # zip stops at the shorter iterator. Extra rows must be rejected.
        if src.readline():
            raise RuntimeError(f'Unexpected extra RAW rows: {source}')
    with dest.open() as f:
        if sum(1 for _ in f)-1 != len(expected_ids):
            raise RuntimeError(f'Not all expected samples written to {dest}')


def symlink_checked(source, dest):
    if not source.is_file():
        raise RuntimeError(f'Missing symlink source: {source}')
    dest.symlink_to(source.resolve())


def extract_selection(opt):
    selected=records(opt/'three_fold_positive_rule/best_extra_all3positive.tsv')
    original=records(opt/'best_extra_settings.tsv')
    if len(selected)!=20 or len(original)!=20:
        raise RuntimeError('Expected exactly 20 old and new best-selection rows')
    orig={}
    for r in original:
        key=r['Outer'],r['Trait']
        if key in orig or key not in EXPECTED: raise RuntimeError(f'Invalid original key {key}')
        orig[key]=r
    out={}
    for r in selected:
        key=r['Outer'],r['Trait']
        if key in out or key not in EXPECTED: raise RuntimeError(f'Invalid new key {key}')
        if key not in orig or r['Best_Chip'] != orig[key]['Best_Chip']:
            raise RuntimeError(f'CHIP setting differs between optimizations: {key}')
        old=int(orig[key]['Best_Extra']); new=int(r['New_Best_Extra'])
        if int(r['Original_Best_Extra'])!=old or old not in OLD_COUNTS or new not in OLD_COUNTS:
            raise RuntimeError(f'Inconsistent old/new best numbers: {key}')
        if new>old:
            raise RuntimeError(f'Cannot subset an enlarged set: {key} new={new} old={old}; requires PLINK')
        out[key]=(r['Best_Chip'],old,new)
    if set(out)!=EXPECTED or set(orig)!=EXPECTED:
        raise RuntimeError('Incomplete outer×trait keys')
    return out


def prepare(work,key,sel):
    outer,trait=key
    chip_setting,old,new=sel
    old_dir=work/OLD_GENO_REL/f'outer_{outer}'/trait
    qtn_dir=work/QTN_REL/f'outer_{outer}'/trait
    new_dir=work/NEW_GENO_REL/f'outer_{outer}'/trait
    old_manifest=records(old_dir/'final_genotype_manifest.tsv')
    if len(old_manifest)!=1 or not (old_dir/'complete.flag').is_file():
        raise RuntimeError(f'Step17 input is incomplete: {key}')
    m=old_manifest[0]
    if (int(m['Best_Extra'])!=old or m['Best_Chip']!=chip_setting or
        m['Trait']!=trait or m['Outer']!=outer):
        raise RuntimeError(f'Step17 manifest does not match old optimization: {key}')
    n_chip=int(m['N_chip']); n_train=int(m['N_train']); n_test=int(m['N_test'])
    if n_chip<1 or n_train+n_test!=904:
        raise RuntimeError(f'Invalid manifest sample/chip counts: {key}')
    chip_ids=nonblank_lines(work/'10_outer_test/01_chip_selection'/f'outer_{outer}'/'selected'/f'{trait}_best_chip.snps')
    if len(chip_ids)!=n_chip or len(chip_ids)!=len(set(chip_ids)):
        raise RuntimeError(f'Invalid CHIP set: {key}')
    cv=work/CV_REL/f'outer_{outer}'
    ids={'train': nonblank_lines(cv/'outer_train.samples'), 'test': nonblank_lines(cv/'outer_test.samples')}
    if len(ids['train'])!=n_train or len(ids['test'])!=n_test:
        raise RuntimeError(f'CV sample count mismatch: {key}')
    selected=(nonblank_lines(qtn_dir/'marker_sets'/f'extra_{new}.snps') if new else [])
    old_selected=(nonblank_lines(qtn_dir/'marker_sets'/f'extra_{old}.snps') if old else [])
    if (len(selected)!=new or len(set(selected))!=new or len(old_selected)!=old or
        len(set(old_selected))!=old or selected!=old_selected[:new] or
        set(selected)&set(chip_ids)):
        raise RuntimeError(f'QTN rank prefix, duplicate or CHIP overlap problem: {key}')
    if not (qtn_dir/'complete.flag').is_file():
        raise RuntimeError(f'Missing complete outer LD filter: {key}')
    for split in ('train','test'):
        h,snps=get_header(old_dir/f'chip_{split}.raw')
        if len(snps)!=n_chip or set(snps)!=set(chip_ids):
            raise RuntimeError(f'CHIP RAW does not match list: {key} {split}')
        # Older Step17 verified all SNPs/rows already. No large RAW reread.
    if new_dir.exists():
        flag=new_dir/'complete.flag'
        if not flag.is_file():
            raise RuntimeError(f'Incomplete existing target directory, inspect manually: {new_dir}')
        prior=records(new_dir/'final_genotype_manifest.tsv')
        if len(prior)!=1 or int(prior[0]['Best_Extra'])!=new or not (new_dir/'chip_train.raw').is_file():
            raise RuntimeError(f'Old output conflicts with current selection: {new_dir}')
        print(f'[SKIP] outer_{outer} {trait} N={new}')
        return
    new_dir.parent.mkdir(parents=True,exist_ok=True)
    staging=new_dir.parent/f'.{trait}.step18_staging_{os.getpid()}'
    if staging.exists(): raise RuntimeError(f'Staging already exists {staging}')
    staging.mkdir()
    try:
        for split in ('train','test'):
            symlink_checked(old_dir/f'chip_{split}.raw',staging/f'chip_{split}.raw')
            if new:
                source=old_dir/f'extra_{split}.raw'
                if new==old:
                    h,snps=get_header(source)
                    if len(snps)!=old or set(snps)!=set(old_selected):
                        raise RuntimeError(f'Original RAW selection differs from marker set: {key}')
                    symlink_checked(source,staging/f'extra_{split}.raw')
                else:
                    subset_extra_raw(source,staging/f'extra_{split}.raw',selected,ids[split],old)
        with (staging/'selected_top_QTN.snps').open('w') as fp:
            for s in selected: fp.write(s+'\n')
        with (staging/'final_genotype_manifest.tsv').open('w',newline='') as fp:
            writer=csv.DictWriter(fp,fieldnames=['Outer','Trait','Best_Chip','N_chip','Best_Extra','N_total','N_train','N_test','Old_Best_Extra','Source'],delimiter='\t',lineterminator='\n')
            writer.writeheader()
            writer.writerow(dict(Outer=outer,Trait=trait,Best_Chip=chip_setting,N_chip=n_chip,
                Best_Extra=new,N_total=n_chip+new,N_train=n_train,N_test=n_test,
                Old_Best_Extra=old,Source=str(old_dir.relative_to(work))))
        (staging/'complete.flag').write_text('DONE Step18b, all-three-positive selection\n')
        staging.rename(new_dir)
        print(f'[PASS] outer_{outer} {trait}: old QTN={old}, new={new}, CHIP={n_chip}, train={n_train}, test={n_test}')
    finally:
        if staging.exists(): shutil.rmtree(staging)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--work',default='.',help='GS-reviewer project root')
    args=ap.parse_args()
    work=Path(args.work).resolve()
    choices=extract_selection(work/OPT_REL)
    n_zero=n_same=n_sub=0
    for key in sorted(choices):
        _,old,new=choices[key]
        n_zero+=(new==0)
        n_same+=(new==old)
        n_sub+=(0<new<old)
        prepare(work,key,choices[key])
    print(f'VERIFIED/PREPARED ALL 20 DATASETS: zero={n_zero}, same={n_same}, smaller_nonzero={n_sub}')

if __name__=='__main__':
    main()
