#!/usr/bin/env python3
"""Outer-training-only phenotype and full chip roster reconstruction.

No outer-test phenotypes and no old 10-fold input files are read.  Each
outer-training individual appears exactly once in the 3 inner-validation
phenotype files; cross-check all 3 inner-training phenotype files.
"""
import argparse
import csv
import os
from pathlib import Path

ALLOWED={0,1,2,3,4,5,10,20,30,40,50,75,100,150,200,300,400,500,600,700}
WORK_SUB='Cross_family'

def require(test, msg):
    if not test:
        raise RuntimeError(msg)

def read_ids(path):
    p=Path(path);require(p.is_file(),f'missing sample/SNP file: {p}')
    vals=[x.strip() for x in p.read_text().splitlines()]
    require(len(vals)>0 and all(x and len(x.split())==1 for x in vals),f'invalid ID list {p}')
    require(len(vals)==len(set(vals)),f'duplicate IDs: {p}')
    return vals

def read_pheno(path):
    p=Path(path);require(p.is_file(),f'missing phenotype {p}')
    vals={}
    for line_no,line in enumerate(p.read_text().splitlines(),1):
        s=line.split()
        require(len(s)==3 and s[0]==s[1],f'bad FID IID Y: {p}:{line_no}')
        require(s[1] not in vals,f'duplicate phenotype ID: {p}:{line_no}')
        import math
        try: v=float(s[2])
        except ValueError as e: raise RuntimeError(f'non-numeric phenotype: {p}:{line_no}') from e
        require(math.isfinite(v),f'nonfinite phenotype: {p}:{line_no}')
        vals[s[1]]=(s[2],v)
    require(vals,f'empty phenotype: {p}')
    return vals

def table(path):
    p=Path(path);require(p.is_file(),f'missing table: {p}')
    with p.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def atomic_lines(path,vals):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(str(p)+'.tmp');tmp.write_text('\n'.join(vals)+'\n');os.replace(tmp,p)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--work',required=True)
    ap.add_argument('--outer',required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args();o=f'{int(a.outer):02d}'
    require(1 <= int(o) <= 5,'outer fold out of range')
    w=Path(a.work);cf=w/WORK_SUB;out=Path(a.out)
    selection=cf/'results/08_inner_qtn_selection'
    select_flag=selection/'complete.flag'
    require(select_flag.is_file() and 'status=PASS' in select_flag.read_text().splitlines(),
            f'Stage08 selection not PASS: {select_flag}')
    require('selection_rule=ALL_THREE_SIGNED_POSITIVE_ELSE_MAX_MEAN_SIGNED_R' in
            select_flag.read_text().splitlines(),
            f'Stage08 selection is not the signed-r fallback version: {select_flag}')
    chips=table(cf/'config/best_chip_by_outer_trait.tsv')
    selected=table(selection/'selected_qtn_counts.tsv')
    chosen=[]
    for trait in ('BB','H6'):
        c=[r for r in chips if r.get('Outer')==o and r.get('Trait')==trait]
        q=[r for r in selected if r.get('outer')==o and r.get('trait')==trait]
        require(len(c)==len(q)==1,f'chip or QTN selection missing/duplicate: {o}/{trait}')
        strategy=c[0]['Best_Chip'];n=int(q[0]['n_qtn'])
        require(strategy.startswith('chip_') and strategy[5:].isdigit() and
                int(strategy[5:])>=1 and n in ALLOWED,
                f'invalid strategy/count: {o}/{trait}')
        # Both selection paths are valid: strict 3/3 SIGNED r improvement,
        # or (only when no count passes 3/3) the maximum mean signed r.
        record=q[0]
        mode=record.get('selection_mode','')
        eligible=record.get('eligible','')
        npositive=int(record.get('n_eligible_positive_counts','-1'))
        nfolds=int(record.get('n_folds_positive','-1'))
        selected_rule=record.get('selection_rule','')
        valid_choice=(
            (mode=='THREE_OF_THREE_SIGNED_IMPROVEMENT' and n>0 and
             eligible=='YES' and npositive>=1 and nfolds==3)
            or (mode=='FALLBACK_MAX_MEAN_SIGNED_R' and
                eligible=='NO' and npositive==0 and 0<=nfolds<=2)
        )
        require(valid_choice and
                selected_rule=='all3_signed_improvement_else_max_mean_signed_Pearson',
                f'QTN selection inconsistent with signed-r fallback rule: '
                f'{o}/{trait}, n={n}, eligible={eligible}, mode={mode}, '
                f'positive={npositive}, nfolds={nfolds}')
        chosen.append((trait,strategy,n))
    cv=cf/'cv_samples_v1'/f'outer_{o}'
    train=read_ids(cv/'outer_train.samples');test=read_ids(cv/'outer_test.samples')
    require(len(train)>=500 and len(test)>=50,f'unexpected outer sample sizes {o}')
    require(not (set(train)&set(test)),f'outer train/test overlap outer_{o}')
    fam=w/'08_imputed_GWAS/00_reference/all_904_genotypes.fam'
    require(fam.is_file(),f'missing master fam {fam}')
    observed={}
    for line in fam.read_text().splitlines():
        t=line.split();require(len(t)>=2 and t[0]==t[1],f'bad FAM entry: {line}')
        require(t[1] not in observed,f'duplicate master FAM: {t[1]}')
        observed[t[1]]=1
    require(len(observed)==904,f'expected 904 individuals in master, found {len(observed)}')
    require(set(train)|set(test)==set(observed),f'outer train/test not an exact master partition outer_{o}')
    fullchips=None
    for i in range(1,4):
        inner=f'inner_{i:02d}'
        cvdir=cv/inner
        tr=read_ids(cvdir/'inner_train.samples')
        va=read_ids(cvdir/'inner_valid.samples')
        require(not (set(tr)&set(va)) and set(tr)|set(va)==set(train),
                f'invalid inner split outer_{o}/{inner}')
        chip_base=cf/'results/03_inner_chip_gwas_clumping'/f'outer_{o}'/inner/'marker_sets'
        for trait in ('BB','H6'):
            c=read_ids(chip_base/f'{trait}_chip_ALL.snps')
            require(len(c)==43267,f'full chip list not 43267 at {chip_base}')
            if fullchips is None: fullchips=c
            require(c==fullchips,f'chip_ALL differs between outer_{o} inner folds/traits')
    atomic_lines(out/'chip_ALL.snps',fullchips)
    atomic_lines(out/'outer_train.keep',[f'{sid}\t{sid}' for sid in train])
    atomic_lines(out/'outer_test.ids',test)   # exclusion audit, no phenotypes
    for trait,strategy,n in chosen:
        # each train individual appears in exactly one of 3 inner-validation files
        all_values={}
        for i in range(1,4):
            inner=f'inner_{i:02d}';p=cf/'results/02_inner_phenotypes'/f'outer_{o}'/inner
            va=read_pheno(p/f'{trait}.valid.pheno')
            expected=set(read_ids(cv/inner/'inner_valid.samples'))
            require(set(va)==expected,f'inner validation phenotype samples mismatch: {p}/{trait}')
            require(not (set(va)&set(all_values)),f'duplicate OOF validation phenotype ID: {o}/{trait}')
            all_values.update(va)
        require(set(all_values)==set(train),f'outer-training phenotype reconstruction incomplete {o}/{trait}')
        for i in range(1,4):
            inner=f'inner_{i:02d}';p=cf/'results/02_inner_phenotypes'/f'outer_{o}'/inner
            tr=read_pheno(p/f'{trait}.train.pheno')
            expected=set(read_ids(cv/inner/'inner_train.samples'))
            require(set(tr)==expected,f'inner-training phenotype samples mismatch: {p}/{trait}')
            for sid,(raw,value) in tr.items():
                require(abs(value-all_values[sid][1])<=1e-9*max(1,abs(value)),
                        f'inconsistent phenotype for {o}/{trait}/{sid}')
        atomic_lines(out/f'{trait}.train.pheno',
                     [f'{sid}\t{sid}\t{all_values[sid][0]}' for sid in train])
    atomic_lines(out/'frozen_inner_selection.tsv',
        ['Outer\tTrait\tBest_Chip\tBest_QTN']+
        [f'{o}\t{trait}\t{strategy}\t{n}' for trait,strategy,n in chosen])
    atomic_lines(out/'prepare.complete',[
        'status=PASS',f'outer={o}',f'outer_training={len(train)}',f'outer_testing={len(test)}',
        'sample_partition_master=PASS','inner_validation_pheno_merged=PASS',
        'inner_training_pheno_consistency=PASS','full_chip_consistency=PASS',
        'outer_test_phenotype_read=NO'])
    print(f'PASS outer_{o}: n_train={len(train)} n_test={len(test)} full_chip={len(fullchips)} '
          f'BB={chosen[0][1]}+{chosen[0][2]}QTN H6={chosen[1][1]}+{chosen[1][2]}QTN')
if __name__=='__main__': main()
