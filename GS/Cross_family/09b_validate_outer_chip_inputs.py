#!/usr/bin/env python3
"""Training-only genotype input/output audits for the outer-chip GWAS stage.

This stage extracts ALL chip SNPs; it does not perform GWAS, rank SNPs,
clump LD, or decide the final SNP identity for a selected chip strategy.
"""
import argparse
import csv
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def ids(path, ncol=1):
    path = Path(path)
    require(path.is_file(), f'missing file: {path}')
    ans=[]
    for j, line in enumerate(path.read_text().splitlines(),1):
        a=line.split()
        require(len(a)==ncol, f'{path}:{j}: expected {ncol} fields, found {len(a)}')
        if ncol==2:
            require(a[0]==a[1], f'{path}:{j}: FID != IID')
        ans.append(a[0])
    require(bool(ans), f'empty IDs: {path}')
    require(len(ans)==len(set(ans)), f'duplicate IDs: {path}')
    return ans


def pheno_ids(path):
    path=Path(path)
    require(path.is_file(), f'missing phenotype: {path}')
    seen=[]
    for j, line in enumerate(path.read_text().splitlines(),1):
        a=line.split()
        require(len(a)==3 and a[0]==a[1],f'invalid phenotype FID IID Y: {path}:{j}')
        v=float(a[2])
        import math
        require(math.isfinite(v),f'nonfinite phenotype: {path}:{j}')
        seen.append(a[0])
    require(len(seen)==len(set(seen)),f'duplicate phenotype IDs: {path}')
    return seen


def bim_snp_ids(path):
    path=Path(path)
    require(path.is_file(), f'missing BIM {path}')
    ans=[]
    for j, line in enumerate(path.read_text().splitlines(),1):
        a=line.split()
        require(len(a)==6, f'invalid BIM columns: {path}:{j}')
        require(a[0] in {str(k) for k in range(1,13)},
                f'non-chr1-12 BIM position: {path}:{j}')
        ans.append(a[1])
    require(len(ans)==len(set(ans)), f'duplicate BIM SNPs: {path}')
    return ans


def fam_ids(path):
    path=Path(path)
    require(path.is_file(),f'missing FAM {path}')
    ans=[]
    for j, line in enumerate(path.read_text().splitlines(),1):
        a=line.split()
        require(len(a)==6 and a[0]==a[1], f'bad FAM: {path}:{j}')
        ans.append(a[0])
    require(len(ans)==len(set(ans)),f'duplicate FAM IDs: {path}')
    return ans


def read_frozen(path,o):
    with Path(path).open(newline='') as f: r=list(csv.DictReader(f,delimiter='\t'))
    require(len(r)==2,f'frozen selection must have exactly 2 rows: {path}')
    require({x['Trait'] for x in r}=={'BB','H6'},'frozen selection traits are not BB/H6')
    for line in r:
        require(line['Outer']==o,f'frozen selection outer != {o}')
        require(line['Best_Chip'].startswith('chip_') and line['Best_Chip'][5:].isdigit(),
                'invalid frozen chip strategy')
        require(int(line['Best_QTN']) in {0,1,2,3,4,5,10,20,30,40,50,75,100,150,200,300,400,500,600,700},
                'invalid frozen QTN number')
    return r


def validate(work,outer,output,prefix,phase):
    w=Path(work); o=f'{int(outer):02d}'; d=w/'Cross_family/results/09_outer_inputs'/f'outer_{o}'
    flag=d/'prepare.complete'
    require(flag.is_file() and 'status=PASS' in flag.read_text().splitlines(),
            f'Stage09a input is not PASS: {flag}')
    require('outer_test_phenotype_read=NO' in flag.read_text().splitlines(),
            f'Stage09a test leakage guard missing: {flag}')
    keep=ids(d/'outer_train.keep',2)
    test=ids(d/'outer_test.ids',1)
    chips=ids(d/'chip_ALL.snps',1)
    require(len(chips)==43267, f'expected 43267 full chip SNPs, got {len(chips)}')
    require(500<=len(keep)<904 and len(keep)+len(test)==904,
            f'outer_train/outer_test count invalid: {len(keep)}+{len(test)}')
    require(set(keep).isdisjoint(test),'outer train/test samples overlap')
    read_frozen(d/'frozen_inner_selection.tsv',o)
    for t in ['BB','H6']:
        p=pheno_ids(d/f'{t}.train.pheno')
        require(p==keep, f'{t} phenotype sample order differs from outer_train.keep')
    master=w/'08_imputed_GWAS/00_reference/all_904_genotypes'
    fam=fam_ids(str(master)+'.fam')
    require(len(fam)==904 and set(keep)|set(test)==set(fam),
            'all training and test individuals must partition master FAM')
    if phase=='post':
        pref=Path(prefix)
        for ext in ['bed','bim','fam']:
            require(Path(str(pref)+'.'+ext).is_file() and Path(str(pref)+'.'+ext).stat().st_size>0,
                    f'missing PLINK output {pref}.{ext}')
        b=bim_snp_ids(str(pref)+'.bim')
        ob=fam_ids(str(pref)+'.fam')
        require(len(b)==43267 and set(b)==set(chips),
                f'PLINK extraction failed to preserve exactly 43267 chip IDs: {len(b)}')
        require(len(ob)==len(keep) and set(ob)==set(keep),
                f'PLINK extraction sample IDs mismatch; observed {len(ob)}, expected {len(keep)}')
        require(set(ob).isdisjoint(test),'leakage: outer-test individual in extraction')
        print(f'PASS post outer_{o}: training={len(ob)} chip_SNPs={len(b)} outer_test_excluded={len(test)}')
    else:
        print(f'PASS pre outer_{o}: training={len(keep)} test={len(test)} chip_ALL={len(chips)}')
    return len(keep),len(test),len(chips)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--work',required=True);p.add_argument('--outer',required=True)
    p.add_argument('--out',required=True);p.add_argument('--prefix',required=True)
    p.add_argument('--phase',choices=['pre','post'],required=True)
    a=p.parse_args()
    validate(a.work,a.outer,a.out,a.prefix,a.phase)

if __name__=='__main__':
    main()
