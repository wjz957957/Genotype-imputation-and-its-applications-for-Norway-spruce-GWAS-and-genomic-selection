#!/usr/bin/env python3
"""Verify every outer-fold final RAW dataset against Step16d and Step17a."""
import csv
from pathlib import Path

WORK = Path.cwd()
OPT = WORK/'09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_allfolds/optimization/best_extra_settings.tsv'
MARKER_ROOT = WORK/'10_outer_test/06_crosschip_LD_r2_0p8_allfolds'
RAW_ROOT = WORK/'10_outer_test/07_final_genotypes_crosschip_r2_0p8'
CV_ROOT = WORK/'02_CV'
CHIP_ROOT = WORK/'10_outer_test/01_chip_selection'


def lines(path):
    if not path.exists():
        raise ValueError(f'Missing file: {path}')
    with path.open() as f:
        return [x.strip() for x in f if x.strip()]


def read_raw(path, n_expected, ids_expected, selected_snps):
    if not path.is_file():
        raise ValueError(f'Missing RAW file: {path}')
    with path.open() as f:
        header = f.readline().split()
        if header[:6] != ['FID','IID','PAT','MAT','SEX','PHENOTYPE']:
            raise ValueError(f'Unexpected PLINK RAW header in {path}: {header[:6]}')
        if len(header) != n_expected + 6:
            raise ValueError(f'Wrong marker count: {path}, expected={n_expected}, got={len(header)-6}')
        raw_snps = [x.rsplit('_', 1)[0] for x in header[6:]]
        if len(raw_snps)!=len(set(raw_snps)) or set(raw_snps)!=set(selected_snps):
            raise ValueError(f'RAW markers do not match selected SNP IDs: {path}')
        observed = []
        for rowno, row in enumerate(f, 2):
            # Restrict RAM usage; only need FID/IID and genotype column count.
            fields=row.split()
            if len(fields)!=len(header):
                raise ValueError(f'Column mismatch {path}:{rowno} expected={len(header)} got={len(fields)}')
            if fields[0]!=fields[1]:
                raise ValueError(f'Inconsistent FID/IID at {path}:{rowno}')
            observed.append(fields[1])
        if len(observed)!=len(ids_expected) or set(observed)!=set(ids_expected) or len(set(observed))!=len(observed):
            raise ValueError(f'Wrong training/test individuals in {path}')
        return header[6:], observed


def main():
    with OPT.open(newline='') as f:
        rows = list(csv.DictReader(f,delimiter='\t'))
    if len(rows)!=20:
        raise ValueError(f'Expected 20 selection records, got {len(rows)}')
    seen=set()
    for r in rows:
        outer,trait,n=r['Outer'],r['Trait'],int(r['Best_Extra'])
        key=(outer,trait)
        if outer not in [f'{i:02d}' for i in range(1,11)] or trait not in ('BB','H6') or key in seen:
            raise ValueError(f'Invalid/duplicate outer/trait: {key}')
        seen.add(key)
        d=RAW_ROOT/f'outer_{outer}'/trait
        q=MARKER_ROOT/f'outer_{outer}'/trait
        sample_dir=CV_ROOT/f'outer_{outer}'
        chip=CHIP_ROOT/f'outer_{outer}'/'selected'/f'{trait}_best_chip.snps'
        if not (d/'complete.flag').is_file() or not (q/'complete.flag').is_file():
            raise ValueError(f'Completion flag missing: {key}')
        cids=lines(chip)
        qids=lines(q/'marker_sets'/f'extra_{n}.snps') if n else []
        if len(qids)!=n or len(cids)!=len(set(cids)) or len(qids)!=len(set(qids)) or set(cids)&set(qids):
            raise ValueError(f'Invalid CHIP/QTN marker sets for {key}')
        train=lines(sample_dir/'outer_train.samples')
        test=lines(sample_dir/'outer_test.samples')
        for split,samples in (('train',train),('test',test)):
            chip_h, chip_s = read_raw(d/f'chip_{split}.raw',len(cids),samples,cids)
            if n:
                extra_h, extra_s = read_raw(d/f'extra_{split}.raw',n,samples,qids)
                if chip_s!=extra_s:
                    raise ValueError(f'CHIP/EXTRA RAW sample orders differ: {key} {split}')
            elif (d/f'extra_{split}.raw').exists():
                raise ValueError(f'Unexpected extra RAW at zero-QTN: {key} {split}')
        with (d/'final_genotype_manifest.tsv').open(newline='') as f:
            mani=list(csv.DictReader(f,delimiter='\t'))
        if (len(mani)!=1 or int(mani[0]['Best_Extra'])!=n or
                int(mani[0]['N_chip'])!=len(cids) or
                int(mani[0]['N_total'])!=len(cids)+n or
                int(mani[0]['N_train'])!=len(train) or
                int(mani[0]['N_test'])!=len(test)):
            raise ValueError(f'Manifest inconsistent: {key}')
        print(f'[PASS] outer_{outer} {trait}: CHIP={len(cids)} QTN={n} N_train={len(train)} N_test={len(test)}')
    assert len(seen)==20
    print('VERIFIED ALL 20 FINAL DATASETS')

if __name__=='__main__':
    main()
