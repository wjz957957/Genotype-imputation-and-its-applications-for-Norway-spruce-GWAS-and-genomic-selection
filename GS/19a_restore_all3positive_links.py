#!/usr/bin/env python3
"""Restore missing Step18 symlinks to already transferred Step17 RAW files.
Never overwrites a file; checks each fold's QTN count and PLINK RAW headers.
Run from ~/GS-904-final. Requires original Step17 directory on same server.
"""
import csv
import os
from pathlib import Path

OLD = Path('07_final_genotypes_crosschip_r2_0p8')
NEW = Path('08_final_genotypes_crosschip_r2_0p8_all3positive')
SPLITS = ('train', 'test')


def one_row(path):
    if not path.is_file():
        raise RuntimeError(f'Missing manifest: {path}')
    with path.open(newline='', encoding='utf-8-sig') as fp:
        items = list(csv.DictReader(fp, delimiter='\t'))
    if len(items) != 1:
        raise RuntimeError(f'Expected one manifest row: {path}')
    return items[0]


def marker_names(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing or empty RAW: {path}')
    with path.open() as fp:
        h = fp.readline().split()
    if h[:6] != ['FID','IID','PAT','MAT','SEX','PHENOTYPE']:
        raise RuntimeError(f'Bad RAW header: {path}')
    names = [x.rsplit('_', 1)[0] for x in h[6:]]
    if len(names) != len(set(names)):
        raise RuntimeError(f'Duplicate RAW SNP IDs: {path}')
    return names


def link_missing(src, dst):
    if dst.is_symlink() and not dst.exists():
        raise RuntimeError(f'Broken symlink, please investigate: {dst}')
    if dst.exists():
        return 'present'
    if not src.is_file():
        raise RuntimeError(f'Missing original Step17 source: {src}')
    # Use relative links: the two top-level directories sit side by side.
    target = os.path.relpath(src.resolve(), dst.parent.resolve())
    dst.symlink_to(target)
    if not dst.is_file():
        raise RuntimeError(f'Created link is not readable: {dst}')
    return 'linked'


def main():
    if not OLD.is_dir():
        raise RuntimeError('Step17 directory missing! Need the old 07_final_genotypes_crosschip_r2_0p8 directory on this server.')
    if not NEW.is_dir():
        raise RuntimeError('Transferred Step18 directory missing.')

    linked = 0
    nzero = nsame = nsmaller = 0
    for outer_no in range(1, 11):
        outer = f'{outer_no:02d}'
        for trait in ('BB', 'H6'):
            old = OLD / f'outer_{outer}' / trait
            new = NEW / f'outer_{outer}' / trait
            if not (new / 'complete.flag').is_file():
                raise RuntimeError(f'Missing Step18 complete.flag: {new}')
            a = one_row(old / 'final_genotype_manifest.tsv')
            b = one_row(new / 'final_genotype_manifest.tsv')
            before = int(a['Best_Extra'])
            after = int(b['Best_Extra'])
            n_chip = int(b['N_chip'])
            if before != int(b['Old_Best_Extra']) or after > before:
                raise RuntimeError(f'Inconsistent QTN counts at {new}: {before} -> {after}')
            if int(a['N_chip']) != n_chip or int(a['N_train']) != int(b['N_train']) or int(a['N_test']) != int(b['N_test']):
                raise RuntimeError(f'CHIP/sample counts disagree between Step17 and Step18 at {new}')

            for split in SPLITS:
                s = old / f'chip_{split}.raw'
                d = new / f'chip_{split}.raw'
                linked += (link_missing(s, d) == 'linked')
                if marker_names(d) != marker_names(s):
                    raise RuntimeError(f'CHIP raw marker order mismatch: {d}')
                if len(marker_names(d)) != n_chip:
                    raise RuntimeError(f'CHIP marker count mismatch: {d}')

            if after == 0:
                nzero += 1
                for split in SPLITS:
                    if (new / f'extra_{split}.raw').exists():
                        raise RuntimeError(f'QTN=0 but extra file exists in {new}')
            else:
                for split in SPLITS:
                    src = old / f'extra_{split}.raw'
                    dst = new / f'extra_{split}.raw'
                    if before == after:
                        linked += (link_missing(src, dst) == 'linked')
                    elif not dst.is_file():
                        raise RuntimeError(f'Missing reduced QTN RAW created by Step18: {dst}')
                    names = marker_names(dst)
                    if len(names) != after:
                        raise RuntimeError(f'QTN marker count mismatch {dst}: {len(names)} != {after}')
                    if before == after and names != marker_names(src):
                        raise RuntimeError(f'Unchanged-QTN marker order mismatch: {dst}')
                if before == after:
                    nsame += 1
                else:
                    nsmaller += 1
            print(f'[PASS] outer_{outer} {trait}: CHIP={n_chip}, old QTN={before}, new QTN={after}')

    if (nzero, nsame, nsmaller) != (6, 5, 9):
        raise RuntimeError(f'Unexpected category counts: zero={nzero}, same={nsame}, smaller={nsmaller}')
    print(f'RESTORED ALL 20 INPUTS: created {linked} missing symlinks; zero=6 same=5 reduced=9')


if __name__ == '__main__':
    main()
