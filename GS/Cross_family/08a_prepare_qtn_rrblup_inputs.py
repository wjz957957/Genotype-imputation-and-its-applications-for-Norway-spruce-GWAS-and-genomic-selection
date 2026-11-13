#!/usr/bin/env python3
"""Stage 08a: validate inner-fold data and verify PLINK RAW extraction.

Never calculate LD or re-rank QTN here. Marker selection is frozen by stages 06/07.
"""
import argparse
import csv
import math
import os
from pathlib import Path

COUNTS = (0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75,
          100, 150, 200, 300, 400, 500, 600, 700)


def ensure(cond, message):
    if not cond:
        raise ValueError(message)


def lines(p, allow_empty=False):
    p = Path(p)
    ensure(p.is_file(), f"missing file: {p}")
    vals = [s.strip() for s in p.read_text().splitlines()]
    if not allow_empty:
        ensure(bool(vals), f"empty file: {p}")
    ensure(all(s and len(s.split()) == 1 for s in vals), f"invalid one-column list: {p}")
    ensure(len(vals) == len(set(vals)), f"duplicate IDs in {p}")
    return vals


def pheno_ids(p):
    p = Path(p)
    ensure(p.is_file() and p.stat().st_size > 0, f"phenotype missing: {p}")
    ids = []
    with p.open() as f:
        for j, row in enumerate(f, 1):
            a = row.split()
            ensure(len(a) == 3, f"phenotype has not 3 columns {p}:{j}")
            ensure(a[0] == a[1], f"FID != IID at {p}:{j}")
            try:
                v = float(a[2])
            except ValueError:
                raise ValueError(f"non-numeric phenotype {p}:{j}")
            ensure(math.isfinite(v), f"nonfinite phenotype {p}:{j}")
            ids.append(a[1])
    ensure(len(ids) > 0 and len(ids) == len(set(ids)), f"invalid phenotype IDs: {p}")
    return ids


def read_config(p, outer, trait):
    with open(p, newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    ensure(rows and {'Outer', 'Trait', 'Best_Chip'}.issubset(rows[0]), f"bad config headers {p}")
    vals = [r['Best_Chip'] for r in rows if r['Outer'] == outer and r['Trait'] == trait]
    ensure(len(vals) == 1, f"missing/duplicate chip selection for {outer}/{trait}")
    ensure(vals[0].startswith('chip_') and vals[0][5:].isdigit(), "invalid chip setting")
    return vals[0]


def read_bim_ids(p):
    s = set()
    with open(p) as f:
        for j, line in enumerate(f, 1):
            c = line.split()
            ensure(len(c) == 6, f"bad BIM line {p}:{j}")
            ensure(c[1] not in s, f"duplicated BIM SNP ID {c[1]}")
            s.add(c[1])
    return s


def atomic_lines(path, vals):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(''.join(s + '\n' for s in vals))
    os.replace(tmp, path)


def flag(p):
    return dict(x.split('=', 1) for x in Path(p).read_text().splitlines() if '=' in x)


def parse(args):
    w = Path(args.work)
    o, i, t = args.outer, args.inner, args.trait
    fold = Path(f'outer_{o}') / f'inner_{i}'
    out = Path(args.out)
    source = w / 'Cross_family'
    master = w / '08_imputed_GWAS/00_reference/all_904_genotypes'
    selection = source / 'config/best_chip_by_outer_trait.tsv'
    size = read_config(selection, o, t)
    chip_p = source / 'results/03_inner_chip_gwas_clumping' / fold / 'marker_sets' / f'{t}_{size}.snps'
    extras = source / 'results/07_inner_crosschip_ld_filtering' / fold / t
    train_p = source / 'results/02_inner_phenotypes' / fold / f'{t}.train.pheno'
    valid_p = source / 'results/02_inner_phenotypes' / fold / f'{t}.valid.pheno'
    cv = source / 'cv_samples_v1' / f'outer_{o}'
    train_ids_p = cv / f'inner_{i}/inner_train.samples'
    valid_ids_p = cv / f'inner_{i}/inner_valid.samples'
    test_ids_p = cv / 'outer_test.samples'
    for suffix in ('.bed', '.bim', '.fam'):
        p = Path(str(master) + suffix)
        ensure(p.is_file() and p.stat().st_size > 0, f"missing master genotype {p}")

    f = flag(extras / 'complete.flag')
    for key, value in [('status', 'PASS'), ('outer', o), ('inner', i),
                       ('trait', t), ('chip_strategy', size),
                       ('all_700_available', 'YES'), ('training_only_ld', 'PASS')]:
        ensure(f.get(key) == value, f"07a flag mismatch {key}={f.get(key)} expected {value}")
    chip = lines(chip_p)
    ensure(len(chip) == int(size[5:]), f"wrong chip count {chip_p}")
    extra = lines(extras / 'marker_sets/extra_700.snps')
    ensure(len(extra) == 700, f"wrong extra 700 count {extras}")
    ensure(not (set(chip) & set(extra)), 'chip and QTN overlap')
    ensure((extras / 'marker_sets/extra_0.snps').is_file() and
           (extras / 'marker_sets/extra_0.snps').stat().st_size == 0,
           'extra_0 list should exist and be empty')
    for k in COUNTS[1:]:
        x = lines(extras / 'marker_sets' / f'extra_{k}.snps')
        ensure(x == extra[:k], f"QTN nested-set order mismatch extra_{k}")
    ranked = extras / 'filtered_ranked_extra_snps.tsv'
    with ranked.open(newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        ensure(header == ['SNP', 'P'], 'incorrect 07a ranked header')
        for k, row in enumerate(reader):
            if k == 700:
                break
            ensure(row[0] == extra[k], f"candidate {k} not first top700 from ranked SNPs")
    train, valid, test = lines(train_ids_p), lines(valid_ids_p), lines(test_ids_p)
    ensure(not(set(train)&set(valid) or set(train)&set(test) or set(valid)&set(test)),
           'train, validation or outer-test sample overlap')
    ensure(set(train) == set(pheno_ids(train_p)), 'training phenotype IDs mismatch CV split')
    ensure(set(valid) == set(pheno_ids(valid_p)), 'validation phenotype IDs mismatch CV split')
    ensure(100 <= len(train) <= 904 and 50 <= len(valid) <= 904,
           'unexpected train or validation sample count')
    in_bim = read_bim_ids(str(master)+'.bim')
    missing = set(chip+extra) - in_bim
    ensure(not missing, f"chip/QTN absent from master BIM: {sorted(missing)[:5]}")
    fam = {}
    with open(str(master)+'.fam') as fh:
        for line in fh:
            ff = line.split()
            ensure(len(ff) >= 2 and ff[0] == ff[1], 'FID/IID mismatch in master FAM')
            ensure(ff[1] not in fam, 'duplicate master FAM sample')
            fam[ff[1]] = 1
    ensure(set(train+valid+test) <= set(fam), 'CV sample absent from master FAM')
    ensure(out.exists(), f"output path missing: {out}")
    atomic_lines(out/'selected_chip.snps', chip)
    atomic_lines(out/'extra_700.snps', extra)
    atomic_lines(out/'union.snps', chip + extra)
    atomic_lines(out/'train.keep', [f'{s}\t{s}' for s in train])
    atomic_lines(out/'valid.keep', [f'{s}\t{s}' for s in valid])
    atomic_lines(out/'metadata.tsv', ['field\tvalue', f'outer\t{o}', f'inner\t{i}',
                                       f'trait\t{t}', f'chip_strategy\t{size}',
                                       f'n_chip\t{len(chip)}', 'n_extra\t700',
                                       f'n_train\t{len(train)}', f'n_valid\t{len(valid)}',
                                       'genotype_source\tall_904_genotypes'])
    print(f'PREP PASS outer={o} inner={i} {t}, chip={len(chip)}, extra=700, train={len(train)}, valid={len(valid)}')


def raw_info(path, expected):
    p = Path(path)
    ensure(p.is_file() and p.stat().st_size > 0, f"missing RAW file {p}")
    with p.open() as f:
        head = f.readline().split()
        ensure(head[:6] == ['FID','IID','PAT','MAT','SEX','PHENOTYPE'], f"invalid RAW header {p}")
        markers = head[6:]
        ensure(len(markers) >= 1, f"RAW genotype has no SNPs {p}")
        snps = [s.rsplit('_',1)[0] for s in markers]
        for name in markers:
            ensure('_' in name, f"RAW marker header lacks counted allele: {name}")
            sid, allele = name.rsplit('_', 1)
            parts = sid.split(':')
            ensure(len(parts) == 4 and allele in parts[2:4],
                   f"RAW counted allele incompatible with SNP ID: {name}")
        ensure(len(snps) == len(set(snps)), f"duplicate RAW markers in {p}")
        ids = []
        for j, line in enumerate(f, 2):
            row = line.split()
            ensure(len(row) == len(head), f"column count mismatch {p}:{j}")
            ensure(row[0] == row[1], f"RAW FID/IID mismatch {p}:{j}")
            ids.append(row[1])
    ensure(len(ids) == len(expected) and set(ids) == set(expected) and len(ids) == len(set(ids)),
           f"RAW sample mismatch in {p}")
    return markers, snps


def verify(args):
    out = Path(args.out)
    chip = lines(out/'selected_chip.snps')
    extra = lines(out/'extra_700.snps')
    train = [s.split()[1] for s in (out/'train.keep').read_text().splitlines()]
    valid = [s.split()[1] for s in (out/'valid.keep').read_text().splitlines()]
    a, aid = raw_info(out/'train.raw', train)
    b, bid = raw_info(out/'valid.raw', valid)
    ensure(a == b, 'training/validation RAW SNP headers or counted alleles differ')
    ensure(set(aid) == set(chip + extra) and len(aid) == len(chip)+700,
           'missing or extra SNPs in PLINK RAW')
    atomic_lines(out/'train_snp_headers.txt', a)
    atomic_lines(out/'valid_snp_headers.txt', b)
    atomic_lines(out/'complete.flag', ['version=1', f'outer={args.outer}', f'inner={args.inner}',
                                       f'trait={args.trait}', f'n_chip={len(chip)}',
                                       'n_extra=700', f'n_train={len(train)}',
                                       f'n_valid={len(valid)}', 'sample_ID_match=PASS',
                                       'train_valid_snp_header_match=PASS',
                                       'outer_test_excluded=PASS',
                                       'training_only_QTN_selection=PASS', 'status=PASS'])
    print(f'VERIFY PASS: {args.outer}/{args.inner}/{args.trait}: {len(train)} train, {len(valid)} valid, {len(aid)} SNPs')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('prepare', 'verify'))
    p.add_argument('--work', required=True)
    p.add_argument('--outer', required=True)
    p.add_argument('--inner', required=True)
    p.add_argument('--trait', required=True, choices=('BB','H6'))
    p.add_argument('--out', required=True)
    a = p.parse_args()
    ensure(a.outer in [f'{x:02d}' for x in range(1,6)], 'invalid outer')
    ensure(a.inner in ['01','02','03'], 'invalid inner')
    {'prepare': parse, 'verify': verify}[a.action](a)

if __name__ == '__main__':
    main()
