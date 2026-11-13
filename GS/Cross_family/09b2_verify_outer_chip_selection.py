#!/usr/bin/env python3
"""Audit outer-training-only chip GWAS and marker sets, based on inner Stage03.

Pre: compare Stage09a/09b1 frozen selection, phenotype, FAM, BIM, GRM safety.
Post: verify 2-trait GWAS P tables, PLINK index-SNP ranking and the exact top-N lists.
The validator never reads outer-test phenotype files.
"""
import argparse
import csv
import math
from pathlib import Path

SETS = (50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000,
        5000, 8000, 10000, 20000, 30000, 40000)
VALID_QTN = {0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75, 100,
             150, 200, 300, 400, 500, 600, 700}


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def lines(path):
    p = Path(path)
    need(p.is_file() and p.stat().st_size > 0, f'missing/empty input: {p}')
    return p.read_text().splitlines()


def ids(path, fields=1):
    items = []
    for i, line in enumerate(lines(path), 1):
        row = line.split()
        need(len(row) == fields, f'{path}:{i}: expected {fields} columns, got {len(row)}')
        if fields == 2:
            need(row[0] == row[1], f'{path}:{i}: FID != IID')
        items.append(row[0])
    need(len(items) == len(set(items)), f'duplicate IDs: {path}')
    return items


def read_tsv(path):
    with Path(path).open(newline='') as handle:
        records = list(csv.DictReader(handle, delimiter='\t'))
    need(bool(records), f'empty TSV: {path}')
    return records


def flag_has(path, expectations):
    values = set(lines(path))
    for item in expectations:
        need(item in values, f'{path}: missing completion field {item}')


def bim(path):
    snps = []
    for i, line in enumerate(lines(path), 1):
        r = line.split()
        need(len(r) == 6 and r[0] in {str(n) for n in range(1, 13)},
             f'invalid chromosome/BIM row {path}:{i}: {line[:120]}')
        snps.append(r[1])
    need(len(snps) == len(set(snps)), f'duplicate BIM SNP IDs: {path}')
    return snps


def fam(path):
    samples = []
    for i, line in enumerate(lines(path), 1):
        r = line.split()
        need(len(r) == 6 and r[0] == r[1], f'bad genotype FAM: {path}:{i}')
        samples.append(r[1])
    need(len(samples) == len(set(samples)), f'duplicate FAM IDs: {path}')
    return samples


def pheno(path):
    people = []
    for i, line in enumerate(lines(path), 1):
        r = line.split()
        need(len(r) == 3 and r[0] == r[1], f'invalid phenotype: {path}:{i}')
        try:
            value = float(r[2])
        except ValueError as e:
            raise RuntimeError(f'non-numeric phenotype: {path}:{i}') from e
        need(math.isfinite(value), f'nonfinite phenotype: {path}:{i}')
        people.append(r[1])
    need(len(people) == len(set(people)), f'duplicate phenotype IDs: {path}')
    return people


def parse_frozen(path, outer):
    selection = read_tsv(path)
    need(len(selection) == 2 and {r['Trait'] for r in selection} == {'BB', 'H6'},
         f'expected two frozen selection rows for outer {outer}')
    result = {}
    for r in selection:
        need(r['Outer'] == outer, f'frozen Outer differs from {outer}')
        best = r['Best_Chip']
        need(best.startswith('chip_') and best[5:].isdigit(), f'invalid Best_Chip: {best}')
        n = int(best[5:])
        need(n in SETS, f'unsupported best chip strategy: {best}')
        need(int(r['Best_QTN']) in VALID_QTN, f'invalid frozen QTN: {r}')
        result[r['Trait']] = (best, n)
    return result


def validate_pre(work, outer):
    c = work / 'Cross_family'
    source = c / 'results/09_outer_inputs' / f'outer_{outer}'
    genotype = c / 'results/09_outer_chip_genotypes' / f'outer_{outer}'
    base = genotype / 'chip_ALL_outer_train'
    flag_has(source / 'prepare.complete', [
        'status=PASS', f'outer={outer}', 'sample_partition_master=PASS',
        'outer_test_phenotype_read=NO', 'full_chip_consistency=PASS'])
    flag_has(genotype / 'complete.flag', [
        'status=PASS', f'outer={outer}', 'training_only=PASS',
        'outer_test_phenotype_read=NO', 'chip_ALL=43267'])
    train = ids(source / 'outer_train.keep', 2)
    test = ids(source / 'outer_test.ids')
    need(len(train) + len(test) == 904 and set(train).isdisjoint(test),
         f'invalid outer-training / outer-testing sample partition for {outer}')
    g = fam(str(base) + '.fam')
    need(set(g) == set(train) and set(g).isdisjoint(test),
         'outer-training genotype identities do not match --keep / test leakage')
    roster = ids(source / 'chip_ALL.snps')
    observed = bim(str(base) + '.bim')
    need(len(observed) == 43267 and set(roster) == set(observed),
         'outer-training genotype does not contain EXACT 43,267 chip SNP IDs')
    for t in ('BB', 'H6'):
        observed_pheno = pheno(source / f'{t}.train.pheno')
        need(observed_pheno == train, f'{t} phenotype sample/order differs from outer_train.keep')
    frozen = parse_frozen(source / 'frozen_inner_selection.tsv', outer)
    return source, base, train, observed, frozen


def validate_trait(outer, trait, out, chip_ids, frozen_best):
    clump = out / 'clumping' / trait
    gwas = out / 'gwas' / f'chip_{trait}_gwas.loco.mlma'
    header = lines(gwas)[0].split()
    need(len(header) >= 9 and header[0].lower() == 'chr' and
         header[1].lower() == 'snp' and header[8].lower() == 'p',
         f'GCTA MLMA-LOCO header is unexpected: {gwas}')
    need(len(lines(gwas)) - 1 == 43267, f'{trait} GWAS must have 43267 SNP rows')
    p_rows = lines(clump / 'gwas_p.txt')
    need(p_rows[0].split() == ['SNP', 'P'], f'invalid chip P file header {trait}')
    pvals = {}
    chip_id_set = set(chip_ids)
    for i, line in enumerate(p_rows[1:], 2):
        r = line.split()
        need(len(r) == 2 and r[0] not in pvals and r[0] in chip_id_set,
             f'invalid/duplicate/nonchip GWAS P row {trait}:{i}')
        try:
            p = float(r[1])
        except ValueError as exc:
            raise RuntimeError(f'bad P value {trait}:{i}') from exc
        need(math.isfinite(p) and 0 < p <= 1, f'nonfinite/zero/out-of-range P {trait}:{i}')
        pvals[r[0]] = p
    need(pvals, f'no valid GWAS P values: {trait}')
    clumped = lines(clump / 'clump/chip.clumped')
    need(clumped[0].split()[0:5] == ['CHR', 'F', 'SNP', 'BP', 'P'],
         f'unexpected original PLINK clump header: {trait}')
    index = {}
    for i, line in enumerate(clumped[1:], 2):
        r = line.split()
        if not r:
            continue
        need(len(r) >= 5 and r[2] not in index and r[2] in pvals,
             f'invalid/duplicate .clumped index {trait}:{i}')
        index[r[2]] = float(r[4])
    need(index, f'no PLINK clumped index SNPs for {trait}')
    ordered = sorted(index.items(), key=lambda x: (x[1], x[0]))
    ranked = lines(clump / 'ranked_chip_snps.tsv')
    ranks = []
    for i, line in enumerate(ranked, 1):
        r = line.split()
        need(len(r) == 2, f'invalid ranked chip SNP {trait}:{i}')
        ranks.append((r[0], float(r[1])))
    need(ranks == ordered, f'clumped-index rank differs from ascending P (ties SNP): {trait}')
    need(len(ranks) == len(index), f'{trait} SNP index count mismatch')
    need(all(math.isclose(p, pvals[s], rel_tol=1e-2, abs_tol=1e-300) for s, p in ranks),
         f'{trait} index SNP P does not agree with GWAS P file')
    rank_ids = [s for s, _ in ranks]
    sets = out / 'marker_sets'
    for n in SETS:
        expected = rank_ids[:n]
        got = ids(sets / f'{trait}_chip_{n}.snps')
        need(got == expected,
             f'{trait} chip_{n} not exactly the first {len(expected)} outer-training index SNPs')
    need(ids(sets / f'{trait}_chip_LD_ALL.snps') == rank_ids,
         f'{trait} chip_LD_ALL SNP order differs from clump index rank')
    need(ids(sets / f'{trait}_chip_ALL.snps') == chip_ids,
         f'{trait} full chip marker IDs must equal outer-training BIM order')
    strategy, chosen_n = frozen_best
    need(len(rank_ids) >= chosen_n,
         f'frozen strategy {outer}/{trait}/{strategy} needs {chosen_n}, but only {len(rank_ids)} index SNPs')
    actual = ids(sets / f'{trait}_{strategy}.snps')
    need(actual == rank_ids[:chosen_n],
         f'outer GWAS-derived best chip SNP IDs do not match frozen strategy: {trait}')
    print(f'PASS outer_{outer}/{trait}: chip={strategy} ({chosen_n}), '
          f'index_SNPs={len(rank_ids)}, GWAS_valid_P={len(pvals)}')
    return len(index)


def validate_post(work, outer, base, chip_ids, train, frozen, no_flag_check):
    out = work / 'Cross_family/results/09_outer_chip_gwas_clumping' / f'outer_{outer}'
    if not no_flag_check:
        flag_has(out / 'complete.flag', ['status=PASS', 'version=1',
                 f'outer={outer}', 'training_only=PASS',
                 'outer_test_phenotype_read=NO',
                 'chip_GWAS=GCTA_MLMA_LOCO',
                 'chip_clumping_r2=0.2', 'chip_clumping_kb=50',
                 'selected_outer_chip_SNPs=PASS'])
    grm = out / 'grm/chip_outer_train_grm'
    for ext in ('.grm.bin', '.grm.N.bin', '.grm.id'):
        need(Path(str(grm) + ext).is_file() and Path(str(grm) + ext).stat().st_size > 0,
             f'missing GRM output: {grm}{ext}')
    g = ids(str(grm) + '.grm.id', 2)
    need(set(g) == set(train), f'GRM sample IDs must match outer training {outer}')
    qc = read_tsv(out / 'gwas_qc.tsv')
    need(len(qc) == 2 and {r['Trait'] for r in qc} == {'BB','H6'},
         'GWAS QC must have BB and H6 rows')
    for r in qc:
        need(int(r['Total_SNPs']) == 43267 and
             int(r['Valid_P']) + int(r['Nonfinite_P']) == 43267 and
             int(r['Valid_P']) > 0, f'incorrect GWAS QC for {r["Trait"]}')
    indexes = {}
    for trait in ('BB', 'H6'):
        indexes[trait] = validate_trait(outer, trait, out, chip_ids, frozen[trait])
    summary = read_tsv(out / 'marker_set_counts.tsv')
    need(len(summary) == 36, f'expected 36 marker count rows (18 per trait), got {len(summary)}')
    for trait in ('BB', 'H6'):
        subset = [r for r in summary if r['Trait'] == trait]
        need(len(subset) == 18, f'expected 18 counts for {trait}')
        actual = {r['Marker_Set']: r for r in subset}
        need(len(actual) == 18, f'duplicate marker set count row {trait}')
        for n in SETS:
            r = actual[f'chip_{n}']
            need(r['Requested_N'] == str(n) and int(r['Actual_N']) == min(n, indexes[trait]),
                 f'marker count table disagrees with {trait}/chip_{n}')
        r = actual['chip_LD_ALL']
        need(r['Requested_N'] == 'ALL_CLUMPED' and int(r['Actual_N']) == indexes[trait],
             'incorrect chip_LD_ALL summary')
        r = actual['chip_ALL']
        need(int(r['Actual_N']) == 43267 and int(r['Requested_N']) == 43267,
             'incorrect chip_ALL summary')
    print(f'OVERALL PASS outer_{outer}: n_train={len(train)}, chip_ALL=43267, '
          f'BB={frozen["BB"][0]}, H6={frozen["H6"][0]}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--work', required=True)
    p.add_argument('--outer', required=True)
    p.add_argument('--phase', choices=('pre', 'post'), required=True)
    p.add_argument('--no-flag-check', action='store_true')
    a = p.parse_args()
    outer = f'{int(a.outer):02d}'
    need(outer in {'01','02','03','04','05'}, f'invalid outer: {outer}')
    work = Path(a.work)
    source, base, train, chip_ids, frozen = validate_pre(work, outer)
    if a.phase == 'post':
        validate_post(work, outer, base, chip_ids, train, frozen, a.no_flag_check)
    else:
        print(f'PASS 09b2 pre outer_{outer}: n_train={len(train)}, '
              f'chip_ALL={len(chip_ids)}, no outer-test phenotype input')


if __name__ == '__main__':
    main()
