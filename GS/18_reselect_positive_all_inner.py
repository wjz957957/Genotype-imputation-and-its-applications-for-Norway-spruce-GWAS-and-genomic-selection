#!/usr/bin/env python3
"""Re-select QTN counts: all three inner-fold Pearson deltas must be >0.

Reads Step16d's all_inner_accuracies.tsv, writes parallel selection tables
without overwriting any old optimization files. Python 3 standard library only.
"""
import argparse
import csv
import math
import os
from collections import defaultdict
from pathlib import Path

COUNTS = (0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75, 100,
          150, 200, 300, 400, 500, 600, 700)
OUTERS = tuple(f"{i:02d}" for i in range(1, 11))
INNERS = ("01", "02", "03")
TRAITS = ("BB", "H6")
EPS_TIE = 1e-8  # match the existing Step16d rule


def load_tsv(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        need = {'Outer', 'Inner', 'Trait', 'Best_Chip', 'Extra_Count',
                'Accuracy', 'RMSE', 'Delta_Accuracy_vs_0'}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'Missing columns in {path}: {sorted(missing)}')
        rows = list(reader)
    return rows


def write_tsv(path, cols, rows):
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    os.replace(temp, path)


def fmt(x):
    if isinstance(x, float):
        return format(x, '.13g')
    return x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default=(
        '09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_allfolds/'
        'optimization/all_inner_accuracies.tsv'))
    parser.add_argument('--outdir', default=(
        '09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_allfolds/'
        'optimization/three_fold_positive_rule'))
    args = parser.parse_args()

    inp = Path(args.input)
    outdir = Path(args.outdir)
    records = load_tsv(inp)
    if len(records) != 1200:
        raise ValueError(f'Expected 1200 rows, got {len(records)}')

    indexed = {}
    chip_by_group = {}
    for row in records:
        outer, inner, trait = row['Outer'].strip(), row['Inner'].strip(), row['Trait'].strip()
        if outer not in OUTERS or inner not in INNERS or trait not in TRAITS:
            raise ValueError(f'Unrecognized fold/trait: {(outer, inner, trait)}')
        count = int(row['Extra_Count'])
        if count not in COUNTS:
            raise ValueError(f'Unexpected count: {count}')
        key = (outer, trait, inner, count)
        if key in indexed:
            raise ValueError(f'Duplicate row: {key}')
        chip_key = (outer, trait)
        chip = row['Best_Chip']
        if chip_key in chip_by_group and chip_by_group[chip_key] != chip:
            raise ValueError(f'Inconsistent CHIP: {chip_key}')
        chip_by_group[chip_key] = chip
        for col in ('Accuracy', 'RMSE', 'Delta_Accuracy_vs_0'):
            v = float(row[col])
            if not math.isfinite(v):
                raise ValueError(f'Non-finite {col} at {key}')
            row[col] = v
        indexed[key] = row

    expected = {(outer, trait, inner, count) for outer in OUTERS
                for trait in TRAITS for inner in INNERS for count in COUNTS}
    if indexed.keys() != expected:
        raise ValueError(f'Incomplete rows: missing={sorted(expected - indexed.keys())[:5]}')

    choices = []
    candidate_rows = []
    changed = 0
    for outer in OUTERS:
        for trait in TRAITS:
            stats = []
            for count in COUNTS:
                xs = [indexed[(outer, trait, inner, count)] for inner in INNERS]
                bases = [indexed[(outer, trait, inner, 0)] for inner in INNERS]
                deltas = [x['Accuracy'] - b['Accuracy'] for x, b in zip(xs, bases)]
                for i, x in enumerate(xs):
                    if abs(x['Delta_Accuracy_vs_0'] - deltas[i]) > 2e-5:
                        raise ValueError(f'Delta differs from recomputed baseline: {(outer, trait, INNERS[i], count)}')
                meanacc = sum(x['Accuracy'] for x in xs) / 3
                meanrmse = sum(x['RMSE'] for x in xs) / 3
                meandelta = sum(deltas) / 3
                eligible = count == 0 or all(d > 0 for d in deltas)
                result = dict(Outer=outer, Trait=trait, Best_Chip=chip_by_group[(outer, trait)],
                              Extra_Count=count, Mean_Accuracy=meanacc, Mean_RMSE=meanrmse,
                              Mean_Delta_Pearson=meandelta, Min_Delta_Pearson=min(deltas),
                              Max_Delta_Pearson=max(deltas), N_positive=sum(d > 0 for d in deltas),
                              Delta_inner_01=deltas[0], Delta_inner_02=deltas[1],
                              Delta_inner_03=deltas[2], Eligible=int(eligible))
                stats.append(result)
                candidate_rows.append(result)

            def choose(cands):
                max_acc = max(r['Mean_Accuracy'] for r in cands)
                tied = [r for r in cands if r['Mean_Accuracy'] >= max_acc - EPS_TIE]
                return min(tied, key=lambda r: r['Extra_Count'])

            old = choose(stats)
            eligible_stats = [r for r in stats if r['Eligible']]
            new = choose(eligible_stats)
            positive_nonzero = [r for r in eligible_stats if r['Extra_Count'] > 0]
            did_change = old['Extra_Count'] != new['Extra_Count']
            changed += did_change
            choices.append(dict(Outer=outer, Trait=trait, Best_Chip=chip_by_group[(outer, trait)],
                                Original_Best_Extra=old['Extra_Count'],
                                New_Best_Extra=new['Extra_Count'],
                                Changed=int(did_change),
                                Original_Best_All3_Positive=int(old['Eligible']) if old['Extra_Count'] else 0,
                                N_Eligible_Positive_Counts=len(positive_nonzero),
                                Baseline_Mean_Pearson=stats[0]['Mean_Accuracy'],
                                New_Mean_Pearson=new['Mean_Accuracy'],
                                New_Mean_Gain=new['Mean_Delta_Pearson'],
                                New_Min_Fold_Gain=new['Min_Delta_Pearson'],
                                New_Mean_RMSE=new['Mean_RMSE']))

    outdir.mkdir(parents=True, exist_ok=True)
    choice_cols = list(choices[0])
    candidate_cols = list(candidate_rows[0])
    write_tsv(outdir / 'best_extra_all3positive.tsv', choice_cols,
              [{k: fmt(v) for k, v in r.items()} for r in choices])
    write_tsv(outdir / 'all_candidate_all3positive.tsv', candidate_cols,
              [{k: fmt(v) for k, v in r.items()} for r in candidate_rows])

    header = 'Outer Trait CHIP Old_Extra New_Extra Changed Eligible_Positive_Counts Mean_Gain Min_Fold_Gain'
    print(header)
    for r in choices:
        print(f"{r['Outer']:>5} {r['Trait']:>5} {r['Best_Chip']:>11} "
              f"{r['Original_Best_Extra']:>9} {r['New_Best_Extra']:>9} "
              f"{r['Changed']:>7} {r['N_Eligible_Positive_Counts']:>24} "
              f"{r['New_Mean_Gain']:+.6f} {r['New_Min_Fold_Gain']:+.6f}")
    print(f'\nCOMPLETE: 20 outer×trait settings. Changed {changed}/20.')
    print(f'NEW BEST: {outdir / "best_extra_all3positive.tsv"}')
    print(f'ALL 400 SETTINGS: {outdir / "all_candidate_all3positive.tsv"}')
    print('Original optimization files were NOT changed.')


if __name__ == '__main__':
    main()
