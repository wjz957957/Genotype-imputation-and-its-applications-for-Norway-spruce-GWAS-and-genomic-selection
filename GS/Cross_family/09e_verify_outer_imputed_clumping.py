#!/usr/bin/env python3
"""Read-only audit of all 60 outer-training imputed clumping tasks."""
import argparse
import csv
import math
from pathlib import Path


def readflag(path):
    if not path.is_file():
        raise RuntimeError(f"Missing {path}")
    lines = path.read_text().splitlines()
    flag = {}
    for line in lines:
        k, sep, v = line.partition('=')
        if sep:
            if k in flag:
                raise RuntimeError(f"Duplicate flag key {k}: {path}")
            flag[k] = v
    return flag


def rows(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty {path}")
    with path.open(newline='') as fd:
        table = csv.reader(fd, delimiter='\t')
        h = next(table, [])
        if h != ['SNP', 'P']:
            raise RuntimeError(f"Invalid header {path}: {h}")
        seen = set()
        previous = None
        count = 0
        for n, row in enumerate(table, 2):
            if len(row)!=2 or row[0] in seen:
                raise RuntimeError(f"Malformed/duplicate SNP {path}:{n}")
            try:
                p = float(row[1])
            except ValueError:
                raise RuntimeError(f"Invalid p {path}:{n}")
            if not math.isfinite(p) or p <= 0 or p > 1:
                raise RuntimeError(f"Invalid P value {path}:{n}")
            order = (p, row[0])
            if previous is not None and order < previous:
                raise RuntimeError(f"Not P-ranked {path}:{n}")
            previous = order
            seen.add(row[0])
            count += 1
    if count < 1:
        raise RuntimeError(f"No index SNPs {path}")
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--work', required=True)
    args = parser.parse_args()
    work = Path(args.work)
    root = work/'Cross_family/results/09_outer_imputed_clumping_by_chr'
    qcroot = work/'Cross_family/results/09_outer_imputed_gwas'
    total_pairs = total_snps = 0
    for j in range(1, 6):
        o = f'{j:02d}'
        nidx = 0
        for c in range(1, 13):
            pad = f'{c:02d}'
            d = root/f'outer_{o}'/f'chr_{pad}'
            flag = readflag(d/'complete.flag')
            for k,v in {'status':'PASS','version':'1','outer':o,'chr':str(c),
                        'r2':'0.2','window_kb':'1000','training_only':'PASS',
                        'outer_test_phenotype_read':'NO','traits':'BB,H6'}.items():
                if flag.get(k)!=v:
                    raise RuntimeError(f'Bad 09e flag {d}: expected {k}={v}')
            with (d/'clumping_summary.tsv').open(newline='') as fd:
                summary = list(csv.DictReader(fd, delimiter='\t'))
            if len(summary)!=2 or {s.get('Trait') for s in summary}!={'BB','H6'}:
                raise RuntimeError(f'Wrong trait count in summary: {d}')
            with (qcroot/f'outer_{o}'/f'chr_{pad}'/'gwas_qc.tsv').open(newline='') as fd:
                gwasqc = {s['Trait']:s for s in csv.DictReader(fd, delimiter='\t')}
            for record in summary:
                t = record['Trait']
                n = rows(d/f'{t}_chr{pad}.index_ranked.tsv')
                n_txt = (d/f'{t}_chr{pad}.index_snps.txt').read_text().splitlines()
                with (d/f'{t}_chr{pad}.index_ranked.tsv').open() as fh:
                    ranked_ids = [line.split('\t')[0] for line in fh.read().splitlines()[1:]]
                if n != len(n_txt) or n_txt != ranked_ids or n != len(set(n_txt)):
                    raise RuntimeError(f'Index SNP text mismatch: {d}/{t}')
                if not (record['Outer']==o and record['Chr']==str(c) and record['R2']=='0.2' and
                        record['Window_kb']=='1000' and record['Clumped_Index_SNPs']==str(n) and
                        record['GWAS_SNPs']==gwasqc[t]['Expected_SNPs']):
                    raise RuntimeError(f'Summary/QC mismatch: {d}/{t}')
                total_pairs+=1
                nidx+=n
        total_snps+=nidx
        print(f'PASS outer_{o}: 12 chromosomes, 24 BB/H6 index SNP lists, combined index rows={nidx}')
    print(f'OVERALL 09e PASS=5/5: 60 chromosome tasks, {total_pairs} trait/chromosome lists, {total_snps} clump index rows')

if __name__=='__main__':
    main()
