#!/usr/bin/env python3
"""Merge outer-training per-chromosome PLINK index SNPs; no outer-test input.

Important: matching the successful inner 06d manual merge layout, output is
  09_outer_imputed_clumping_merged/outer_XX/{BB,H6}_index_ranked.tsv
and one global merge_summary.tsv with outer,trait,n_candidate_snps.
"""
from __future__ import annotations
import argparse
import csv
import math
import os
from pathlib import Path

TRAITS = ("BB", "H6")


def require(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty file: {path}")
    return path


def read_flag(path: Path) -> dict[str, str]:
    result = {}
    for line in require(path).read_text().splitlines():
        key, sep, value = line.partition("=")
        if sep and key in result:
            raise RuntimeError(f"Duplicate completion key in {path}: {key}")
        if sep:
            result[key] = value
    return result


def atomic_write(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f".tmp.{os.getpid()}")
    temp.write_text("\n".join(lines) + "\n")
    os.replace(temp, path)


def parse_ranked(path: Path) -> list[tuple[str, float]]:
    with require(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["SNP", "P"]:
            raise RuntimeError(f"Invalid SNP/P header {path}: {header}")
        rows = []
        seen = set()
        for n, row in enumerate(reader, start=2):
            if len(row) != 2 or not row[0] or row[0] in seen:
                raise RuntimeError(f"Malformed or duplicate SNP at {path}:{n}")
            try:
                p = float(row[1])
            except ValueError as exc:
                raise RuntimeError(f"Non-numeric P at {path}:{n}") from exc
            if not math.isfinite(p) or not 0 < p <= 1:
                raise RuntimeError(f"P out of range at {path}:{n}: {p}")
            seen.add(row[0])
            rows.append((row[0], p))
        if not rows:
            raise RuntimeError(f"Zero index SNPs in {path}")
        return rows


def merge_outer(work: Path, outer: str) -> None:
    inp = work / "Cross_family/results/09_outer_imputed_clumping_by_chr" / f"outer_{outer}"
    output = work / "Cross_family/results/09_outer_imputed_clumping_merged" / f"outer_{outer}"
    # Invalidate a previous completion marker before validating or rewriting any outputs.
    (output / "complete.flag").unlink(missing_ok=True)
    pending: dict[str, list[tuple[str, float]]] = {t: [] for t in TRAITS}
    per_chr_counts: list[tuple[int, str, int]] = []

    for chrnum in range(1, 13):
        chrom = f"{chrnum:02d}"
        chrom_dir = inp / f"chr_{chrom}"
        flag = read_flag(chrom_dir / "complete.flag")
        required = {"status": "PASS", "outer": outer, "chr": str(chrnum),
                    "r2": "0.2", "window_kb": "1000", "training_only": "PASS",
                    "outer_test_phenotype_read": "NO", "traits": "BB,H6"}
        for key, value in required.items():
            if flag.get(key) != value:
                raise RuntimeError(f"Invalid 09e completion flag: {chrom_dir}, {key}={flag.get(key)!r}")
        for trait in TRAITS:
            rank_path = chrom_dir / f"{trait}_chr{chrom}.index_ranked.tsv"
            rows = parse_ranked(rank_path)
            pending[trait].extend(rows)
            per_chr_counts.append((chrnum, trait, len(rows)))

    counts = {}
    files = {}
    # Validate all data *before* marking this outer fold complete.
    for trait in TRAITS:
        rows = pending[trait]
        names = [s for s, _ in rows]
        if len(set(names)) != len(names):
            raise RuntimeError(f"Duplicate clump SNPs between chromosomes for outer_{outer}/{trait}")
        if not 700 <= len(rows) <= 3200000:
            raise RuntimeError(f"Unexpected genome-wide index SNP count: outer_{outer}/{trait}={len(rows)}")
        rows.sort(key=lambda pair: (pair[1], pair[0]))
        files[trait] = rows
        counts[trait] = len(rows)

    for trait in TRAITS:
        rows = files[trait]
        atomic_write(output / f"{trait}_index_ranked.tsv",
                     ["SNP\tP"] + [f"{snp}\t{p:.15g}" for snp, p in rows])
    atomic_write(output / "per_chr_counts.tsv", ["chr\ttrait\tn_index_snps"] +
                 [f"{chrnum}\t{trait}\t{n}" for chrnum, trait, n in per_chr_counts])
    atomic_write(output / "merge_summary.tsv", ["outer\ttrait\tn_candidate_snps"] +
                 [f"{outer}\t{trait}\t{counts[trait]}" for trait in TRAITS])
    atomic_write(output / "complete.flag", ["status=PASS", "version=1", f"outer={outer}",
                 "chromosomes=12/12", "traits=BB,H6", "training_only=PASS",
                 "outer_test_phenotype_read=NO", "clumping_r2=0.2",
                 "clumping_kb=1000", "ranking=GWAS_P_ascending"])
    print(f"PASS outer_{outer}: BB={counts['BB']} H6={counts['H6']} genome-wide clump index SNPs")


def verify_all(work: Path) -> None:
    merged = work / "Cross_family/results/09_outer_imputed_clumping_merged"
    lines = ["outer\ttrait\tn_candidate_snps"]
    for num in range(1, 6):
        outer = f"{num:02d}"
        root = merged / f"outer_{outer}"
        flag = read_flag(root / "complete.flag")
        if any(flag.get(k) != v for k, v in {"status": "PASS", "outer": outer,
                                              "chromosomes": "12/12", "training_only": "PASS"}.items()):
            raise RuntimeError(f"Invalid 09f fold flag outer_{outer}")
        with require(root / "merge_summary.tsv").open(newline="") as fh:
            entries = list(csv.DictReader(fh, delimiter="\t"))
        if len(entries) != 2 or {entry["trait"] for entry in entries} != set(TRAITS):
            raise RuntimeError(f"Invalid trait summary outer_{outer}")
        for entry in entries:
            trait = entry["trait"]
            n = int(entry["n_candidate_snps"])
            if entry["outer"] != outer or n < 700:
                raise RuntimeError(f"Invalid count for outer_{outer}/{trait}")
            rows = parse_ranked(root / f"{trait}_index_ranked.tsv")
            if len(rows) != n or rows != sorted(rows, key=lambda pair: (pair[1], pair[0])):
                raise RuntimeError(f"Merged SNP list count/order mismatch outer_{outer}/{trait}")
            lines.append(f"{outer}\t{trait}\t{n}")
            print(f"PASS outer_{outer}/{trait}: {n} merged index SNPs")
    atomic_write(merged / "merge_summary.tsv", lines)
    atomic_write(merged / "complete.flag", ["status=PASS", "version=1", "outer_folds=5/5",
              "trait_combinations=10/10", "training_only=PASS", "outer_test_phenotype_read=NO"])
    print("OVERALL 09f PASS=5/5; 10 trait candidate ranked lists ready")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True)
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--outer", choices=[f"{i:02d}" for i in range(1, 6)])
    g.add_argument("--verify-all", action="store_true")
    args = parser.parse_args()
    if args.verify_all:
        verify_all(Path(args.work))
    else:
        merge_outer(Path(args.work), args.outer)


if __name__ == "__main__":
    main()
