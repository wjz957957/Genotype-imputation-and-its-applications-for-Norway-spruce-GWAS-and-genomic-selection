#!/usr/bin/env python3
"""Preflight all 10 outer-fold/trait cross-chip LD filtering tasks."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

VALID_CHIPS = {50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000,
               5000, 8000, 10000, 20000, 30000, 40000}
VALID_QTN = {0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75, 100,
             150, 200, 300, 400, 500, 600, 700}
TRAITS = ("BB", "H6")


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError(msg)


def required(path: Path) -> Path:
    need(path.is_file() and path.stat().st_size > 0, f"missing/empty input: {path}")
    return path


def read_flag(path: Path) -> dict[str, str]:
    out = {}
    for line in required(path).read_text().splitlines():
        k, sep, v = line.partition("=")
        need(bool(sep) and k not in out, f"malformed/duplicate flag line in {path}: {line}")
        out[k] = v
    return out


def check_flag(path: Path, expected: dict[str, str]) -> None:
    observed = read_flag(path)
    for k, v in expected.items():
        need(observed.get(k) == v, f"{path}: expected {k}={v}, found {observed.get(k)!r}")


def read_ids(path: Path) -> list[str]:
    vals = []
    for n, line in enumerate(required(path).read_text().splitlines(), 1):
        f = line.split()
        need(len(f) == 1, f"expected one ID per line: {path}:{n}")
        vals.append(f[0])
    need(bool(vals) and len(vals) == len(set(vals)), f"empty/duplicate IDs: {path}")
    return vals


def read_keep(path: Path) -> list[str]:
    vals = []
    for n, line in enumerate(required(path).read_text().splitlines(), 1):
        f = line.split()
        need(len(f) == 2 and f[0] == f[1], f"invalid FID/IID: {path}:{n}")
        vals.append(f[1])
    need(bool(vals) and len(vals) == len(set(vals)), f"empty/duplicate train IDs: {path}")
    return vals


def read_ranked(path: Path) -> list[tuple[str, float]]:
    rows = []
    seen = set()
    prev = -1.0
    with required(path).open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        need(header == ["SNP", "P"], f"invalid ranking header in {path}: {header}")
        for n, r in enumerate(reader, 2):
            need(len(r) == 2 and bool(r[0]) and r[0] not in seen,
                 f"malformed/duplicate ranking row {path}:{n}")
            try:
                p = float(r[1])
            except ValueError as exc:
                raise RuntimeError(f"non-numeric P at {path}:{n}") from exc
            need(math.isfinite(p) and 0 < p <= 1 and p + 1e-15 >= prev,
                 f"invalid or unsorted P at {path}:{n}")
            rows.append((r[0], p))
            seen.add(r[0])
            prev = p
    need(bool(rows), f"no clumping index SNPs: {path}")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    work = Path(a.work)
    cf = work / "Cross_family"
    master = work / "08_imputed_GWAS/00_reference/all_904_genotypes"
    for ext in ("bed", "bim", "fam"):
        required(Path(str(master) + f".{ext}"))
    merged_root = cf / "results/09_outer_imputed_clumping_merged"
    check_flag(merged_root / "complete.flag", {
        "status": "PASS", "outer_folds": "5/5", "trait_combinations": "10/10",
        "training_only": "PASS", "outer_test_phenotype_read": "NO"})
    with required(merged_root / "merge_summary.tsv").open(newline="") as fh:
        global_rows = list(csv.DictReader(fh, delimiter="\t"))
    need(len(global_rows) == 10, f"expected 10 rows in global merge_summary.tsv, got {len(global_rows)}")
    merged_summary = {(r["outer"], r["trait"]): int(r["n_candidate_snps"]) for r in global_rows}

    passed = 0
    for oi in range(1, 6):
        outer = f"{oi:02d}"
        inp = cf / f"results/09_outer_inputs/outer_{outer}"
        imp = cf / f"results/09_outer_imputed_gwas_inputs/outer_{outer}"
        chip_fold = cf / f"results/09_outer_chip_gwas_clumping/outer_{outer}"
        merged_fold = merged_root / f"outer_{outer}"
        check_flag(inp / "prepare.complete", {
            "status": "PASS", "outer": outer, "full_chip_consistency": "PASS",
            "outer_test_phenotype_read": "NO"})
        check_flag(imp / "complete.flag", {
            "status": "PASS", "outer": outer, "training_only": "PASS",
            "outer_test_phenotype_read": "NO"})
        check_flag(chip_fold / "complete.flag", {
            "status": "PASS", "outer": outer, "training_only": "PASS",
            "outer_test_phenotype_read": "NO", "selected_outer_chip_SNPs": "PASS"})
        check_flag(merged_fold / "complete.flag", {
            "status": "PASS", "outer": outer, "chromosomes": "12/12",
            "training_only": "PASS", "outer_test_phenotype_read": "NO"})

        train = read_keep(inp / "outer_train.keep")
        test = read_ids(inp / "outer_test.ids")
        need(len(train) + len(test) == 904 and set(train).isdisjoint(test),
             f"outer train/test partition invalid outer_{outer}")
        with required(inp / "frozen_inner_selection.tsv").open(newline="") as fh:
            frozen = list(csv.DictReader(fh, delimiter="\t"))
        need(len(frozen) == 2 and {r.get("Trait") for r in frozen} == set(TRAITS),
             f"frozen selection must have BB/H6 rows outer_{outer}")

        for row in frozen:
            trait = row["Trait"]
            best_chip = row["Best_Chip"]
            target = int(row["Best_QTN"])
            need(best_chip.startswith("chip_") and best_chip[5:].isdigit(),
                 f"invalid chip strategy {outer}/{trait}: {best_chip}")
            chip_n = int(best_chip[5:])
            need(chip_n in VALID_CHIPS and target in VALID_QTN,
                 f"unsupported frozen chip/QTN choice {outer}/{trait}: {best_chip}/{target}")
            chip_path = chip_fold / f"marker_sets/{trait}_{best_chip}.snps"
            chip_ids = read_ids(chip_path)
            need(len(chip_ids) == chip_n, f"chip list count mismatch {chip_path}")
            ranked_path = merged_fold / f"{trait}_index_ranked.tsv"
            rows = read_ranked(ranked_path)
            need(len(rows) == merged_summary[(outer, trait)],
                 f"merged candidate count differs from global summary outer_{outer}/{trait}")
            chip_set = set(chip_ids)
            n_overlap = sum(snp in chip_set for snp, _ in rows)
            n_nonchip = len(rows) - n_overlap
            need(target == 0 or n_nonchip >= target,
                 f"too few non-chip candidates before LD filtering outer_{outer}/{trait}: {n_nonchip} < {target}")
            print(f"PASS outer_{outer}/{trait}: chip={best_chip} ({chip_n}), frozen_QTN={target}, "
                  f"index_SNPs={len(rows)}, exact_overlaps_to_remove={n_overlap}, nonchip={n_nonchip}, "
                  f"training={len(train)} test_excluded={len(test)}")
            passed += 1

    print(f"OVERALL 09g PREFLIGHT PASS={passed}/10 FAIL={10-passed}/10")
    need(passed == 10, "not all 10 outer fold/trait tasks passed preflight")


if __name__ == "__main__":
    main()
