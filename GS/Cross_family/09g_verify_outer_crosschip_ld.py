#!/usr/bin/env python3
"""Verify outer-training-only cross-chip LD filtering and frozen QTN counts."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

TRAITS = ("BB", "H6")


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def file(path: Path, allow_empty: bool = False) -> Path:
    require(path.is_file(), f"missing file: {path}")
    if not allow_empty:
        require(path.stat().st_size > 0, f"empty file: {path}")
    return path


def lines(path: Path, allow_empty: bool = False) -> list[str]:
    p = file(path, allow_empty=allow_empty)
    return p.read_text().splitlines()


def flag(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in lines(path):
        k, sep, v = line.partition("=")
        require(bool(sep) and k not in out, f"malformed/duplicate flag key in {path}: {line}")
        out[k] = v
    return out


def read_ids(path: Path, allow_empty: bool = False) -> list[str]:
    vals: list[str] = []
    for n, line in enumerate(lines(path, allow_empty=allow_empty), 1):
        f = line.split()
        require(len(f) == 1, f"expected one ID per line: {path}:{n}")
        vals.append(f[0])
    require(allow_empty or bool(vals), f"no IDs in {path}")
    require(len(vals) == len(set(vals)), f"duplicate IDs in {path}")
    return vals


def read_train_keep(path: Path) -> list[str]:
    vals: list[str] = []
    for n, line in enumerate(lines(path), 1):
        f = line.split()
        require(len(f) == 2 and f[0] == f[1], f"invalid FID/IID in {path}:{n}")
        vals.append(f[1])
    require(len(vals) == len(set(vals)), f"duplicate train IDs in {path}")
    return vals


def read_ranked(path: Path, allow_empty: bool = False) -> list[tuple[str, str]]:
    with file(path, allow_empty=allow_empty).open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        require(header == ["SNP", "P"], f"invalid SNP/P header in {path}: {header}")
        out: list[tuple[str, str]] = []
        seen: set[str] = set()
        prev = -1.0
        for n, row in enumerate(reader, 2):
            require(len(row) == 2 and bool(row[0]), f"bad ranked row {path}:{n}")
            require(row[0] not in seen, f"duplicate ranked SNP {row[0]} in {path}")
            try:
                p = float(row[1])
            except ValueError as exc:
                raise RuntimeError(f"non-numeric P at {path}:{n}") from exc
            require(math.isfinite(p) and 0 <= p <= 1, f"P out of range at {path}:{n}")
            require(p + 1e-15 >= prev, f"P values not ascending at {path}:{n}")
            seen.add(row[0])
            out.append((row[0], row[1]))
            prev = p
    require(allow_empty or bool(out), f"no ranked SNP rows in {path}")
    return out


def read_pairs(path: Path) -> list[dict[str, str]]:
    with file(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kv_table(path: Path) -> dict[str, str]:
    rows = lines(path)
    require(rows[0].split("\t") == ["Item", "Value"], f"invalid Item/Value header: {path}")
    out: dict[str, str] = {}
    for n, line in enumerate(rows[1:], 2):
        f = line.split("\t")
        if len(f) != 2:
            continue
        require(f[0] not in out, f"duplicate summary key {f[0]} in {path}:{n}")
        out[f[0]] = f[1]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    work = Path(a.work)
    cf = work / "Cross_family"
    outer_input_root = cf / "results/09_outer_inputs"
    chip_root = cf / "results/09_outer_chip_gwas_clumping"
    merged_root = cf / "results/09_outer_imputed_clumping_merged"
    out_root = cf / "results/09_outer_crosschip_ld_filtering"

    global_flag = flag(merged_root / "complete.flag")
    for k, v in {"status": "PASS", "outer_folds": "5/5", "trait_combinations": "10/10",
                  "training_only": "PASS", "outer_test_phenotype_read": "NO"}.items():
        require(global_flag.get(k) == v, f"09f global flag mismatch: {k}={global_flag.get(k)!r}")

    passed = 0
    failed: list[str] = []
    for oi in range(1, 6):
        outer = f"{oi:02d}"
        fold_input = outer_input_root / f"outer_{outer}"
        frozen_path = fold_input / "frozen_inner_selection.tsv"
        with file(frozen_path).open(newline="") as fh:
            frozen_rows = list(csv.DictReader(fh, delimiter="\t"))
        require(len(frozen_rows) == 2 and {r.get("Trait") for r in frozen_rows} == set(TRAITS),
                f"invalid frozen selection rows for outer_{outer}")
        frozen = {r["Trait"]: r for r in frozen_rows}
        train_ids = read_train_keep(fold_input / "outer_train.keep")
        test_ids = read_ids(fold_input / "outer_test.ids")
        require(len(train_ids) + len(test_ids) == 904 and set(train_ids).isdisjoint(test_ids),
                f"outer train/test partition mismatch outer_{outer}")

        merge_fold_flag = flag(merged_root / f"outer_{outer}/complete.flag")
        require(merge_fold_flag.get("status") == "PASS" and merge_fold_flag.get("outer") == outer and
                merge_fold_flag.get("training_only") == "PASS" and
                merge_fold_flag.get("outer_test_phenotype_read") == "NO",
                f"09f fold flag mismatch outer_{outer}")
        with file(merged_root / f"outer_{outer}/merge_summary.tsv").open(newline="") as fh:
            merge_rows = list(csv.DictReader(fh, delimiter="\t"))
        require(len(merge_rows) == 2 and {r.get("trait") for r in merge_rows} == set(TRAITS),
                f"invalid 09f merge summary outer_{outer}")
        merge_counts = {r["trait"]: int(r["n_candidate_snps"]) for r in merge_rows}

        for trait in TRAITS:
            try:
                row = frozen[trait]
                best_chip = row["Best_Chip"]
                target = int(row["Best_QTN"])
                require(best_chip.startswith("chip_") and best_chip[5:].isdigit(),
                        f"invalid frozen chip strategy outer_{outer}/{trait}: {best_chip}")
                expected_chip_n = int(best_chip[5:])
                out = out_root / f"outer_{outer}/{trait}"
                f = flag(out / "complete.flag")
                expected = {
                    "status": "PASS", "version": "1", "outer": outer, "trait": trait,
                    "chip_strategy": best_chip, "chip_markers": str(expected_chip_n),
                    "selected_qtn": str(target), "ld_threshold": "0.8",
                    "n_training": str(len(train_ids)), "n_outer_test": str(len(test_ids)),
                    "n_clump_indices": str(merge_counts[trait]),
                    "training_only": "PASS", "outer_test_phenotype_read": "NO",
                }
                for k, v in expected.items():
                    require(f.get(k) == v, f"{out}/complete.flag mismatch: {k}={f.get(k)!r}, expected {v!r}")
                require((out / "input_signature.txt").is_file(), f"missing input signature: {out}")

                chip_path = chip_root / f"outer_{outer}/marker_sets/{trait}_{best_chip}.snps"
                chip = read_ids(chip_path)
                require(len(chip) == expected_chip_n, f"chip marker count mismatch {chip_path}")
                chip_set = set(chip)

                original = read_ranked(merged_root / f"outer_{outer}/{trait}_index_ranked.tsv")
                require(len(original) == merge_counts[trait], f"09f ranking count mismatch outer_{outer}/{trait}")
                nonchip = read_ranked(out / "ranked_nonchip.tsv")
                require(not (set(s for s, _ in nonchip) & chip_set), f"exact chip overlap remains in {out}/ranked_nonchip.tsv")
                expected_nonchip = [(s, p) for s, p in original if s not in chip_set]
                require([s for s, _ in nonchip] == [s for s, _ in expected_nonchip],
                        f"non-chip ranking does not preserve original P order outer_{outer}/{trait}")
                rank_summary = kv_table(out / "rank_nonchip_summary.tsv")
                require(rank_summary.get("Input_clump_indices") == str(len(original)) and
                        rank_summary.get("Nonchip_candidates_retained") == str(len(nonchip)) and
                        rank_summary.get("Exact_chip_overlaps_removed") == str(len(original) - len(nonchip)),
                        f"rank non-chip summary mismatch outer_{outer}/{trait}")

                marker = out / f"marker_sets/extra_{target}.snps"
                selected = read_ids(marker, allow_empty=(target == 0))
                selected_top = read_ids(out / "selected_top_QTN.snps", allow_empty=(target == 0))
                require(len(selected) == target and selected == selected_top,
                        f"selected QTN count/list mismatch outer_{outer}/{trait}")
                require(not (set(selected) & chip_set), f"selected QTN overlaps chip markers outer_{outer}/{trait}")
                require((out / "best_extra.txt").read_text().strip() == str(target),
                        f"best_extra.txt mismatch outer_{outer}/{trait}")

                if target == 0:
                    require(f.get("ld_filter_performed") == "NO_QTN_REQUIRED",
                            f"zero-QTN task should skip LD: outer_{outer}/{trait}")
                    filtered = read_ranked(out / "filtered_ranked_extra_snps.tsv", allow_empty=True)
                    require(len(filtered) == 0, f"zero-QTN filtered list should be empty outer_{outer}/{trait}")
                    retained_count = 0
                else:
                    require(f.get("ld_filter_performed") == "YES", f"LD not performed outer_{outer}/{trait}")
                    filtered = read_ranked(out / "filtered_ranked_extra_snps.tsv")
                    require(len(filtered) >= target, f"not enough filtered candidates outer_{outer}/{trait}")
                    require(selected == [s for s, _ in filtered[:target]],
                            f"selected SNPs are not the top frozen count of filtered ranking outer_{outer}/{trait}")
                    require(set(s for s, _ in filtered).issubset(set(s for s, _ in nonchip)),
                            f"filtered candidates not a subset of ranked non-chip candidates outer_{outer}/{trait}")
                    summary = kv_table(out / "filter_summary.tsv")
                    require(summary.get("LD_threshold") == "0.8" and
                            summary.get("Best_Extra") == str(target) and
                            summary.get("Can_make_selected_count") == "1",
                            f"LD filter summary mismatch outer_{outer}/{trait}")
                    retained_count = len(filtered)
                    require(f.get("n_retained_after_ld") == str(retained_count),
                            f"retained count in completion flag mismatch outer_{outer}/{trait}")
                    annotated = read_pairs(out / "candidate_ld_annotated.tsv")
                    require(len(annotated) == len(nonchip),
                            f"LD annotation row count mismatch outer_{outer}/{trait}")
                    n_keep = sum(r.get("Keep") == "1" for r in annotated)
                    n_filtered = sum(r.get("LD_status") == "FILTERED_HIGH_LD" for r in annotated)
                    require(n_keep == retained_count and n_filtered == int(summary.get("Filtered_high_LD", "-1")),
                            f"LD annotation/summary counts mismatch outer_{outer}/{trait}")
                    availability = read_pairs(out / "candidate_count_availability.tsv")
                    target_avail = [r for r in availability if r.get("N_QTN") == str(target)]
                    require(len(target_avail) == 1 and target_avail[0].get("Available") == "1",
                            f"frozen target unavailable in count table outer_{outer}/{trait}")

                print(f"PASS outer_{outer}/{trait}: chip={best_chip}, target_QTN={target}, "
                      f"nonchip={len(nonchip)}, retained={retained_count}, selected={len(selected)}")
                passed += 1
            except Exception as exc:
                failed.append(f"outer_{outer}/{trait}: {exc}")
                print(f"FAIL outer_{outer}/{trait}: {exc}")

    print(f"OVERALL 09g PASS={passed}/10 FAIL={len(failed)}/10")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
