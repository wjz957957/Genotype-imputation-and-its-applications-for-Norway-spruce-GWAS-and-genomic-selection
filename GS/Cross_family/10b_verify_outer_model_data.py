#!/usr/bin/env python3
"""Verify all 10 outer-fold/trait model-data bundles and optionally package transfer archives."""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import tarfile
from pathlib import Path

TRAITS = ("BB", "H6")
ROOT_REL = Path("Cross_family/results/10_outer_model_data")


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError(msg)


def read_flag(path: Path) -> dict[str, str]:
    need(path.is_file() and path.stat().st_size > 0, f"missing/empty flag: {path}")
    d = {}
    for line in path.read_text().splitlines():
        k, sep, v = line.partition("=")
        need(bool(sep) and k not in d, f"malformed flag: {path}: {line}")
        d[k] = v
    return d


def read_list(path: Path, allow_empty: bool = False) -> list[str]:
    need(path.is_file(), f"missing list: {path}")
    vals = []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        f = line.split()
        need(len(f) == 1, f"invalid list row {path}:{n}")
        vals.append(f[0])
    need(allow_empty or bool(vals), f"empty list: {path}")
    need(len(vals) == len(set(vals)), f"duplicate entries: {path}")
    return vals


def raw_header(path: Path) -> tuple[list[str], list[str]]:
    need(path.is_file() and path.stat().st_size > 0, f"missing/empty RAW: {path}")
    with path.open() as fh:
        h = fh.readline().split()
    need(h[:6] == ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] and len(h) >= 7,
         f"invalid RAW header: {path}")
    snps = [x.rsplit("_", 1)[0] for x in h[6:]]
    need(len(snps) == len(set(snps)), f"duplicate marker IDs in RAW: {path}")
    return h, snps


def raw_rows(path: Path, expected_ids: list[str], expected_cols: int, deep_check: bool) -> None:
    with path.open() as fh:
        h = fh.readline().split()
        need(len(h) == expected_cols, f"RAW header column count mismatch: {path}")
        n = 0
        for line_no, line in enumerate(fh, 2):
            f = line.split()
            need(len(f) == expected_cols, f"RAW row width mismatch at {path}:{line_no}")
            need(n < len(expected_ids), f"extra sample row in {path}:{line_no}")
            sid = expected_ids[n]
            need(f[0] == sid and f[1] == sid, f"sample order mismatch at {path}:{line_no}; expected {sid}")
            need(f[5] == "-9", f"RAW PHENOTYPE metadata must be -9 at {path}:{line_no}")
            if deep_check:
                for j, value in enumerate(f[6:], 7):
                    try:
                        x = float(value)
                    except ValueError:
                        raise RuntimeError(f"non-numeric/missing genotype at {path}:{line_no}, column {j}: {value}")
                    need(math.isfinite(x) and x in (0.0, 1.0, 2.0),
                         f"genotype is not finite 0/1/2 at {path}:{line_no}, column {j}: {value}")
            n += 1
        need(n == len(expected_ids), f"sample count mismatch in {path}: {n} != {len(expected_ids)}")


def check_pheno(path: Path, ids: list[str], header: bool, first_col: int) -> None:
    need(path.is_file() and path.stat().st_size > 0, f"missing phenotype file: {path}")
    rows = path.read_text().splitlines()
    if header:
        need(bool(rows), f"empty phenotype file: {path}")
        rows = rows[1:]
    observed = []
    for n, line in enumerate(rows, 2 if header else 1):
        f = line.split()
        need(len(f) >= first_col + 1, f"bad phenotype row {path}:{n}")
        sid = f[first_col]
        try:
            y = float(f[-1])
        except ValueError:
            raise RuntimeError(f"non-numeric phenotype {path}:{n}")
        need(math.isfinite(y), f"non-finite phenotype {path}:{n}")
        observed.append(sid)
    need(observed == ids, f"phenotype sample order mismatch: {path}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_one(work: Path, outer: str, trait: str, deep_check: bool) -> dict[str, str]:
    d = work / ROOT_REL / f"outer_{outer}" / trait
    flag = read_flag(d / "complete.flag")
    need(flag.get("status") == "PASS" and flag.get("outer") == outer and flag.get("trait") == trait,
         f"invalid completion flag: {d}")
    for k in ("test_phenotype_used_for_training", "test_phenotype_in_raw", "raw_missing_genotypes"):
        expected = {"test_phenotype_used_for_training":"NO", "test_phenotype_in_raw":"NO", "raw_missing_genotypes":"0"}[k]
        need(flag.get(k) == expected, f"invalid {k} in {d}/complete.flag")
    with (d / "final_genotype_manifest.tsv").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    need(len(rows) == 1, f"manifest must have exactly one row: {d}")
    m = rows[0]
    n_chip, n_qtn = int(m["N_chip"]), int(m["N_QTN"])
    n_train, n_test = int(m["N_train"]), int(m["N_test"])
    need(m["Best_Chip"] == flag.get("chip_strategy") and str(n_qtn) == flag.get("selected_qtn"),
         f"manifest/flag selection mismatch: {d}")
    need(n_chip + n_qtn == int(m["N_total"]) and n_train + n_test == 904,
         f"manifest count mismatch: {d}")

    train_ids = read_list(d / "train.ids")
    test_ids = read_list(d / "test.ids")
    need(len(train_ids) == n_train and len(test_ids) == n_test and not (set(train_ids) & set(test_ids)),
         f"sample partition mismatch: {d}")
    chip_ids = read_list(d / "selected_chip.snps")
    qtn_ids = read_list(d / "selected_qtn.snps", allow_empty=(n_qtn == 0))
    chip_order = read_list(d / "chip_model_order.snps")
    qtn_order = read_list(d / "qtn_model_order.snps", allow_empty=(n_qtn == 0))
    combined_order = read_list(d / "combined_model_order.snps")
    need(len(chip_ids) == n_chip and len(qtn_ids) == n_qtn and not (set(chip_ids) & set(qtn_ids)),
         f"selected marker count/overlap mismatch: {d}")
    need(set(chip_order) == set(chip_ids) and len(chip_order) == n_chip,
         f"chip model order mismatch: {d}")
    need(qtn_order == [] if n_qtn == 0 else (set(qtn_order) == set(qtn_ids) and len(qtn_order) == n_qtn),
         f"QTN model order mismatch: {d}")
    need(combined_order == chip_order + qtn_order, f"combined marker order mismatch: {d}")

    expected_sets = {
        "chip_train.raw": chip_order, "chip_test.raw": chip_order,
        "combined_train.raw": combined_order, "combined_test.raw": combined_order
    }
    if n_qtn:
        expected_sets.update({"extra_train.raw": qtn_order, "extra_test.raw": qtn_order})
    else:
        need(not (d / "extra_train.raw").exists() and not (d / "extra_test.raw").exists(),
             f"extra RAW should not exist for zero-QTN fold: {d}")

    for name, exp in expected_sets.items():
        h, snps = raw_header(d / name)
        need(snps == exp, f"RAW SNP order mismatch: {d/name}")
        expected_ids = train_ids if name.endswith("_train.raw") else test_ids
        raw_rows(d / name, expected_ids, len(h), deep_check)
    check_pheno(d / "train.pheno.tsv", train_ids, header=True, first_col=1)
    check_pheno(d / "deepgs_train_pheno.tsv", train_ids, header=True, first_col=0)
    check_pheno(d / "train.pheno", train_ids, header=False, first_col=1)
    check_pheno(d / "evaluation_only/test.pheno.tsv", test_ids, header=True, first_col=0)

    if n_qtn == 0:
        need((d / "chip_train.raw").read_bytes() == (d / "combined_train.raw").read_bytes(),
             f"zero-QTN training combined file differs from CHIP: {d}")
        need((d / "chip_test.raw").read_bytes() == (d / "combined_test.raw").read_bytes(),
             f"zero-QTN test combined file differs from CHIP: {d}")

    print(f"PASS outer_{outer}/{trait}: chip={n_chip}, QTN={n_qtn}, total={n_chip+n_qtn}, "
          f"train={n_train}, test={n_test}")
    return {
        "outer": outer, "trait": trait, "best_chip": m["Best_Chip"],
        "n_chip": str(n_chip), "n_qtn": str(n_qtn), "n_total": str(n_chip+n_qtn),
        "n_train": str(n_train), "n_test": str(n_test), "status": "PASS"
    }


def make_bundles(work: Path) -> None:
    root = work / ROOT_REL
    transfer = work / "Cross_family/transfer"
    transfer.mkdir(parents=True, exist_ok=True)
    model_tar = transfer / "DeepGS_model_inputs_5outer.tar.gz"
    eval_tar = transfer / "DeepGS_evaluation_only_test_phenotypes_5outer.tar.gz"
    readme = transfer / "README_DeepGS_transfer.txt"
    readme.write_text(
        "DeepGS transfer bundles for current 5-outer-fold analysis.\n"
        "1) Transfer/extract DeepGS_model_inputs_5outer.tar.gz first. It contains training phenotypes,\n"
        "   train/test genotypes, sample IDs, SNP lists, marker order, manifests and QC flags.\n"
        "   It EXCLUDES all outer-test phenotype values.\n"
        "2) Keep DeepGS_evaluation_only_test_phenotypes_5outer.tar.gz on the old server until DeepGS\n"
        "   predictions for all outer-test individuals have been written. Transfer it only for final\n"
        "   independent evaluation. Never use these values for fitting, validation, feature selection\n"
        "   or tuning.\n"
        "3) The RAW PHENOTYPE metadata field is set to -9; the real training phenotype is separate.\n"
    )
    with tarfile.open(model_tar, "w:gz") as tf:
        for p in sorted(root.rglob("*")):
            if not p.is_file() or "evaluation_only" in p.parts:
                continue
            if p.name.startswith("plink_extract_"):
                continue
            tf.add(p, arcname=str(p.relative_to(work)))
        tf.add(readme, arcname=str(readme.relative_to(work)))
    with tarfile.open(eval_tar, "w:gz") as tf:
        for p in sorted(root.glob("outer_*/**/evaluation_only/test.pheno.tsv")):
            tf.add(p, arcname=str(p.relative_to(work)))
        tf.add(readme, arcname=str(readme.relative_to(work)))
    sums = transfer / "SHA256SUMS.txt"
    sums.write_text(
        f"{sha256(model_tar)}  {model_tar.name}\n"
        f"{sha256(eval_tar)}  {eval_tar.name}\n"
    )
    print(f"Created model-input bundle: {model_tar}")
    print(f"Created separate evaluation-only bundle: {eval_tar}")
    print(f"SHA256: {sums}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--deep-check", action="store_true",
                    help="scan every genotype cell for missing/non-0/1/2 values")
    ap.add_argument("--make-transfer-bundles", action="store_true",
                    help="after successful 10/10 verification, create separate model/evaluation archives")
    a = ap.parse_args()
    work = Path(a.work).resolve()
    rows = []
    for oi in range(1, 6):
        outer = f"{oi:02d}"
        for trait in TRAITS:
            rows.append(verify_one(work, outer, trait, a.deep_check))
    summary = work / ROOT_REL / "model_data_summary.tsv"
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["outer","trait","best_chip","n_chip","n_qtn","n_total","n_train","n_test","status"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("OVERALL 10a VERIFY PASS=10/10 FAIL=0/10")
    if a.make_transfer_bundles:
        make_bundles(work)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        import sys
        print(f"ERROR 10b verify: {exc}", file=sys.stderr)
        sys.exit(1)
