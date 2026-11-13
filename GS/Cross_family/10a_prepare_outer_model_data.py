#!/usr/bin/env python3
"""Prepare leakage-safe outer-fold genotype inputs for RRBLUP and DeepGS.

This stage only exports already-frozen CHIP and QTN sets. It does not rerun GWAS,
clumping, LD filtering, or feature selection. PLINK RAW metadata PHENOTYPE is
forced to -9; phenotypes are stored in separate files. Test phenotypes are
written only to an evaluation_only subdirectory, excluded from the DeepGS
training-input archive.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

TRAITS = ("BB", "H6")
VALID_QTN = {0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500, 600, 700}
WORK_REL = Path("Cross_family")
MASTER_REL = Path("08_imputed_GWAS/00_reference/all_904_genotypes")
INPUT_REL = Path("Cross_family/results/09_outer_inputs")
CHIP_REL = Path("Cross_family/results/09_outer_chip_gwas_clumping")
IMPUTED_REL = Path("Cross_family/results/09_outer_imputed_gwas_inputs")
QTN_REL = Path("Cross_family/results/09_outer_crosschip_ld_filtering")
SELECTION_REL = Path("Cross_family/results/08_inner_qtn_selection")
CONFIG_REL = Path("Cross_family/config/best_chip_by_outer_trait.tsv")
OUT_REL = Path("Cross_family/results/10_outer_model_data")


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError(msg)


def read_flag(path: Path) -> dict[str, str]:
    need(path.is_file() and path.stat().st_size > 0, f"missing/empty flag: {path}")
    out: dict[str, str] = {}
    for n, line in enumerate(path.read_text().splitlines(), 1):
        key, sep, value = line.partition("=")
        need(bool(sep) and key not in out, f"malformed/duplicate flag {path}:{n}: {line}")
        out[key] = value
    return out


def check_flag(path: Path, expected: dict[str, str]) -> dict[str, str]:
    got = read_flag(path)
    for k, v in expected.items():
        need(got.get(k) == v, f"{path}: expected {k}={v}, found {got.get(k)!r}")
    return got


def read_tsv(path: Path) -> list[dict[str, str]]:
    need(path.is_file() and path.stat().st_size > 0, f"missing/empty TSV: {path}")
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_onecol(path: Path, allow_empty: bool = False) -> list[str]:
    need(path.is_file(), f"missing list: {path}")
    vals = []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        need(len(fields) == 1, f"expected one field at {path}:{n}")
        vals.append(fields[0])
    need(allow_empty or bool(vals), f"empty list: {path}")
    need(len(vals) == len(set(vals)), f"duplicate IDs in {path}")
    return vals


def read_keep(path: Path) -> list[str]:
    need(path.is_file() and path.stat().st_size > 0, f"missing keep file: {path}")
    ids = []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        f = line.split()
        need(len(f) == 2 and f[0] == f[1], f"invalid FID/IID at {path}:{n}")
        ids.append(f[1])
    need(len(ids) == len(set(ids)) and bool(ids), f"duplicate/empty IDs: {path}")
    return ids


def read_threecol_pheno(path: Path) -> dict[str, float]:
    need(path.is_file() and path.stat().st_size > 0, f"missing phenotype: {path}")
    out: dict[str, float] = {}
    for n, line in enumerate(path.read_text().splitlines(), 1):
        f = line.split()
        need(len(f) == 3 and f[0] == f[1], f"expected FID IID Y at {path}:{n}")
        try:
            y = float(f[2])
        except ValueError as exc:
            raise RuntimeError(f"non-numeric phenotype at {path}:{n}") from exc
        need(math.isfinite(y) and y != -9.0, f"missing/non-finite phenotype at {path}:{n}")
        need(f[1] not in out, f"duplicate phenotype ID at {path}:{n}")
        out[f[1]] = y
    need(bool(out), f"empty phenotype: {path}")
    return out


def read_full_pheno(path: Path) -> dict[str, float]:
    """Read the original three-column FID IID phenotype file, with or without header."""
    need(path.is_file() and path.stat().st_size > 0, f"missing full phenotype source: {path}")
    out: dict[str, float] = {}
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        f = line.split()
        need(len(f) == 3, f"expected 3 phenotype columns at {path}:{n}")
        try:
            y = float(f[2])
        except ValueError:
            if n == 1:
                continue
            raise RuntimeError(f"non-numeric phenotype at {path}:{n}")
        need(f[0] == f[1], f"FID/IID mismatch at {path}:{n}")
        need(math.isfinite(y) and y != -9.0, f"missing/non-finite phenotype at {path}:{n}")
        need(f[1] not in out, f"duplicate phenotype ID at {path}:{n}")
        out[f[1]] = y
    need(bool(out), f"no phenotype records: {path}")
    return out


def read_master_fam(path: Path) -> list[str]:
    need(path.is_file() and path.stat().st_size > 0, f"missing master FAM: {path}")
    ids = []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        f = line.split()
        need(len(f) >= 2 and f[0] == f[1], f"master FAM must have FID=IID at {path}:{n}")
        ids.append(f[1])
    need(len(ids) == 904, f"expected 904 master samples, found {len(ids)}")
    need(len(ids) == len(set(ids)), "duplicate master FAM sample IDs")
    return ids


def fold2(value: str) -> str:
    return f"{int(value):02d}"


def load_selection(work: Path, outer: str, trait: str) -> tuple[str, int, list[str], list[str], dict[str, str]]:
    cf = work / WORK_REL
    inp = cf / "results/09_outer_inputs" / f"outer_{outer}"
    chip_fold = cf / "results/09_outer_chip_gwas_clumping" / f"outer_{outer}"
    imp_fold = cf / "results/09_outer_imputed_gwas_inputs" / f"outer_{outer}"
    qtn_fold = cf / "results/09_outer_crosschip_ld_filtering" / f"outer_{outer}" / trait

    check_flag(inp / "prepare.complete", {
        "status": "PASS", "outer": outer, "outer_test_phenotype_read": "NO"
    })
    check_flag(chip_fold / "complete.flag", {
        "status": "PASS", "outer": outer, "training_only": "PASS",
        "outer_test_phenotype_read": "NO", "selected_outer_chip_SNPs": "PASS"
    })
    check_flag(imp_fold / "complete.flag", {
        "status": "PASS", "outer": outer, "training_only": "PASS",
        "outer_test_phenotype_read": "NO"
    })
    gflag = check_flag(qtn_fold / "complete.flag", {
        "status": "PASS", "outer": outer, "trait": trait,
        "training_only": "PASS", "outer_test_phenotype_read": "NO"
    })

    frozen_rows = read_tsv(inp / "frozen_inner_selection.tsv")
    fr = [r for r in frozen_rows if fold2(r.get("Outer", "0")) == outer and r.get("Trait") == trait]
    need(len(fr) == 1, f"expected one frozen selection row for outer_{outer}/{trait}")
    best_chip = fr[0]["Best_Chip"]
    target = int(fr[0]["Best_QTN"])
    need(best_chip.startswith("chip_") and best_chip[5:].isdigit(), f"invalid frozen chip: {best_chip}")
    need(int(best_chip[5:]) in {50,100,200,300,400,500,1000,2000,3000,4000,5000,8000,10000,20000,30000,40000},
         f"unsupported frozen chip strategy: {best_chip}")
    need(target in VALID_QTN, f"unsupported frozen QTN count: {target}")
    need(gflag.get("chip_strategy") == best_chip, f"09g chip strategy mismatch for outer_{outer}/{trait}")
    need(int(gflag.get("selected_qtn", "-1")) == target, f"09g QTN count mismatch for outer_{outer}/{trait}")

    # Cross-check the authoritative frozen selection against the upstream selection tables.
    config_rows = read_tsv(work / CONFIG_REL)
    config = [r for r in config_rows if fold2(r.get("Outer", "0")) == outer and r.get("Trait") == trait]
    need(len(config) == 1 and config[0].get("Best_Chip") == best_chip,
         f"best-chip config does not match frozen selection for outer_{outer}/{trait}")
    qtn_rows = read_tsv(work / SELECTION_REL / "selected_qtn_counts.tsv")
    qrows = [r for r in qtn_rows if fold2(r.get("outer", "0")) == outer and r.get("trait") == trait]
    need(len(qrows) == 1 and int(qrows[0].get("n_qtn", "-1")) == target,
         f"QTN summary does not match frozen selection for outer_{outer}/{trait}")

    chip_file = chip_fold / "marker_sets" / f"{trait}_{best_chip}.snps"
    qtn_file = qtn_fold / "marker_sets" / f"extra_{target}.snps"
    chip_ids = read_onecol(chip_file)
    qtn_ids = read_onecol(qtn_file, allow_empty=(target == 0))
    need(len(chip_ids) == int(best_chip[5:]), f"chip marker count mismatch: {chip_file}")
    need(len(qtn_ids) == target, f"QTN marker count mismatch: {qtn_file}")
    need(not (set(chip_ids) & set(qtn_ids)), f"chip/QTN overlap in outer_{outer}/{trait}")
    return best_chip, target, chip_ids, qtn_ids, {
        "input_dir": str(inp), "chip_file": str(chip_file), "qtn_file": str(qtn_file),
        "full_pheno_BB": str(work / "00_input/Trait-Adj_Budburst.txt"),
        "full_pheno_H6": str(work / "00_input/Trait-Adj_H6.txt")
    }


def validate_partition(work: Path, outer: str, master_ids: list[str]) -> tuple[list[str], list[str]]:
    inp = work / "Cross_family/results/09_outer_inputs" / f"outer_{outer}"
    train_set = set(read_keep(inp / "outer_train.keep"))
    test_set = set(read_onecol(inp / "outer_test.ids"))
    need(not (train_set & test_set), f"train/test overlap outer_{outer}")
    need(train_set | test_set == set(master_ids), f"train/test are not a complete master partition outer_{outer}")
    train_order = [sid for sid in master_ids if sid in train_set]
    test_order = [sid for sid in master_ids if sid in test_set]
    need(len(train_order) + len(test_order) == 904, f"bad partition count outer_{outer}")
    return train_order, test_order


def phenotype_paths(work: Path, trait: str) -> tuple[Path, Path]:
    inp = work / "Cross_family/results/09_outer_inputs"
    if trait == "BB":
        return inp, work / "00_input/Trait-Adj_Budburst.txt"
    return inp, work / "00_input/Trait-Adj_H6.txt"


def validate_phenotypes(work: Path, outer: str, trait: str, train_order: list[str],
                        test_order: list[str]) -> tuple[dict[str, float], dict[str, float], Path, Path]:
    inp = work / "Cross_family/results/09_outer_inputs" / f"outer_{outer}"
    train_path = inp / f"{trait}.train.pheno"
    train_map = read_threecol_pheno(train_path)
    need(set(train_map) == set(train_order), f"training phenotype IDs differ from outer-training IDs: {train_path}")
    _, full_path = phenotype_paths(work, trait)
    full_map = read_full_pheno(full_path)
    need(set(train_order) <= set(full_map), f"original phenotype source lacks training IDs: {outer}/{trait}")
    need(set(test_order) <= set(full_map), f"original phenotype source lacks test IDs: {outer}/{trait}")
    for sid in train_order:
        a, b = train_map[sid], full_map[sid]
        need(abs(a-b) <= 1e-8 * max(1.0, abs(a), abs(b)),
             f"training phenotype differs from original source for {outer}/{trait}/{sid}: {a} vs {b}")
    return train_map, full_map, train_path, full_path


def check_markers_in_bim(bim: Path, wanted: set[str]) -> None:
    found: Counter[str] = Counter()
    with bim.open() as fh:
        for n, line in enumerate(fh, 1):
            f = line.split()
            need(len(f) == 6, f"invalid BIM row {bim}:{n}")
            if f[1] in wanted:
                found[f[1]] += 1
    missing = sorted(wanted - set(found))
    duplicated = sorted(s for s, n in found.items() if n != 1)
    need(not missing, f"{len(missing)} selected SNPs absent from master BIM; examples={missing[:10]}")
    need(not duplicated, f"selected SNP IDs duplicated in master BIM; examples={duplicated[:10]}")


def preflight(work: Path) -> None:
    master = work / MASTER_REL
    for ext in ("bed", "bim", "fam"):
        need(Path(str(master) + f".{ext}").is_file(), f"missing master {ext}: {master}.{ext}")
    master_ids = read_master_fam(Path(str(master) + ".fam"))
    all_wanted: set[str] = set()
    for oi in range(1, 6):
        outer = f"{oi:02d}"
        train_order, test_order = validate_partition(work, outer, master_ids)
        for trait in TRAITS:
            best_chip, target, chip_ids, qtn_ids, _ = load_selection(work, outer, trait)
            train_map, full_map, train_path, full_path = validate_phenotypes(
                work, outer, trait, train_order, test_order
            )
            all_wanted.update(chip_ids)
            all_wanted.update(qtn_ids)
            print(f"PASS preflight outer_{outer}/{trait}: chip={best_chip} ({len(chip_ids)}), "
                  f"QTN={target}, train={len(train_order)}, test={len(test_order)}, "
                  "phenotype_match=PASS")
    check_markers_in_bim(Path(str(master) + ".bim"), all_wanted)
    print(f"PASS master BIM: all {len(all_wanted)} distinct selected chip/QTN SNP IDs exist exactly once")
    print("OVERALL 10a PREFLIGHT PASS=10/10 FAIL=0/10")


def marker_id(raw_column: str) -> str:
    need("_" in raw_column, f"PLINK RAW marker column lacks allele suffix: {raw_column}")
    return raw_column.rsplit("_", 1)[0]


def write_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w") as fh:
        for line in lines:
            fh.write(line + "\n")
    os.replace(tmp, path)


def write_pheno_files(stage: Path, train_order: list[str], test_order: list[str],
                      train_map: dict[str, float], full_map: dict[str, float],
                      train_source: Path, test_source: Path) -> None:
    (stage / "evaluation_only").mkdir(parents=True, exist_ok=True)
    with (stage / "train.pheno.tsv").open("w") as fh, \
         (stage / "deepgs_train_pheno.tsv").open("w") as dh, \
         (stage / "train.pheno").open("w") as rh, \
         (stage / "evaluation_only/test.pheno.tsv").open("w") as th:
        fh.write("FID\tIID\tphenotype\n")
        dh.write("ID\tphenotype\n")
        th.write("ID\tphenotype\n")
        for sid in train_order:
            y = train_map[sid]
            fh.write(f"{sid}\t{sid}\t{y:.17g}\n")
            dh.write(f"{sid}\t{y:.17g}\n")
            rh.write(f"{sid}\t{sid}\t{y:.17g}\n")
        for sid in test_order:
            th.write(f"{sid}\t{full_map[sid]:.17g}\n")
    write_lines(stage / "phenotype_sources.txt", [
        f"training_phenotype_source={train_source}",
        f"evaluation_only_test_phenotype_source={test_source}",
        "test_phenotype_used_for_training=NO",
        "test_phenotype_used_for_feature_selection=NO"
    ])


def extract_with_plink(plink: str, master: Path, keep: Path, extract: Path,
                       prefix: Path, cpus: int, memory_mb: int) -> Path:
    cmd = [
        plink, "--bfile", str(master), "--keep", str(keep), "--extract", str(extract),
        "--keep-allele-order", "--allow-extra-chr", "--recode", "A",
        "--threads", str(cpus), "--memory", str(memory_mb), "--out", str(prefix)
    ]
    print("[PLINK]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    raw = Path(str(prefix) + ".raw")
    need(raw.is_file() and raw.stat().st_size > 0, f"PLINK did not create RAW: {raw}")
    return raw


def split_raw(raw: Path, stage: Path, split: str, expected_ids: list[str],
              chip_ids: list[str], qtn_ids: list[str]) -> tuple[list[str], list[str]]:
    with raw.open() as fh:
        first = fh.readline()
    header = first.split()
    need(header[:6] == ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] and len(header) >= 7,
         f"invalid PLINK RAW header: {raw}")
    ids_from_header = [marker_id(x) for x in header[6:]]
    need(len(ids_from_header) == len(set(ids_from_header)), f"duplicate SNP IDs in RAW header: {raw}")
    expected_union = set(chip_ids) | set(qtn_ids)
    need(set(ids_from_header) == expected_union,
         f"RAW marker set mismatch {raw}: expected {len(expected_union)}, got {len(ids_from_header)}")
    chip_idx = [i for i, sid in enumerate(ids_from_header) if sid in set(chip_ids)]
    qtn_idx = [i for i, sid in enumerate(ids_from_header) if sid in set(qtn_ids)]
    chip_order = [ids_from_header[i] for i in chip_idx]
    qtn_order = [ids_from_header[i] for i in qtn_idx]
    chip_header = header[:6] + [header[6+i] for i in chip_idx]
    qtn_header = header[:6] + [header[6+i] for i in qtn_idx]
    combined_header = header[:6] + chip_header[6:] + qtn_header[6:]

    chip_out = stage / f"chip_{split}.raw"
    combined_out = stage / f"combined_{split}.raw"
    extra_out = stage / f"extra_{split}.raw"
    missing_count = 0
    first_missing = None
    row_count = 0

    with raw.open() as src, chip_out.open("w") as cf, combined_out.open("w") as combf:
        if qtn_ids:
            ef = extra_out.open("w")
        else:
            ef = None
        try:
            src.readline()
            cf.write("\t".join(chip_header) + "\n")
            combf.write("\t".join(combined_header) + "\n")
            if ef is not None:
                ef.write("\t".join(qtn_header) + "\n")
            for line_no, line in enumerate(src, 2):
                f = line.split()
                need(len(f) == len(header), f"RAW column count mismatch at {raw}:{line_no}")
                need(row_count < len(expected_ids), f"unexpected extra sample row at {raw}:{line_no}")
                expected_id = expected_ids[row_count]
                need(f[0] == expected_id and f[1] == expected_id,
                     f"sample order mismatch at {raw}:{line_no}; expected {expected_id}, got {f[0]}/{f[1]}")
                meta = f[:6]
                meta[5] = "-9"  # prevent genotype RAW metadata from carrying phenotype information
                vals = f[6:]
                for j, value in enumerate(vals):
                    if value.upper() in {"NA", "NAN", ".", "-9"}:
                        missing_count += 1
                        if first_missing is None:
                            first_missing = (line_no, ids_from_header[j], expected_id, value)
                    else:
                        try:
                            x = float(value)
                        except ValueError:
                            raise RuntimeError(f"non-numeric genotype {value!r} at {raw}:{line_no}")
                        need(math.isfinite(x) and x in (0.0, 1.0, 2.0),
                             f"unexpected genotype dosage {value!r} at {raw}:{line_no}, SNP={ids_from_header[j]}")
                chip_vals = [vals[i] for i in chip_idx]
                qtn_vals = [vals[i] for i in qtn_idx]
                cf.write("\t".join(meta + chip_vals) + "\n")
                combf.write("\t".join(meta + chip_vals + qtn_vals) + "\n")
                if ef is not None:
                    ef.write("\t".join(meta + qtn_vals) + "\n")
                row_count += 1
        finally:
            if ef is not None:
                ef.close()

    need(row_count == len(expected_ids), f"sample count mismatch in {raw}: {row_count} != {len(expected_ids)}")
    if missing_count:
        raise RuntimeError(
            f"{missing_count} missing/invalid genotype cells in {raw}; first={first_missing}. "
            "No data were published. Do not impute silently; choose a training-only imputation policy first."
        )
    write_lines(stage / f"{split}.ids", expected_ids)
    if split == "train":
        write_lines(stage / "chip_model_order.snps", chip_order)
        write_lines(stage / "qtn_model_order.snps", qtn_order)
        write_lines(stage / "combined_model_order.snps", chip_order + qtn_order)
    return chip_order, qtn_order


def validate_existing(out: Path, outer: str, trait: str, best_chip: str, target: int) -> bool:
    flag_path = out / "complete.flag"
    if not flag_path.is_file():
        return False
    f = read_flag(flag_path)
    expected = {"status":"PASS", "outer":outer, "trait":trait,
                "chip_strategy":best_chip, "selected_qtn":str(target)}
    for k, v in expected.items():
        need(f.get(k) == v, f"existing output selection mismatch in {flag_path}: {k}")
    for name in ("chip_train.raw","chip_test.raw","combined_train.raw","combined_test.raw",
                 "train.pheno.tsv","deepgs_train_pheno.tsv","train.ids","test.ids",
                 "selected_chip.snps","selected_qtn.snps","final_genotype_manifest.tsv"):
        need((out/name).is_file(), f"complete flag exists but required output is missing: {out/name}")
    return True


def prepare_one(work: Path, outer: str, trait: str, plink: str, cpus: int, memory_mb: int) -> None:
    master = work / MASTER_REL
    master_ids = read_master_fam(Path(str(master) + ".fam"))
    train_order, test_order = validate_partition(work, outer, master_ids)
    best_chip, target, chip_ids, qtn_ids, paths = load_selection(work, outer, trait)
    train_map, full_map, train_pheno_path, full_pheno_path = validate_phenotypes(
        work, outer, trait, train_order, test_order
    )
    wanted = set(chip_ids) | set(qtn_ids)
    check_markers_in_bim(Path(str(master) + ".bim"), wanted)

    out = work / OUT_REL / f"outer_{outer}" / trait
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        if validate_existing(out, outer, trait, best_chip, target):
            print(f"[SKIP] verified existing outer_{outer}/{trait}")
            return
        raise RuntimeError(f"output directory exists without a matching complete flag: {out}; inspect manually")

    stage = out.parent / f".{trait}.stage10a.{os.getpid()}"
    need(not stage.exists(), f"staging directory already exists: {stage}")
    stage.mkdir(parents=True)
    try:
        write_lines(stage / "selected_chip.snps", chip_ids)
        write_lines(stage / "selected_qtn.snps", qtn_ids)
        write_lines(stage / "union.snps", chip_ids + qtn_ids)
        write_lines(stage / "train.keep", [f"{x}\t{x}" for x in train_order])
        write_lines(stage / "test.keep", [f"{x}\t{x}" for x in test_order])
        write_pheno_files(stage, train_order, test_order, train_map, full_map, train_pheno_path, full_pheno_path)

        orders = {}
        for split, ids in (("train", train_order), ("test", test_order)):
            prefix = stage / f".plink_{split}_union"
            raw = extract_with_plink(plink, master, stage / f"{split}.keep",
                                     stage / "union.snps", prefix, cpus, memory_mb)
            chip_order, qtn_order = split_raw(raw, stage, split, ids, chip_ids, qtn_ids)
            orders[split] = (chip_order, qtn_order)
            raw.unlink()
            log = Path(str(prefix) + ".log")
            if log.exists():
                log.rename(stage / f"plink_extract_{split}.log")
            nosex = Path(str(prefix) + ".nosex")
            if nosex.exists():
                nosex.unlink()

        need(orders["train"] == orders["test"], "train/test CHIP/QTN marker order differs")
        need((stage / "chip_train.raw").read_text().splitlines()[0] ==
             (stage / "chip_test.raw").read_text().splitlines()[0],
             "chip train/test header mismatch")
        need((stage / "combined_train.raw").read_text().splitlines()[0] ==
             (stage / "combined_test.raw").read_text().splitlines()[0],
             "combined train/test header mismatch")
        if target == 0:
            need((stage / "combined_train.raw").read_bytes() == (stage / "chip_train.raw").read_bytes(),
                 "zero-QTN combined training RAW is not identical to CHIP")
            need((stage / "combined_test.raw").read_bytes() == (stage / "chip_test.raw").read_bytes(),
                 "zero-QTN combined test RAW is not identical to CHIP")

        # Ordered SNP manifest uses the exact column order written into model matrices.
        with (stage / "marker_order.tsv").open("w") as fh:
            fh.write("Set\tOrder\tSNP\tRAW_column\n")
            for set_name, raw_name in (("chip","chip_train.raw"),("qtn","extra_train.raw"),("combined","combined_train.raw")):
                raw_path = stage / raw_name
                if not raw_path.exists():
                    continue
                hdr = raw_path.open().readline().split()
                for idx, col in enumerate(hdr[6:], 1):
                    fh.write(f"{set_name}\t{idx}\t{marker_id(col)}\t{col}\n")

        manifest = {
            "Outer": outer, "Trait": trait, "Best_Chip": best_chip,
            "N_chip": str(len(chip_ids)), "Best_QTN": str(target), "N_QTN": str(len(qtn_ids)),
            "N_total": str(len(chip_ids) + len(qtn_ids)),
            "N_train": str(len(train_order)), "N_test": str(len(test_order)),
            "Raw_missing_genotypes": "0", "Genotype_encoding": "0/1/2",
            "Phenotype_source_match": "PASS", "Test_phenotype_in_raw": "NO",
            "Test_phenotype_use": "evaluation_only"
        }
        with (stage / "final_genotype_manifest.tsv").open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(manifest), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(manifest)

        write_lines(stage / "complete.flag", [
            "status=PASS", "version=1", f"outer={outer}", f"trait={trait}",
            f"chip_strategy={best_chip}", f"selected_qtn={target}",
            f"n_training={len(train_order)}", f"n_outer_test={len(test_order)}",
            f"n_chip={len(chip_ids)}", f"n_qtn={target}", f"n_total={len(chip_ids)+target}",
            "training_phenotype_match=PASS", "test_phenotype_used_for_training=NO",
            "test_phenotype_in_raw=NO", "test_phenotype_saved_separately=YES",
            "raw_missing_genotypes=0", "status_detail=all genotype cells are finite 0/1/2"
        ])
        stage.rename(out)
        print(f"PASS outer_{outer}/{trait}: chip={best_chip} ({len(chip_ids)}), "
              f"QTN={target}, combined={len(chip_ids)+target}, "
              f"train={len(train_order)}, test={len(test_order)}")
    except Exception:
        print(f"ERROR: staging data retained for diagnosis at {stage}", file=sys.stderr)
        raise


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True, help="GS-reviewer project root")
    ap.add_argument("--plink", default="plink", help="PLINK 1.9 executable")
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--memory-mb", type=int, default=24000)
    ap.add_argument("--preflight", action="store_true", help="validate all 10 frozen input combinations; do not extract")
    ap.add_argument("--outer", help="outer fold 01..05 for a single task")
    ap.add_argument("--trait", choices=TRAITS, help="BB or H6 for a single task")
    args = ap.parse_args()
    work = Path(args.work).resolve()
    need(args.cpus >= 1 and args.memory_mb >= 1000, "invalid CPU/memory setting")
    if args.preflight:
        need(not args.outer and not args.trait, "--preflight cannot be combined with --outer/--trait")
        preflight(work)
        return
    need(args.outer is not None and args.trait is not None, "provide --outer and --trait, or use --preflight")
    outer = fold2(args.outer)
    need(outer in {f"{i:02d}" for i in range(1,6)}, f"invalid outer fold: {args.outer}")
    prepare_one(work, outer, args.trait, args.plink, args.cpus, args.memory_mb)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR 10a: {exc}", file=sys.stderr)
        sys.exit(1)
