#!/usr/bin/env python3

import argparse
import csv
import math
import os
from pathlib import Path

EXPECTED_COUNTS = {
    "01": {"01": (483, 246), "02": (493, 236), "03": (482, 247)},
    "02": {"01": (483, 243), "02": (492, 234), "03": (477, 249)},
    "03": {"01": (494, 234), "02": (479, 249), "03": (483, 245)},
    "04": {"01": (473, 246), "02": (482, 237), "03": (483, 236)},
    "05": {"01": (469, 245), "02": (474, 240), "03": (485, 229)},
}

TRAITS = {
    "BB": ("Trait-Adj_Budburst.txt", "Adj_Budburst"),
    "H6": ("Trait-Adj_H6.txt", "Adj_H6"),
}


def read_ids(path):
    if not path.is_file():
        raise ValueError(f"Missing sample list: {path}")

    lines = path.read_text(encoding="utf-8-sig").splitlines()
    ids = []

    for line_no, line in enumerate(lines, 1):
        sample = line.strip()

        if not sample or len(sample.split()) != 1:
            raise ValueError(
                f"Invalid sample ID at {path}:{line_no}"
            )

        ids.append(sample)

    if not ids:
        raise ValueError(f"Empty sample list: {path}")

    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate sample IDs: {path}")

    return ids


def read_phenotypes(path, trait_column):
    if not path.is_file():
        raise ValueError(f"Missing phenotype file: {path}")

    values = {}

    with path.open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"FID", "IID", trait_column}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"{path}: missing columns {sorted(missing)}"
            )

        for line_no, row in enumerate(reader, 2):
            fid = (row["FID"] or "").strip()
            iid = (row["IID"] or "").strip()
            value = (row[trait_column] or "").strip()

            if not fid or fid != iid:
                raise ValueError(
                    f"{path}:{line_no}: FID/IID mismatch"
                )

            if iid in values:
                raise ValueError(
                    f"{path}:{line_no}: duplicate IID {iid}"
                )

            try:
                number = float(value)
            except ValueError:
                raise ValueError(
                    f"{path}:{line_no}: invalid phenotype {value}"
                )

            if not math.isfinite(number):
                raise ValueError(
                    f"{path}:{line_no}: non-finite phenotype"
                )

            # Preserve the original numeric text.
            values[iid] = value

    if len(values) != 904:
        raise ValueError(
            f"{path}: expected 904 records, found {len(values)}"
        )

    return values


def write_pheno(path, ids, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")

    with temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(
            handle, delimiter="\t", lineterminator="\n"
        )

        for iid in ids:
            writer.writerow((iid, iid, values[iid]))

    os.replace(temp, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--work",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    args = parser.parse_args()

    work = args.work.resolve()
    sample_root = work / "Cross_family" / "cv_samples_v1"
    phenotype_root = work / "00_input"
    output_root = (
        work / "Cross_family" / "results" / "02_inner_phenotypes"
    )

    # 1. Validate every outer/inner sample partition before writing.
    folds = {}
    reference_ids = None

    for outer in EXPECTED_COUNTS:
        outer_dir = sample_root / f"outer_{outer}"

        outer_train = read_ids(
            outer_dir / "outer_train.samples"
        )
        outer_test = read_ids(
            outer_dir / "outer_test.samples"
        )

        train_set = set(outer_train)
        test_set = set(outer_test)

        if train_set & test_set:
            raise ValueError(f"Outer {outer}: train/test overlap")

        if len(train_set | test_set) != 904:
            raise ValueError(
                f"Outer {outer}: expected 904 unique samples"
            )

        if reference_ids is None:
            reference_ids = train_set | test_set
        elif train_set | test_set != reference_ids:
            raise ValueError(
                f"Outer {outer}: sample universe differs"
            )

        for inner, expected in EXPECTED_COUNTS[outer].items():
            inner_dir = outer_dir / f"inner_{inner}"

            inner_train = read_ids(
                inner_dir / "inner_train.samples"
            )
            inner_valid = read_ids(
                inner_dir / "inner_valid.samples"
            )

            itr = set(inner_train)
            iva = set(inner_valid)

            if itr & iva:
                raise ValueError(
                    f"Outer {outer}, inner {inner}: train/valid overlap"
                )

            if itr | iva != train_set:
                raise ValueError(
                    f"Outer {outer}, inner {inner}: "
                    "inner train + valid do not equal outer train"
                )

            if itr & test_set or iva & test_set:
                raise ValueError(
                    f"Outer {outer}, inner {inner}: outer-test leakage"
                )

            observed = (len(inner_train), len(inner_valid))

            if observed != expected:
                raise ValueError(
                    f"Outer {outer}, inner {inner}: "
                    f"expected {expected}, found {observed}"
                )

            folds[(outer, inner)] = (inner_train, inner_valid)

    print("PASS: all 15 inner-fold partitions validated")

    # 2. Read and validate both complete phenotype files.
    all_values = {}

    for trait, (filename, column) in TRAITS.items():
        values = read_phenotypes(phenotype_root / filename, column)

        if set(values) != reference_ids:
            missing = reference_ids - set(values)
            extra = set(values) - reference_ids
            raise ValueError(
                f"{trait}: phenotype/sample ID mismatch; "
                f"missing={len(missing)}, extra={len(extra)}"
            )

        all_values[trait] = values
        print(f"PASS: {trait}: 904 unique, finite phenotypes")

    # 3. Prepare training and validation phenotypes.
    manifest = []

    for (outer, inner), (train_ids, valid_ids) in folds.items():
        outdir = output_root / f"outer_{outer}" / f"inner_{inner}"

        for trait, values in all_values.items():
            train_file = outdir / f"{trait}.train.pheno"
            valid_file = outdir / f"{trait}.valid.pheno"

            write_pheno(train_file, train_ids, values)
            write_pheno(valid_file, valid_ids, values)

            manifest.append({
                "Outer": outer,
                "Inner": inner,
                "Trait": trait,
                "N_train": len(train_ids),
                "N_valid": len(valid_ids),
                "Train_file": str(train_file.relative_to(work)),
                "Valid_file": str(valid_file.relative_to(work)),
            })

    # 4. Write a manifest for auditing.
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_file = output_root / "phenotype_manifest.tsv"
    temp = manifest_file.with_suffix(".tsv.tmp")

    with temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Outer", "Inner", "Trait",
                "N_train", "N_valid",
                "Train_file", "Valid_file",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest)

    os.replace(temp, manifest_file)

    if len(manifest) != 30:
        raise ValueError(f"Expected 30 manifest rows, got {len(manifest)}")

    print()
    print("OVERALL PASS")
    print("15 inner folds validated")
    print("30 trait/fold combinations prepared")
    print("60 phenotype files written")
    print(f"Manifest: {manifest_file}")
    print("No GWAS or prediction models were run.")


if __name__ == "__main__":
    main()
