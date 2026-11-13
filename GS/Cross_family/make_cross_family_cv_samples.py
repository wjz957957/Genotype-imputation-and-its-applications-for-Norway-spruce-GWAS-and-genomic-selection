#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

FOLDS = [f"{i:02d}" for i in range(1, 6)]
INNER_FOLDS = ["1", "2", "3"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header: {path}")
        return [
            {k: (v or "").strip() for k, v in row.items()}
            for row in reader
        ]


def read_sample_ids(path, has_header=False):
    ids = []
    with open(path, encoding="utf-8-sig") as f:
        for i, line in enumerate(f):
            fields = line.split()
            if not fields:
                continue
            if has_header and i == 0:
                continue
            if len(fields) < 2:
                raise ValueError(f"Malformed row in {path}: {line.rstrip()}")
            ids.append(fields[1])

    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate sample IDs: {path}")

    return set(ids)


def write_ids(path, ids):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for sid in sorted(ids):
            f.write(sid + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", default=".")
    parser.add_argument(
        "--outdir", default="Cross_family/cv_samples_v1"
    )
    args = parser.parse_args()

    work = Path(args.work).resolve()
    cross = work / "Cross_family"
    outdir = work / args.outdir

    assignment_file = (
        cross / "progeny_parent_connected_outer_fold_assignment.tsv"
    )
    nested_file = (
        cross / "nested_parent_connected_cv_assignments.tsv"
    )
    fam_file = (
        work / "08_imputed_GWAS/00_reference/"
        "all_904_genotypes.fam"
    )
    pheno_files = [
        work / "00_input/Trait-Adj_Budburst.txt",
        work / "00_input/Trait-Adj_H6.txt",
    ]

    for path in [
        assignment_file, nested_file, fam_file, *pheno_files
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"Missing input: {path}")

    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(
            f"Output directory is not empty: {outdir}\n"
            "Choose a new --outdir to avoid overwriting results."
        )

    assignment = read_tsv(assignment_file)
    nested = read_tsv(nested_file)

    required_assignment = {
        "Genotype_name", "Family", "Sample.Filename",
        "Parent_connected_group", "Outer_test_fold"
    }
    required_nested = {
        "Outer_fold", "Parent_connected_group", "Family",
        "Outer_role", "Inner_validation_fold"
    }

    if not required_assignment.issubset(
        assignment[0].keys() if assignment else set()
    ):
        raise ValueError("Unexpected progeny assignment columns")

    if not required_nested.issubset(
        nested[0].keys() if nested else set()
    ):
        raise ValueError("Unexpected nested assignment columns")

    if len(assignment) != 904:
        raise ValueError(
            f"Expected 904 progeny rows, got {len(assignment)}"
        )

    if len(nested) != 160:
        raise ValueError(
            f"Expected 160 nested records, got {len(nested)}"
        )

    # Use Sample.Filename: these are the genotype/phenotype IDs.
    samples = {}
    family_to_group = {}
    family_to_outer = {}
    samples_by_family = defaultdict(list)

    for row in assignment:
        sid = row["Sample.Filename"]
        family = row["Family"]
        group = row["Parent_connected_group"]
        outer = row["Outer_test_fold"]

        if not all([sid, family, group, outer]):
            raise ValueError(f"Incomplete assignment row: {row}")

        if sid in samples:
            raise ValueError(f"Duplicate sample ID: {sid}")

        if family in family_to_group:
            if family_to_group[family] != group:
                raise ValueError(
                    f"Family assigned to multiple groups: {family}"
                )
            if family_to_outer[family] != outer:
                raise ValueError(
                    f"Family split across outer folds: {family}"
                )

        samples[sid] = row
        family_to_group[family] = group
        family_to_outer[family] = outer
        samples_by_family[family].append(sid)

    families = set(family_to_group)

    if len(families) != 32:
        raise ValueError(f"Expected 32 families, got {len(families)}")

    # Check genotype and phenotype sample IDs.
    fam_ids = read_sample_ids(fam_file)

    if fam_ids != set(samples):
        raise ValueError(
            "Assignment IDs do not match PLINK FAM IDs: "
            f"missing={len(fam_ids - set(samples))}, "
            f"extra={len(set(samples) - fam_ids)}"
        )

    for pheno_file in pheno_files:
        pheno_ids = read_sample_ids(pheno_file, has_header=True)
        if pheno_ids != set(samples):
            raise ValueError(
                f"Phenotype IDs do not match assignment: {pheno_file}; "
                f"missing={len(set(samples) - pheno_ids)}, "
                f"extra={len(pheno_ids - set(samples))}"
            )

    # Index nested assignments by outer fold and family.
    nested_by_outer = defaultdict(dict)

    for row in nested:
        outer = row["Outer_fold"]
        family = row["Family"]

        if outer not in FOLDS:
            raise ValueError(f"Unexpected outer fold: {outer}")

        if family not in families:
            raise ValueError(f"Unknown family in nested file: {family}")

        if family in nested_by_outer[outer]:
            raise ValueError(
                f"Duplicate nested record: outer={outer}, family={family}"
            )

        nested_by_outer[outer][family] = row

    if set(nested_by_outer) != set(FOLDS):
        raise ValueError("Not all five outer folds are present")

    for outer in FOLDS:
        records = nested_by_outer[outer]

        if set(records) != families:
            raise ValueError(
                f"Outer {outer}: nested file does not contain "
                "every family exactly once"
            )

        # Verify outer and inner roles at family/group level.
        group_roles = defaultdict(set)
        group_inner_folds = defaultdict(set)

        for family, row in records.items():
            role = row["Outer_role"]
            group = row["Parent_connected_group"]
            inner = row["Inner_validation_fold"]

            if group != family_to_group[family]:
                raise ValueError(
                    f"Group mismatch: outer={outer}, family={family}"
                )

            expected_role = (
                "TEST" if family_to_outer[family] == outer else "TRAIN"
            )
            if role != expected_role:
                raise ValueError(
                    f"Outer role mismatch: outer={outer}, "
                    f"family={family}, observed={role}, "
                    f"expected={expected_role}"
                )

            group_roles[group].add(role)

            if role == "TEST":
                if inner:
                    raise ValueError(
                        f"Test family has inner fold: {family}"
                    )
            elif role == "TRAIN":
                if inner not in INNER_FOLDS:
                    raise ValueError(
                        f"Invalid inner fold for {family}: {inner}"
                    )
                group_inner_folds[group].add(inner)
            else:
                raise ValueError(f"Unknown outer role: {role}")

        for group, roles in group_roles.items():
            if len(roles) != 1:
                raise ValueError(
                    f"Parent-connected group split across outer roles: "
                    f"{outer}, {group}"
                )

        for group, inner_folds in group_inner_folds.items():
            if len(inner_folds) != 1:
                raise ValueError(
                    f"Parent-connected group split across inner folds: "
                    f"{outer}, {group}, {inner_folds}"
                )

        test_families = {
            fam for fam, row in records.items()
            if row["Outer_role"] == "TEST"
        }
        train_families = families - test_families

        outer_test_ids = {
            sid for sid, row in samples.items()
            if row["Outer_test_fold"] == outer
        }
        outer_train_ids = set(samples) - outer_test_ids

        ids_from_test_families = {
            sid for fam in test_families
            for sid in samples_by_family[fam]
        }
        ids_from_train_families = {
            sid for fam in train_families
            for sid in samples_by_family[fam]
        }

        if outer_test_ids != ids_from_test_families:
            raise ValueError(f"Outer {outer}: test ID mismatch")

        if outer_train_ids != ids_from_train_families:
            raise ValueError(f"Outer {outer}: train ID mismatch")

        outer_dir = outdir / f"outer_{outer}"
        write_ids(
            outer_dir / "outer_train.samples", outer_train_ids
        )
        write_ids(
            outer_dir / "outer_test.samples", outer_test_ids
        )

        print(
            f"Outer {outer}: train={len(outer_train_ids)}, "
            f"test={len(outer_test_ids)}, "
            f"train_families={len(train_families)}, "
            f"test_families={len(test_families)}"
        )

        for inner in INNER_FOLDS:
            valid_families = {
                fam for fam in train_families
                if records[fam]["Inner_validation_fold"] == inner
            }
            inner_train_families = train_families - valid_families

            valid_ids = {
                sid for fam in valid_families
                for sid in samples_by_family[fam]
            }
            inner_train_ids = {
                sid for fam in inner_train_families
                for sid in samples_by_family[fam]
            }

            if inner_train_ids & valid_ids:
                raise ValueError(
                    f"Outer {outer}, inner {inner}: train/valid overlap"
                )

            if inner_train_ids | valid_ids != outer_train_ids:
                raise ValueError(
                    f"Outer {outer}, inner {inner}: "
                    "inner train + validation != outer training set"
                )

            if (inner_train_ids | valid_ids) & outer_test_ids:
                raise ValueError(
                    f"Outer {outer}, inner {inner}: outer-test leakage"
                )

            inner_dir = outer_dir / f"inner_{int(inner):02d}"

            write_ids(
                inner_dir / "inner_train.samples", inner_train_ids
            )
            write_ids(
                inner_dir / "inner_valid.samples", valid_ids
            )

            print(
                f"  Inner {inner}: train={len(inner_train_ids)}, "
                f"valid={len(valid_ids)}"
            )

    print("\nPASS: all sample IDs and nested group partitions validated.")
    print(f"Sample lists saved under: {outdir}")
    print("No GWAS or prediction models were run.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
