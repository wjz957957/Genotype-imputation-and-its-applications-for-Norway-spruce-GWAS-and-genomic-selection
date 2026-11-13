#!/usr/bin/env python3
"""
Preflight validation for family/pedigree-grouped nested CV inputs.

Run from the GS-904/GS-reviewer working directory:
  python3 Cross_family/check_cross_family_inputs.py \
    --work "$PWD" --cross-dir Cross_family --pheno-dir 00_input

Expected phenotype format: whitespace-delimited, header + exactly 3 columns:
FID IID phenotype_value
FID must equal IID. The script compares phenotype IID with Genotype_name
in progeny_parent_connected_outer_fold_assignment.tsv.
"""
import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np

def fail(msg):
    print(f"[FAIL] {msg}")
    raise SystemExit(1)

def passmsg(msg):
    print(f"[PASS] {msg}")

def read_tsv(path, required):
    if not path.is_file():
        fail(f"Missing file: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in required if c not in df.columns]
    if missing:
        fail(f"{path.name} missing required columns: {missing}; found {list(df.columns)}")
    return df

def norm_fold(s):
    s = str(s).strip()
    try:
        return f"{int(float(s)):02d}"
    except Exception:
        return s.zfill(2)

def read_pheno(path):
    if not path.is_file():
        fail(f"Missing phenotype file: {path}")
    try:
        df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)
    except Exception as e:
        fail(f"Cannot read phenotype file {path}: {e}")
    if df.shape[1] != 3:
        fail(f"{path} must have exactly 3 columns (FID IID value); got {df.shape[1]}: {list(df.columns)}")
    df.columns = ["FID", "IID", "VALUE"]
    for c in ["FID", "IID"]:
        df[c] = df[c].astype(str).str.strip()
    if (df["FID"] == "").any() or (df["IID"] == "").any():
        fail(f"{path} contains blank FID/IID")
    if df["IID"].duplicated().any():
        dups = df.loc[df["IID"].duplicated(keep=False), "IID"].head(10).tolist()
        fail(f"{path} has duplicate IID(s), examples: {dups}")
    if not (df["FID"] == df["IID"]).all():
        n = int((df["FID"] != df["IID"]).sum())
        ex = df.loc[df["FID"] != df["IID"], ["FID", "IID"]].head(5).to_dict("records")
        fail(f"{path} has {n} rows where FID != IID; examples: {ex}")
    vals = pd.to_numeric(df["VALUE"], errors="coerce")
    if vals.isna().any() or not np.isfinite(vals.to_numpy(dtype=float)).all():
        fail(f"{path} contains missing/non-numeric/non-finite phenotype values")
    df["VALUE"] = vals
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=".", help="GS-904/GS-reviewer working directory")
    ap.add_argument("--cross-dir", default="Cross_family")
    ap.add_argument("--pheno-dir", default="00_input")
    args = ap.parse_args()

    work = Path(args.work).resolve()
    cross = (work / args.cross_dir).resolve()
    pheno = (work / args.pheno_dir).resolve()

    assignment = read_tsv(
        cross / "progeny_parent_connected_outer_fold_assignment.tsv",
        ["Genotype_id", "Genotype_name", "Family", "Sample.Filename",
         "Parent_connected_group", "Outer_test_fold"]
    )
    famfold = read_tsv(
        cross / "parent_connected_outer_fold_assignment.tsv",
        ["Family", "N_progeny", "Parent_connected_group", "Outer_test_fold",
         "Shared_parent_ids_in_group"]
    )
    groups = read_tsv(
        cross / "parent_connected_group_assignment.tsv",
        ["Parent_connected_group", "N_families", "N_progeny", "Outer_test_fold",
         "Families", "Known_parents"]
    )
    nested = read_tsv(
        cross / "nested_parent_connected_cv_assignments.tsv",
        ["Outer_fold", "Parent_connected_group", "Family", "Outer_role", "Inner_validation_fold"]
    )

    # Normalize IDs and fold labels.
    for df, cols in [
        (assignment, ["Genotype_name", "Family", "Parent_connected_group", "Sample.Filename"]),
        (famfold, ["Family", "Parent_connected_group", "Shared_parent_ids_in_group"]),
        (groups, ["Parent_connected_group", "Families", "Known_parents"]),
        (nested, ["Parent_connected_group", "Family", "Outer_role", "Inner_validation_fold"])
    ]:
        for c in cols:
            df[c] = df[c].astype(str).str.strip()
    assignment["Outer_test_fold"] = assignment["Outer_test_fold"].map(norm_fold)
    famfold["Outer_test_fold"] = famfold["Outer_test_fold"].map(norm_fold)
    groups["Outer_test_fold"] = groups["Outer_test_fold"].map(norm_fold)
    nested["Outer_fold"] = nested["Outer_fold"].map(norm_fold)
    nested["Inner_validation_fold"] = nested["Inner_validation_fold"].replace("", np.nan)

    # Assignment integrity.
    if len(assignment) != 904:
        fail(f"Expected 904 progeny assignment rows, found {len(assignment)}")
    if assignment["Genotype_name"].duplicated().any():
        fail("Genotype_name is duplicated in progeny assignment")
    if assignment["Sample.Filename"].duplicated().any():
        fail("Sample.Filename is duplicated in progeny assignment")
    if assignment["Family"].nunique() != 32:
        fail(f"Expected 32 Family labels, found {assignment['Family'].nunique()}")
    if assignment["Parent_connected_group"].nunique() != 17:
        fail(f"Expected 17 parent-connected groups, found {assignment['Parent_connected_group'].nunique()}")
    if assignment["Outer_test_fold"].nunique() != 5:
        fail(f"Expected 5 outer folds, found {assignment['Outer_test_fold'].nunique()}")
    if assignment.groupby("Family")["Outer_test_fold"].nunique().max() != 1:
        fail("At least one Family label is split across outer folds")
    if assignment.groupby("Parent_connected_group")["Outer_test_fold"].nunique().max() != 1:
        fail("At least one parent-connected group is split across outer folds")
    if set(assignment["Family"]) != set(famfold["Family"]):
        fail("Family labels differ between progeny assignment and family-fold file")
    map_a = assignment.drop_duplicates("Family").set_index("Family")["Parent_connected_group"].to_dict()
    map_f = famfold.set_index("Family")["Parent_connected_group"].to_dict()
    if any(map_a.get(f) != map_f.get(f) for f in map_f):
        fail("Family-to-parent-connected-group mapping differs between assignment files")
    fold_map_a = assignment.drop_duplicates("Family").set_index("Family")["Outer_test_fold"].to_dict()
    fold_map_f = famfold.set_index("Family")["Outer_test_fold"].to_dict()
    if any(fold_map_a.get(f) != fold_map_f.get(f) for f in fold_map_f):
        fail("Outer-fold assignments differ between assignment files")
    passmsg("904 progeny, 32 Family labels, 17 parent-connected groups and 5 outer folds are structurally consistent")

    # Parent IDs should not cross outer folds.
    parent_folds = {}
    for _, row in famfold.iterrows():
        fold = row["Outer_test_fold"]
        raw = row["Shared_parent_ids_in_group"]
        for parent in [x.strip() for x in raw.split(";") if x.strip()]:
            parent_folds.setdefault(parent, set()).add(fold)
    crossfold_parents = {p: fs for p, fs in parent_folds.items() if len(fs) > 1}
    if crossfold_parents:
        fail(f"Recorded parent(s) appear in more than one outer fold: {list(crossfold_parents.items())[:10]}")
    passmsg("No recorded shared parent crosses outer test folds")

    # Outer-fold role and inner-fold integrity.
    for outer in sorted(nested["Outer_fold"].unique()):
        sub = nested[nested["Outer_fold"] == outer].copy()
        # Every group must have one role and, if TRAIN, one inner fold.
        for grp, g in sub.groupby("Parent_connected_group"):
            if g["Outer_role"].nunique() != 1:
                fail(f"Outer {outer}: group {grp} has inconsistent outer roles")
            role = g["Outer_role"].iloc[0]
            fam = g["Family"].iloc[0]
            expected_fold = fold_map_a[fam]
            if role == "TEST":
                if expected_fold != outer:
                    fail(f"Outer {outer}: TEST group {grp} is assigned to outer fold {expected_fold}")
                if g["Inner_validation_fold"].notna().any():
                    fail(f"Outer {outer}: TEST group {grp} incorrectly has an inner fold")
            elif role == "TRAIN":
                if expected_fold == outer:
                    fail(f"Outer {outer}: group {grp} should be TEST, not TRAIN")
                vals = g["Inner_validation_fold"].dropna().astype(str).unique()
                if len(vals) != 1 or vals[0] not in {"1", "2", "3", "1.0", "2.0", "3.0"}:
                    fail(f"Outer {outer}: TRAIN group {grp} has invalid/inconsistent inner fold {vals}")
            else:
                fail(f"Outer {outer}: unknown Outer_role {role}")
        # Test groups should exactly match family outer-fold assignments.
        expected_test = {g for g, f in assignment.drop_duplicates("Parent_connected_group")[["Parent_connected_group","Outer_test_fold"]].itertuples(index=False, name=None) if f == outer}
        observed_test = set(sub.loc[sub["Outer_role"] == "TEST", "Parent_connected_group"])
        if expected_test != observed_test:
            fail(f"Outer {outer}: test group set does not match outer assignment")
    passmsg("Nested outer/inner group roles are consistent; test groups are held out and training groups have one inner fold")

    # Check family counts in group table.
    if groups["Parent_connected_group"].nunique() != 17:
        fail("Parent-connected group table does not contain 17 unique groups")
    if set(groups["Parent_connected_group"]) != set(assignment["Parent_connected_group"]):
        fail("Parent-connected group labels differ between group table and progeny assignment")
    fam_counts = assignment.groupby("Family").size().to_dict()
    for _, row in famfold.iterrows():
        fam = row["Family"]
        try:
            expected_n = int(float(row["N_progeny"]))
        except Exception:
            fail(f"Invalid N_progeny for family {fam}")
        if fam_counts.get(fam) != expected_n:
            fail(f"Family {fam} has {fam_counts.get(fam)} assignment rows but N_progeny says {expected_n}")
    if int(pd.to_numeric(groups["N_progeny"]).sum()) != 904:
        fail("Parent-connected group N_progeny values do not sum to 904")
    passmsg("Family and group sample counts match the 904 assigned progeny")

    # Validate both actual server phenotype files against Genotype_name / IID.
    expected_ids = set(assignment["Sample.Filename"])
    trait_files = [
        ("BB", pheno / "Trait-Adj_Budburst.txt"),
        ("H6", pheno / "Trait-Adj_H6.txt"),
    ]
    for trait, path in trait_files:
        p = read_pheno(path)
        observed = set(p["IID"])
        missing = sorted(expected_ids - observed)
        extra = sorted(observed - expected_ids)
        if missing or extra:
            fail(
                f"{trait} phenotype IDs do not exactly match the 904 Genotype_name IDs: "
                f"{len(missing)} expected IDs missing, {len(extra)} unexpected IDs; "
                f"missing examples={missing[:10]}, extra examples={extra[:10]}"
            )
        print(f"[PASS] {trait}: {len(p)} unique phenotype IDs; exact match to the 904 assigned progeny")
    passmsg("Both trait phenotype files match the family-fold assignment exactly")
    print("\nOVERALL PASS: Cross-family assignment and phenotype-ID prechecks succeeded.")
    print("Note: this script checks assignment/phenotype consistency; genotype matrix/raw-genotype ID matching must be checked separately before model fitting.")

if __name__ == "__main__":
    main()
