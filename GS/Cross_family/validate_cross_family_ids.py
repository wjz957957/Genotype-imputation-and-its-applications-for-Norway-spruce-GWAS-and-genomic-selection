#!/usr/bin/env python3
import csv
import sys
from pathlib import Path
from collections import Counter, defaultdict

WORK = Path.cwd()
FAM = WORK / "08_imputed_GWAS/00_reference/all_904_genotypes.fam"
ASSIGN = WORK / "Cross_family/progeny_parent_connected_outer_fold_assignment.tsv"
NESTED = WORK / "Cross_family/nested_parent_connected_cv_assignments.tsv"
BB = WORK / "00_input/Trait-Adj_Budburst.txt"
H6 = WORK / "00_input/Trait-Adj_H6.txt"

errors = []

def check(condition, message):
    if condition:
        print(f"[PASS] {message}")
    else:
        print(f"[FAIL] {message}")
        errors.append(message)

def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = [
            {k.strip(): (v.strip() if v is not None else "")
             for k, v in row.items()}
            for row in reader
        ]
        return reader.fieldnames, rows

def read_space_table(path):
    with open(path, encoding="utf-8-sig") as f:
        lines = [line.split() for line in f if line.strip()]
    if not lines:
        raise ValueError(f"Empty file: {path}")
    header = lines[0]
    rows = []
    for i, fields in enumerate(lines[1:], start=2):
        if len(fields) != len(header):
            raise ValueError(
                f"{path}: line {i} has {len(fields)} fields; "
                f"expected {len(header)}"
            )
        rows.append(dict(zip(header, fields)))
    return header, rows

def nonblank(value):
    return str(value).strip().lower() not in {
        "", "na", "nan", "none"
    }

# --------------------------------------------------
# 1. Read all files
# --------------------------------------------------
for p in [FAM, ASSIGN, NESTED, BB, H6]:
    if not p.is_file():
        print(f"[FAIL] Missing file: {p}")
        sys.exit(1)

with open(FAM, encoding="utf-8-sig") as f:
    fam_rows = [line.split() for line in f if line.strip()]

assign_header, assign = read_tsv(ASSIGN)
nested_header, nested = read_tsv(NESTED)
bb_header, bb = read_space_table(BB)
h6_header, h6 = read_space_table(H6)

print("===== 1. INPUT COUNTS =====")
print(f"PLINK FAM rows: {len(fam_rows)}")
print(f"Family assignment rows: {len(assign)}")
print(f"Nested assignment rows: {len(nested)}")
print(f"BB phenotype rows: {len(bb)}")
print(f"H6 phenotype rows: {len(h6)}")

# --------------------------------------------------
# 2. Check PLINK and phenotype IDs
# --------------------------------------------------
check(all(len(x) >= 2 for x in fam_rows),
      "Every PLINK FAM row has at least FID and IID")

fam_fid = [x[0] for x in fam_rows if len(x) >= 2]
fam_iid = [x[1] for x in fam_rows if len(x) >= 2]

check(len(fam_iid) == 904, "PLINK FAM contains 904 individuals")
check(len(set(fam_iid)) == len(fam_iid),
      "PLINK IID values are unique")
check(fam_fid == fam_iid,
      "PLINK FID and IID are identical in every row")

sample_ids = set(fam_iid)

assignment_ids = [r["Sample.Filename"] for r in assign]
check(len(assign) == 904,
      "Family assignment contains 904 individuals")
check(len(set(assignment_ids)) == len(assignment_ids),
      "Sample.Filename is unique in the family assignment")

check(sample_ids == set(assignment_ids),
      "PLINK sample IDs exactly match family-assignment IDs")

for trait, rows, expected_col in [
    ("BB", bb, "Adj_Budburst"),
    ("H6", h6, "Adj_H6"),
]:
    check({"FID", "IID", expected_col}.issubset(rows[0].keys())
          if rows else False,
          f"{trait} phenotype contains expected columns")

    ids = [r["IID"] for r in rows]
    fids = [r["FID"] for r in rows]

    check(len(rows) == 904,
          f"{trait} contains 904 phenotype records")
    check(len(set(ids)) == len(ids),
          f"{trait} IID values are unique")
    check(fids == ids,
          f"{trait} FID and IID are identical")
    check(set(ids) == sample_ids,
          f"{trait} IDs exactly match PLINK IDs")

# --------------------------------------------------
# 3. Check family/group/fold assignments
# --------------------------------------------------
family_sizes = Counter(r["Family"] for r in assign)
family_group = {}
family_fold = {}
group_fold = {}
group_families = defaultdict(set)

for r in assign:
    fam = r["Family"]
    group = r["Parent_connected_group"]
    fold = r["Outer_test_fold"]

    if fam in family_group and family_group[fam] != group:
        errors.append(f"Family {fam} occurs in multiple groups")
    family_group[fam] = group

    if fam in family_fold and family_fold[fam] != fold:
        errors.append(f"Family {fam} occurs in multiple outer folds")
    family_fold[fam] = fold

    if group in group_fold and group_fold[group] != fold:
        errors.append(f"Group {group} occurs in multiple outer folds")
    group_fold[group] = fold
    group_families[group].add(fam)

folds = sorted(set(family_fold.values()))
families = set(family_group)
groups = set(group_fold)

print("\n===== 2. OUTER-FOLD STRUCTURE =====")
print(f"Families: {len(families)}")
print(f"Parent-connected groups: {len(groups)}")
print(f"Outer folds: {folds}")

check(len(families) == 32, "Exactly 32 families")
check(len(groups) == 17, "Exactly 17 parent-connected groups")
check(len(folds) == 5, "Exactly 5 outer folds")
check(not errors,
      "No family or parent-connected group is split across outer folds")

# --------------------------------------------------
# 4. Check nested assignment
# --------------------------------------------------
nested_by_key = {}
for r in nested:
    key = (r["Outer_fold"], r["Family"])
    if key in nested_by_key:
        errors.append(f"Duplicate nested assignment: {key}")
    nested_by_key[key] = r

expected_keys = {
    (fold, fam)
    for fold in folds
    for fam in families
}

check(len(nested) == 5 * 32,
      "Nested assignment contains 160 records")
check(set(nested_by_key) == expected_keys,
      "Every family has exactly one record in each outer fold")

nested_group_state = defaultdict(set)

for (fold, fam), r in nested_by_key.items():
    if fam not in family_group:
        errors.append(f"Unknown family in nested assignment: {fam}")
        continue

    expected_group = family_group[fam]
    expected_outer_fold = family_fold[fam]
    expected_role = (
        "TEST" if fold == expected_outer_fold else "TRAIN"
    )

    if r["Parent_connected_group"] != expected_group:
        errors.append(f"Group mismatch: fold={fold}, family={fam}")

    if r["Outer_role"] != expected_role:
        errors.append(f"Outer role mismatch: fold={fold}, family={fam}")

    inner = r["Inner_validation_fold"]
    if expected_role == "TEST":
        if nonblank(inner):
            errors.append(
                f"Outer-test family has inner-fold assignment: {fold}, {fam}"
            )
    elif inner not in {"1", "2", "3"}:
        errors.append(
            f"Invalid inner fold: outer={fold}, family={fam}, inner={inner}"
        )

    nested_group_state[(fold, expected_group)].add(
        (r["Outer_role"], inner if nonblank(inner) else "")
    )

for key, states in nested_group_state.items():
    if len(states) != 1:
        errors.append(
            f"Parent-connected group has inconsistent nested roles: {key}"
        )

# --------------------------------------------------
# 5. Verify sample counts and inner-fold partitioning
# --------------------------------------------------
print("\n===== 3. SAMPLE COUNTS BY OUTER/INNER FOLD =====")

for fold in folds:
    test_families = {
        fam for fam in families if family_fold[fam] == fold
    }
    train_families = families - test_families

    test_n = sum(family_sizes[f] for f in test_families)
    train_n = sum(family_sizes[f] for f in train_families)

    inner_families = {}
    for k in ("1", "2", "3"):
        inner_families[k] = {
            fam for fam in train_families
            if nested_by_key.get((fold, fam), {}).get(
                "Inner_validation_fold", ""
            ) == k
        }

    inner_counts = {
        k: sum(family_sizes[f] for f in fs)
        for k, fs in inner_families.items()
    }

    inner_union = set().union(*inner_families.values())
    inner_disjoint = (
        sum(len(fs) for fs in inner_families.values())
        == len(inner_union)
    )

    check(inner_union == train_families and inner_disjoint,
          f"Outer {fold}: inner validation folds partition training families")

    print(
        f"Outer {fold}: test={test_n}, train={train_n}; "
        f"inner validation counts={inner_counts}"
    )

# --------------------------------------------------
# 6. Final result
# --------------------------------------------------
print("\n===== 4. FINAL STATUS =====")

if errors:
    print(f"VALIDATION FAILED: {len(errors)} issue(s)")
    for e in errors:
        print(" -", e)
    sys.exit(1)

print("OVERALL PASS")
print("All sample IDs, family assignments, outer folds,")
print("and nested inner-fold roles passed validation.")
print("No GWAS or prediction models were run.")
