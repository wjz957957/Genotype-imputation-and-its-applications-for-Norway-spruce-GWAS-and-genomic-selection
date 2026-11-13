#!/usr/bin/env python3
"""Verify the 5-fold matrix-baseline outputs before assembling seven-model results."""
import argparse
import csv
import math
from pathlib import Path

TRAITS = ("BB", "H6")
MODELS = ("ABLUP", "RRBLUP-array", "RRBLUP-imputed")
EXPECTED_TEST_N = {1: 175, 2: 178, 3: 176, 4: 185, 5: 190}


def read_tsv(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing/empty file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_flag(path):
    if not path.is_file():
        raise ValueError(f"Missing flag: {path}")
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            raise ValueError(f"Malformed flag line in {path}: {line}")
        key, value = line.split("=", 1)
        if key in result:
            raise ValueError(f"Duplicate key {key} in {path}")
        result[key] = value
    return result


def number_or_nan(text):
    return float("nan") if str(text).strip().upper() in ("NA", "NAN", "") else float(text)


def corr(x, y):
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    vx = sum((v - mx) ** 2 for v in x)
    vy = sum((v - my) ** 2 for v in y)
    if vx <= 0 or vy <= 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(vx * vy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=".")
    ap.add_argument("--out", default="Cross_family/results/20_outer_matrix_baselines")
    args = ap.parse_args()
    work = Path(args.work).resolve()
    out = work / args.out
    assignment_path = work / "Cross_family/config/progeny_parent_connected_outer_fold_assignment.tsv"
    assignment = read_tsv(assignment_path)
    if len(assignment) != 904:
        raise ValueError(f"Expected 904 assignment rows, found {len(assignment)}")
    expected_ids = {}
    for outer in range(1, 6):
        expected_ids[outer] = {
            row["Sample.Filename"] for row in assignment
            if int(row["Outer_test_fold"]) == outer
        }
        if len(expected_ids[outer]) != EXPECTED_TEST_N[outer]:
            raise ValueError(f"Fold {outer:02d} assignment count mismatch")

    metrics_seen = set()
    npass = 0
    for trait in TRAITS:
        for outer in range(1, 6):
            d = out / "per_fold" / trait / f"outer_{outer:02d}"
            for model in MODELS:
                key = (trait, outer, model)
                if key in metrics_seen:
                    raise ValueError(f"Duplicate task key: {key}")
                flag = read_flag(d / f"{model}.complete.flag")
                if not (flag.get("status") == "PASS" and flag.get("outer") == f"{outer:02d}"
                        and flag.get("trait") == trait and flag.get("model") == model
                        and flag.get("test_phenotype_read_after_prediction") == "YES"):
                    raise ValueError(f"Completion flag mismatch: {d}/{model}.complete.flag")
                metric_rows = read_tsv(d / f"{model}_metrics.tsv")
                pred_rows = read_tsv(d / f"{model}_predictions.tsv")
                if len(metric_rows) != 1:
                    raise ValueError(f"Expected one metric row: {d}/{model}_metrics.tsv")
                if len(pred_rows) != EXPECTED_TEST_N[outer]:
                    raise ValueError(f"Prediction count mismatch: {key}, rows={len(pred_rows)}")
                pred_ids = [r["ID"] for r in pred_rows]
                if len(set(pred_ids)) != len(pred_ids) or set(pred_ids) != expected_ids[outer]:
                    raise ValueError(f"Prediction IDs differ from held-out fold: {key}")
                obs = [float(r["Observed_raw"]) for r in pred_rows]
                pred = [float(r["Predicted_raw"]) for r in pred_rows]
                if not all(math.isfinite(v) for v in obs + pred):
                    raise ValueError(f"Non-finite prediction/observation: {key}")
                met = metric_rows[0]
                r = number_or_nan(met["Pearson"])
                rho = number_or_nan(met["Spearman"])
                recomputed = corr(obs, pred)
                status = met.get("Prediction_status")
                flag_status = flag.get("prediction_status")
                if status != flag_status:
                    raise ValueError(f"Prediction status differs between metrics and flag: {key}")
                cross = float(met["Max_abs_K_test_train"])
                pred_range = max(pred) - min(pred)
                if abs(pred_range - float(met["Prediction_range"])) > 1e-7:
                    raise ValueError(f"Prediction range mismatch: {key}")
                if status == "UNDEFINED_NO_PEDIGREE_LINK":
                    if model != "ABLUP" or cross > 1e-12 or pred_range > 1e-10:
                        raise ValueError(f"Undefined prediction correlation not explained by disconnected pedigree: {key}")
                    if not math.isnan(r) or not math.isnan(rho):
                        raise ValueError(f"ABLUP undefined correlation must be NA, not zero: {key}")
                    if not math.isnan(recomputed) and pred_range > 1e-12:
                        raise ValueError(f"Unexpected positive prediction variation: {key}")
                elif status == "DEFINED":
                    if not math.isfinite(r) or not math.isfinite(rho) or not math.isfinite(recomputed) or abs(r - recomputed) > 1e-7:
                        raise ValueError(f"Pearson mismatch: {key}, file={r}, recomputed={recomputed}")
                    if pred_range <= 1e-10:
                        raise ValueError(f"Defined status with constant predictions: {key}")
                else:
                    raise ValueError(f"Unknown prediction status: {key}: {status}")
                for col in ("RMSE_raw", "MAE_raw", "Bias_raw"):
                    if not math.isfinite(float(met[col])):
                        raise ValueError(f"Non-finite {col}: {key}")
                if int(met["N_Train"]) != 904 - EXPECTED_TEST_N[outer] or int(met["N_Test"]) != EXPECTED_TEST_N[outer]:
                    raise ValueError(f"Sample count mismatch in metrics: {key}")
                metrics_seen.add(key)
                npass += 1
                print(f"PASS outer_{outer:02d}/{trait}/{model}: n_test={len(pred_rows)}, "
                      f"Pearson={'NA' if math.isnan(r) else f'{r:.6f}'}, status={status}")

    expected = {(t, o, m) for t in TRAITS for o in range(1, 6) for m in MODELS}
    if metrics_seen != expected or npass != 30:
        raise ValueError(f"Expected 30 unique tasks, found {len(metrics_seen)}")

    joined = read_tsv(out / "seven_models_baseline_outer5.tsv")
    if len(joined) != 30:
        raise ValueError(f"Expected 30 rows in join file, found {len(joined)}")
    join_keys = [(r["Trait"], int(r["Outer"]), r["Model"]) for r in joined]
    if len(set(join_keys)) != 30 or set(join_keys) != expected:
        raise ValueError("Join-file keys are missing, duplicated, or unexpected")
    for row in joined:
        for col in ("RMSE_original", "MAE_original", "RMSE_trainSD", "MAE_trainSD"):
            if not math.isfinite(float(row[col])):
                raise ValueError(f"Non-finite {col} in join file: {row}")
        status = row.get("Prediction_status")
        pearson, spearman = number_or_nan(row["Pearson"]), number_or_nan(row["Spearman"])
        if status == "UNDEFINED_NO_PEDIGREE_LINK":
            if row["Model"] != "ABLUP" or not (math.isnan(pearson) and math.isnan(spearman)):
                raise ValueError(f"Invalid undefined ABLUP join values: {row}")
        elif status == "DEFINED":
            if not (math.isfinite(pearson) and math.isfinite(spearman)):
                raise ValueError(f"Invalid defined join values: {row}")
        else:
            raise ValueError(f"Unknown joined prediction status: {row}")
    flag = read_flag(out / "complete.flag")
    if flag.get("status") != "PASS" or flag.get("n_metric_records") != "30":
        raise ValueError("Aggregate completion flag is not PASS/30 records")
    print("OVERALL PASS=30/30 model tasks")
    print("NOTE: ABLUP Pearson/Spearman=NA where held-out pedigree groups have zero cross-relatedness.")
    print("Undefined correlation must remain NA (not zero); exclude it from Pearson comparisons.")
    print(f"Seven-model baseline rows ready: {out / 'seven_models_baseline_outer5.tsv'}")


if __name__ == "__main__":
    main()
