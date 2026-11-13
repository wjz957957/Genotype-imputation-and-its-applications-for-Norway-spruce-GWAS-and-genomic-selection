#!/usr/bin/env python3
import csv
import math
import statistics
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

ROOT = Path("Cross_family/results/05_inner_rrblup")
OUT = Path("Cross_family/results/06_inner_rrblup_summary")

TRAITS = ["BB", "H6"]
MARKER_SETS = [
    "chip_50", "chip_100", "chip_200", "chip_300",
    "chip_400", "chip_500", "chip_1000", "chip_2000",
    "chip_3000", "chip_4000", "chip_5000", "chip_8000",
    "chip_10000", "chip_20000", "chip_30000", "chip_40000",
    "chip_LD_ALL", "chip_ALL",
]
OUTER_TRAIN_N = {
    "01": 729, "02": 726, "03": 728, "04": 719, "05": 714
}
EXPECTED_KEYS = set(product(TRAITS, MARKER_SETS))


def stop(message):
    raise SystemExit("ERROR: " + message)


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fields, delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def mean_finite(values):
    values = [x for x in values if x is not None]
    return statistics.mean(values) if values else None


def pearson(y, p):
    if len(y) < 2:
        return None
    ym = statistics.mean(y)
    pm = statistics.mean(p)
    a = sum((x - ym) ** 2 for x in y)
    b = sum((x - pm) ** 2 for x in p)
    if a <= 0 or b <= 0:
        return None
    return sum((x - ym) * (z - pm) for x, z in zip(y, p)) / math.sqrt(a * b)


OUT.mkdir(parents=True, exist_ok=True)

all_metrics = []
all_predictions = []

# 1. Validate every fold and read the outputs.
for oi in range(1, 6):
    outer = f"{oi:02d}"

    for ii in range(1, 4):
        inner = f"{ii:02d}"
        fold = ROOT / f"outer_{outer}" / f"inner_{inner}"

        required_files = [
            fold / "complete.flag",
            fold / "metrics.tsv",
            fold / "validation_predictions.tsv",
        ]
        for path in required_files:
            if not path.is_file() or path.stat().st_size == 0:
                stop(f"Missing/empty file: {path}")

        flags = {}
        for line in (fold / "complete.flag").read_text().splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                flags[k] = v

        required_flags = {
            "version": "1",
            "outer": outer,
            "inner": inner,
            "traits": "BB,H6",
            "marker_sets_per_trait": "18",
            "models": "36",
            "model_status": "PASS",
            "sample_ID_check": "PASS",
            "marker_ID_check": "PASS",
            "training_only_imputation_and_centering": "PASS",
            "validation_prediction_check": "PASS",
        }
        for k, v in required_flags.items():
            if flags.get(k) != v:
                stop(f"{fold}: flag {k} expected {v}, got {flags.get(k)}")

        metrics = read_tsv(fold / "metrics.tsv")
        predictions = read_tsv(fold / "validation_predictions.tsv")

        if len(metrics) != 36:
            stop(f"{fold}: expected 36 metrics rows, got {len(metrics)}")

        keys = [(r["trait"], r["marker_set"]) for r in metrics]
        if len(set(keys)) != 36 or set(keys) != EXPECTED_KEYS:
            stop(f"{fold}: marker-set/trait combinations are incorrect")

        if any(r.get("status") != "PASS" for r in metrics):
            stop(f"{fold}: at least one model is not PASS")

        pred_groups = defaultdict(list)
        for r in predictions:
            if r["outer"] != outer or r["inner"] != inner:
                stop(f"{fold}: incorrect fold ID in prediction file")
            pred_groups[(r["trait"], r["marker_set"])].append(r)

        expected_rows = 0
        for m in metrics:
            key = (m["trait"], m["marker_set"])
            rows = pred_groups[key]
            n_valid = int(m["n_valid"])
            expected_rows += n_valid

            if len(rows) != n_valid:
                stop(f"{fold}: prediction count mismatch for {key}")

            ids = [r["sample_id"] for r in rows]
            if len(ids) != len(set(ids)):
                stop(f"{fold}: duplicate validation IDs for {key}")

            for r in rows:
                if number(r["observed"]) is None:
                    stop(f"{fold}: invalid observed value for {key}")
                if number(r["predicted"]) is None:
                    stop(f"{fold}: invalid predicted value for {key}")

        if len(predictions) != expected_rows:
            stop(f"{fold}: incorrect total prediction count")

        all_metrics.extend(metrics)
        all_predictions.extend(predictions)
        print(f"PASS outer_{outer}/inner_{inner}")

# 2. Validate global counts.
if len(all_metrics) != 540:
    stop(f"Expected 540 metric rows, got {len(all_metrics)}")

expected_predictions = 36 * sum(OUTER_TRAIN_N.values())
if len(all_predictions) != expected_predictions:
    stop(
        f"Expected {expected_predictions} prediction rows, "
        f"got {len(all_predictions)}"
    )

metric_groups = defaultdict(list)
prediction_groups = defaultdict(list)

for r in all_metrics:
    metric_groups[(r["outer"], r["trait"], r["marker_set"])].append(r)

for r in all_predictions:
    prediction_groups[(r["outer"], r["trait"], r["marker_set"])].append(r)

# 3. Pool the three inner-validation folds for each outer/trait/set.
pooled_rows = []

for outer in [f"{i:02d}" for i in range(1, 6)]:
    for trait, marker_set in sorted(EXPECTED_KEYS):
        key = (outer, trait, marker_set)
        preds = prediction_groups[key]
        metrics = metric_groups[key]

        expected_n = OUTER_TRAIN_N[outer]
        ids = [r["sample_id"] for r in preds]

        if len(preds) != expected_n or len(set(ids)) != expected_n:
            stop(
                f"{key}: OOF sample coverage is {len(set(ids))}/"
                f"{expected_n}, rows={len(preds)}"
            )

        if len(metrics) != 3:
            stop(f"{key}: expected metrics from 3 inner folds")

        y = [float(r["observed"]) for r in preds]
        p = [float(r["predicted"]) for r in preds]

        pooled_r = pearson(y, p)
        pooled_rmse = math.sqrt(
            sum((a - b) ** 2 for a, b in zip(p, y)) / len(y)
        )
        pooled_mae = sum(abs(a - b) for a, b in zip(p, y)) / len(y)
        pooled_bias = sum(a - b for a, b in zip(p, y)) / len(y)

        inner_r = [number(r["pearson_r"]) for r in metrics]
        inner_rmse = [number(r["RMSE"]) for r in metrics]
        inner_mae = [number(r["MAE"]) for r in metrics]
        inner_bias = [number(r["bias"]) for r in metrics]
        n_used = [int(r["n_markers_used"]) for r in metrics]

        pooled_rows.append({
            "outer": outer,
            "trait": trait,
            "marker_set": marker_set,
            "n_inner_folds": 3,
            "n_oof": len(preds),
            "n_markers_requested": metrics[0]["n_markers_requested"],
            "mean_n_markers_used": statistics.mean(n_used),
            "min_n_markers_used": min(n_used),
            "max_n_markers_used": max(n_used),
            "pooled_pearson_r": pooled_r,
            "pooled_RMSE": pooled_rmse,
            "pooled_MAE": pooled_mae,
            "pooled_bias": pooled_bias,
            "mean_inner_pearson_r": mean_finite(inner_r),
            "sd_inner_pearson_r": (
                statistics.stdev([x for x in inner_r if x is not None])
                if sum(x is not None for x in inner_r) >= 2 else None
            ),
            "mean_inner_RMSE": mean_finite(inner_rmse),
            "mean_inner_MAE": mean_finite(inner_mae),
            "mean_inner_bias": mean_finite(inner_bias),
        })

# 4. Rank within each outer fold and trait.
# Primary: higher pooled OOF Pearson r.
# Tie-breaker: lower pooled OOF RMSE.
ranked_rows = []

for outer in [f"{i:02d}" for i in range(1, 6)]:
    for trait in TRAITS:
        group = [
            r for r in pooled_rows
            if r["outer"] == outer and r["trait"] == trait
        ]

        def rank_key(r):
            # Primary criterion: Mean fold Pearson PA.
            # Arithmetic mean of the three inner-fold Pearson correlations.
            corr = r["mean_inner_pearson_r"]

            # Tie-breaker: lower arithmetic mean of inner-fold RMSE.
            rmse = r["mean_inner_RMSE"]

            return (
                0 if corr is not None else 1,
                -(corr if corr is not None else 0),
                rmse if rmse is not None else math.inf,
                r["marker_set"],
            )

        group.sort(key=rank_key)

        for rank, r in enumerate(group, start=1):
            item = dict(r)
            item["rank"] = rank
            item["selected"] = "YES" if rank == 1 else "NO"
            ranked_rows.append(item)

# 5. Write outputs.
metric_fields = list(all_metrics[0].keys())
all_metrics.sort(
    key=lambda r: (r["outer"], r["inner"], r["trait"], r["marker_set"])
)
write_tsv(OUT / "all_inner_fold_metrics.tsv", all_metrics, metric_fields)

pooled_fields = list(pooled_rows[0].keys())
pooled_rows.sort(key=lambda r: (r["outer"], r["trait"], r["marker_set"]))
write_tsv(
    OUT / "inner_cv_pooled_oof_metrics.tsv", pooled_rows, pooled_fields
)

ranked_fields = ["rank", "selected"] + pooled_fields
write_tsv(
    OUT / "inner_cv_ranked_marker_sets.tsv", ranked_rows, ranked_fields
)

selected = [r for r in ranked_rows if r["selected"] == "YES"]
if len(selected) != 10:
    stop(f"Expected 10 selected strategies, got {len(selected)}")

write_tsv(OUT / "selected_marker_sets.tsv", selected, ranked_fields)

print("\n===== SUMMARY PASS =====")
print("Folds:", 15)
print("Model metric rows:", len(all_metrics))
print("Prediction rows:", len(all_predictions))
print("Pooled OOF groups:", len(pooled_rows))
print("Selected strategies:", len(selected))
print("Output directory:", OUT)
print("OVERALL PASS")
