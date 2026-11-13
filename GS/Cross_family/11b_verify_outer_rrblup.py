#!/usr/bin/env python3
"""Verify all 10 outer-fold/trait RRBLUP outputs against Stage 10 model inputs."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError(msg)

def tsv(path: Path):
    need(path.is_file() and path.stat().st_size > 0, f"missing/empty TSV: {path}")
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def flag(path: Path):
    need(path.is_file() and path.stat().st_size > 0, f"missing/empty flag: {path}")
    d = {}
    for n, line in enumerate(path.read_text().splitlines(), 1):
        k, sep, v = line.partition("=")
        need(bool(sep) and k not in d, f"malformed/duplicate flag at {path}:{n}")
        d[k] = v
    return d

def ids(path: Path):
    need(path.is_file() and path.stat().st_size > 0, f"missing IDs: {path}")
    out = path.read_text().splitlines()
    need(out and all(out) and len(out) == len(set(out)), f"invalid/duplicate IDs: {path}")
    return out

def num(x, label):
    try:
        v = float(x)
    except Exception:
        raise RuntimeError(f"non-numeric {label}: {x!r}")
    need(math.isfinite(v), f"non-finite {label}: {x!r}")
    return v

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True, help="GS-reviewer project root")
    ap.add_argument("--write-summary", action="store_true",
                    help="write a combined 10-row summary TSV after successful verification")
    args = ap.parse_args()
    work = Path(args.work).resolve()
    base_root = work / "Cross_family/results/10_outer_model_data"
    out_root = work / "Cross_family/results/11_outer_rrblup"
    good, bad, summary = 0, 0, []
    for oi in range(1, 6):
        outer = f"{oi:02d}"
        for trait in ("BB", "H6"):
            try:
                inp = base_root / f"outer_{outer}" / trait
                out = out_root / f"outer_{outer}" / trait
                im = tsv(inp / "final_genotype_manifest.tsv")
                need(len(im) == 1, f"expected one input manifest row: {inp}")
                im = im[0]
                n_chip, n_qtn = int(im["N_chip"]), int(im["N_QTN"])
                n_train, n_test = int(im["N_train"]), int(im["N_test"])
                need(len(ids(inp / "train.ids")) == n_train, f"train ID count mismatch {outer}/{trait}")
                test_ids = ids(inp / "test.ids")
                need(len(test_ids) == n_test, f"test ID count mismatch {outer}/{trait}")
                fl = flag(out / "complete.flag")
                for k, v in {
                    "status":"PASS", "outer":outer, "trait":trait,
                    "N_chip":str(n_chip), "N_qtn":str(n_qtn),
                    "N_train":str(n_train), "N_test":str(n_test),
                    "test_phenotype_used_for_fitting":"NO"
                }.items():
                    need(fl.get(k) == v, f"output flag mismatch {k}={v} at {out}")
                pred = tsv(out / "outer_test_predictions.tsv")
                need(len(pred) == n_test, f"prediction row count mismatch {outer}/{trait}")
                need([r["IID"] for r in pred] == test_ids, f"prediction ID/order mismatch {outer}/{trait}")
                for r in pred:
                    for col in ("Observed_raw","RRBLUP_CHIP_raw","RRBLUP_COMBINED_raw"):
                        num(r[col], f"{outer}/{trait}/{col}")
                test_pheno = tsv(inp / "evaluation_only/test.pheno.tsv")
                need(len(test_pheno) == n_test, f"test phenotype count mismatch {outer}/{trait}")
                need([r["ID"] for r in test_pheno] == test_ids, f"test phenotype ID/order mismatch {outer}/{trait}")
                for a, b in zip(pred, test_pheno):
                    need(abs(num(a["Observed_raw"],"observed") - num(b["phenotype"],"test phenotype")) <=
                         1e-10 * max(1.0, abs(float(b["phenotype"]))),
                         f"observed phenotype mismatch for {outer}/{trait}/{a['IID']}")
                pred_only = tsv(out / "outer_test_predictions.predictions_only.tsv")
                need(len(pred_only) == n_test and [r["IID"] for r in pred_only] == test_ids,
                     f"prediction-only output row/ID mismatch {outer}/{trait}")
                for r in pred_only:
                    num(r["RRBLUP_CHIP_raw"], "CHIP prediction")
                    num(r["RRBLUP_COMBINED_raw"], "combined prediction")
                met = tsv(out / "outer_metrics.tsv")
                need(len(met) == 1, f"expected one metrics row {outer}/{trait}")
                m = met[0]
                need(m["Outer"] == outer and m["Trait"] == trait, f"metrics fold/trait mismatch {outer}/{trait}")
                need(int(m["N_chip_selected"]) == n_chip and int(m["N_QTN_selected"]) == n_qtn,
                     f"metrics marker count mismatch {outer}/{trait}")
                for col in ("RMSE_raw_CHIP","RMSE_raw_COMBINED","MAE_raw_CHIP","MAE_raw_COMBINED",
                            "Bias_raw_CHIP","Bias_raw_COMBINED"):
                    num(m[col], f"{outer}/{trait}/{col}")
                summary.append(m)
                good += 1
                print(f"PASS outer_{outer}/{trait}: chip={n_chip}, QTN={n_qtn}, "
                      f"train={n_train}, test={n_test}, predictions={n_test}")
            except Exception as e:
                bad += 1
                print(f"FAIL outer_{outer}/{trait}: {e}")
    print(f"OVERALL 11a VERIFY PASS={good}/10 FAIL={bad}/10")
    if bad:
        raise SystemExit(1)
    if args.write_summary:
        out_root.mkdir(parents=True, exist_ok=True)
        path = out_root / "outer_rrblup_summary.tsv"
        cols = list(summary[0].keys())
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(summary)
        print(f"SUMMARY: {path}")

if __name__ == "__main__":
    main()
