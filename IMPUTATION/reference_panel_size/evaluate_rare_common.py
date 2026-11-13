#!/usr/bin/env python3

import gzip
import os
import math
import numpy as np

BASE = "/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/reference-population-10random/IMPUTE-population"
SITE_DIR = "/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/reference-population-10random/site_folds-no-phasing"

PANEL_SIZES = [50, 100, 200, 300, 400, 700, 1047]
FOLDS = range(1, 11)
INFO_THRESHOLD = 0.8


def parse_gt(gt):
    if gt is None or gt in (".", "./.", ".|."):
        return None

    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)

    if len(parts) != 2 or "." in parts:
        return None

    try:
        return int(parts[0]) + int(parts[1])
    except:
        return None


def pearson_r2(x, y):
    if len(x) < 2:
        return np.nan

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan

    r = np.corrcoef(x, y)[0, 1]
    return r * r


def read_truth_vcf(path):
    """
    Returns:
      samples: sample IDs
      data:
        snpid -> {
            maf,
            gt[]
        }
    """

    samples = []
    data = {}

    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                cols = line.rstrip().split("\t")
                samples = cols[9:]
                continue

            cols = line.rstrip().split("\t")

            chrom = cols[0]
            pos = cols[1]
            snpid = f"{chrom}:{pos}"

            info = {}
            for x in cols[7].split(";"):
                if "=" in x:
                    k, v = x.split("=", 1)
                    info[k] = v

            try:
                ac = float(info["AC"].split(",")[0])
                an = float(info["AN"])
                af = ac / an
                maf = min(af, 1.0 - af)
            except:
                maf = np.nan

            fmt = cols[8].split(":")
            gt_idx = fmt.index("GT")

            gt_values = []

            for s in cols[9:]:
                parts = s.split(":")
                if gt_idx < len(parts):
                    gt_values.append(parse_gt(parts[gt_idx]))
                else:
                    gt_values.append(None)

            data[snpid] = {
                "maf": maf,
                "gt": gt_values
            }

    return samples, data


def read_imputed_vcf(path):
    """
    Returns:
      samples
      snpid -> {
          info,
          ds[],
          gp[]
      }
    """

    samples = []
    data = {}

    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                cols = line.rstrip().split("\t")
                samples = cols[9:]
                continue

            cols = line.rstrip().split("\t")

            chrom = cols[0]
            pos = cols[1]
            snpid = f"{chrom}:{pos}"

            info_dict = {}
            for x in cols[7].split(";"):
                if "=" in x:
                    k, v = x.split("=", 1)
                    info_dict[k] = v

            try:
                info_score = float(info_dict["INFO"].split(",")[0])
            except:
                info_score = np.nan

            fmt = cols[8].split(":")

            try:
                ds_idx = fmt.index("DS")
            except:
                ds_idx = None

            try:
                gp_idx = fmt.index("GP")
            except:
                gp_idx = None

            ds_values = []
            gp_values = []

            for s in cols[9:]:
                parts = s.split(":")

                # DS
                try:
                    ds = float(parts[ds_idx]) if ds_idx is not None else np.nan
                except:
                    ds = np.nan

                ds_values.append(ds)

                # GP
                try:
                    gp = [float(x) for x in parts[gp_idx].split(",")]
                    if len(gp) != 3:
                        gp = [np.nan, np.nan, np.nan]
                except:
                    gp = [np.nan, np.nan, np.nan]

                gp_values.append(gp)

            data[snpid] = {
                "info": info_score,
                "ds": ds_values,
                "gp": gp_values
            }

    return samples, data


def new_bucket():
    return {
        "n_validation_snps": 0,
        "n_retained_snps": 0,

        "per_snp_r2_raw": [],
        "per_snp_conc_raw": [],

        "per_snp_r2_filtered": [],
        "per_snp_conc_filtered": [],

        "pooled_ds_raw": [],
        "pooled_gt_raw": [],

        "pooled_ds_filtered": [],
        "pooled_gt_filtered": [],

        "gp_correct_raw_sum": 0.0,
        "gp_correct_raw_n": 0,

        "gp_correct_filtered_sum": 0.0,
        "gp_correct_filtered_n": 0
    }


def analyze_snp(gt_values, ds_values, gp_values):
    valid_ds = []
    valid_gt = []

    gp_correct = []

    for gt, ds, gp in zip(gt_values, ds_values, gp_values):

        if gt is None:
            continue

        if not np.isnan(ds):
            valid_gt.append(float(gt))
            valid_ds.append(float(ds))

        if (
            gp is not None
            and len(gp) == 3
            and not np.isnan(gp[0])
            and gt in (0, 1, 2)
        ):
            gp_correct.append(gp[gt])

    r2 = pearson_r2(valid_ds, valid_gt)

    if len(gp_correct) > 0:
        concordance = float(np.mean(gp_correct))
    else:
        concordance = np.nan

    return (
        r2,
        concordance,
        valid_ds,
        valid_gt,
        gp_correct
    )


results = []

for panel in PANEL_SIZES:

    print(f"\n===== Reference panel size: {panel} =====")

    buckets = {
        "Rare": new_bucket(),
        "Common": new_bucket(),
        "All": new_bucket()
    }

    for fold in FOLDS:

        truth_vcf = (
            f"{SITE_DIR}/evaluation_snps_fold{fold}_AC.vcf.gz"
        )

        imputed_vcf = (
            f"{BASE}/population-{panel}/repeat-1/fold{fold}/"
            f"impute5_output_{panel}_1_fold{fold}_real.vcf.gz"
        )

        if not os.path.exists(truth_vcf):
            raise FileNotFoundError(truth_vcf)

        if not os.path.exists(imputed_vcf):
            raise FileNotFoundError(imputed_vcf)

        truth_samples, truth = read_truth_vcf(truth_vcf)
        imp_samples, imputed = read_imputed_vcf(imputed_vcf)

        if truth_samples != imp_samples:
            raise RuntimeError(
                f"Sample mismatch: panel={panel}, fold={fold}"
            )

        common_snps = set(truth) & set(imputed)

        print(
            f"panel={panel}, fold={fold}, "
            f"truth={len(truth)}, "
            f"imputed={len(imputed)}, "
            f"matched={len(common_snps)}"
        )

        for snpid in common_snps:

            maf = truth[snpid]["maf"]

            if np.isnan(maf):
                continue

            if maf < 0.05:
                maf_class = "Rare"
            else:
                maf_class = "Common"

            classes = [maf_class, "All"]

            gt = truth[snpid]["gt"]
            ds = imputed[snpid]["ds"]
            gp = imputed[snpid]["gp"]
            info = imputed[snpid]["info"]

            (
                r2,
                concordance,
                valid_ds,
                valid_gt,
                gp_correct
            ) = analyze_snp(gt, ds, gp)

            for cls in classes:

                b = buckets[cls]

                b["n_validation_snps"] += 1

                # RAW
                if not np.isnan(r2):
                    b["per_snp_r2_raw"].append(r2)

                if not np.isnan(concordance):
                    b["per_snp_conc_raw"].append(concordance)

                b["pooled_ds_raw"].extend(valid_ds)
                b["pooled_gt_raw"].extend(valid_gt)

                b["gp_correct_raw_sum"] += sum(gp_correct)
                b["gp_correct_raw_n"] += len(gp_correct)

                # INFO > 0.8
                if not np.isnan(info) and info > INFO_THRESHOLD:

                    b["n_retained_snps"] += 1

                    if not np.isnan(r2):
                        b["per_snp_r2_filtered"].append(r2)

                    if not np.isnan(concordance):
                        b["per_snp_conc_filtered"].append(concordance)

                    b["pooled_ds_filtered"].extend(valid_ds)
                    b["pooled_gt_filtered"].extend(valid_gt)

                    b["gp_correct_filtered_sum"] += sum(gp_correct)
                    b["gp_correct_filtered_n"] += len(gp_correct)

    for cls in ["Rare", "Common", "All"]:

        b = buckets[cls]

        n_total = b["n_validation_snps"]
        n_ret = b["n_retained_snps"]

        retention = (
            n_ret / n_total
            if n_total > 0
            else np.nan
        )

        mean_r2_raw = (
            np.mean(b["per_snp_r2_raw"])
            if b["per_snp_r2_raw"]
            else np.nan
        )

        mean_conc_raw = (
            np.mean(b["per_snp_conc_raw"])
            if b["per_snp_conc_raw"]
            else np.nan
        )

        pooled_r2_raw = pearson_r2(
            b["pooled_ds_raw"],
            b["pooled_gt_raw"]
        )

        pooled_conc_raw = (
            b["gp_correct_raw_sum"] /
            b["gp_correct_raw_n"]
            if b["gp_correct_raw_n"] > 0
            else np.nan
        )

        mean_r2_filt = (
            np.mean(b["per_snp_r2_filtered"])
            if b["per_snp_r2_filtered"]
            else np.nan
        )

        mean_conc_filt = (
            np.mean(b["per_snp_conc_filtered"])
            if b["per_snp_conc_filtered"]
            else np.nan
        )

        pooled_r2_filt = pearson_r2(
            b["pooled_ds_filtered"],
            b["pooled_gt_filtered"]
        )

        pooled_conc_filt = (
            b["gp_correct_filtered_sum"] /
            b["gp_correct_filtered_n"]
            if b["gp_correct_filtered_n"] > 0
            else np.nan
        )

        results.append([
            panel,
            cls,
            n_total,
            n_ret,
            retention,
            mean_r2_raw,
            pooled_r2_raw,
            mean_conc_raw,
            pooled_conc_raw,
            mean_r2_filt,
            pooled_r2_filt,
            mean_conc_filt,
            pooled_conc_filt
        ])


header = [
    "Reference_panel_size",
    "MAF_class",
    "N_validation_SNPs",
    "N_INFOgt0.8",
    "Retention_INFOgt0.8",
    "Mean_perSNP_R2_raw",
    "Pooled_R2_raw",
    "Mean_perSNP_Concordance_raw",
    "Pooled_Concordance_raw",
    "Mean_perSNP_R2_INFOgt0.8",
    "Pooled_R2_INFOgt0.8",
    "Mean_perSNP_Concordance_INFOgt0.8",
    "Pooled_Concordance_INFOgt0.8"
]

outfile = "rare_common_panel_size_summary.tsv"

with open(outfile, "w") as out:

    out.write("\t".join(header) + "\n")

    for row in results:

        formatted = []

        for x in row:
            if isinstance(x, float):
                formatted.append(
                    "NA" if np.isnan(x) else f"{x:.6f}"
                )
            else:
                formatted.append(str(x))

        out.write("\t".join(formatted) + "\n")


print("\n==========================================")
print("Analysis complete")
print(f"Output: {outfile}")
print("==========================================")
