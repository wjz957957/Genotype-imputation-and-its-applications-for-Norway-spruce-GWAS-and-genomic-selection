#!/usr/bin/env python3

import gzip
import os
import numpy as np


# ============================================================
# Paths and analysis settings
# ============================================================

BASE = (
    "/cfs/klemming/projects/snic/snic2020-6-176/"
    "chen2017112/jingz-transfer/IMPUTATION/"
    "reference-population-10random/IMPUTE-population"
)

SITE_DIR = (
    "/cfs/klemming/projects/snic/snic2020-6-176/"
    "chen2017112/jingz-transfer/IMPUTATION/"
    "reference-population-10random/site_folds-no-phasing"
)

PANEL_SIZES = [50, 100, 200, 300, 400, 700, 1047]
FOLDS = range(1, 11)

# IMPORTANT:
# Keep this consistent with the previous Reviewer #7 analysis.
# The original script used INFO > 0.8, not >= 0.8.
INFO_THRESHOLD = 0.8


# ============================================================
# Basic genotype / correlation functions
# ============================================================

def parse_gt(gt):
    """
    Convert diploid GT to alternate-allele dosage:
      0/0 or 0|0 -> 0
      0/1, 1/0, 0|1, 1|0 -> 1
      1/1 or 1|1 -> 2

    Missing or invalid genotypes return None.
    """

    if gt is None or gt in (".", "./.", ".|."):
        return None

    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)

    if len(parts) != 2 or "." in parts:
        return None

    try:
        return int(parts[0]) + int(parts[1])
    except Exception:
        return None


def pearson_r2(x, y):
    """
    Squared Pearson correlation.

    In this analysis:
      x = imputed dosage (DS)
      y = observed genotype dosage (0/1/2)

    When x and y contain all valid SNP-individual observations
    within a class, this is the pooled dosage R².
    """

    if len(x) < 2:
        return np.nan

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan

    r = np.corrcoef(x, y)[0, 1]

    return r * r


# ============================================================
# NRC functions
# ============================================================

def valid_gp(gp):
    """
    Check whether GP contains three finite genotype probabilities:
      P(0/0), P(0/1), P(1/1)
    """

    if gp is None or len(gp) != 3:
        return False

    try:
        return all(np.isfinite(x) for x in gp)
    except Exception:
        return False


def nrc_from_matrix(matrix):
    """
    Calculate non-reference concordance (NRC).

    Matrix orientation:
        rows    = observed / true genotype
        columns = imputed hard genotype

                    Imputed
                  0     1     2
        True 0   N00   N01   N02
             1   N10   N11   N12
             2   N20   N21   N22

    Correct homozygous-reference calls (0 -> 0) are excluded
    from both numerator and denominator.

    NRC =
        (N11 + N22)
        -----------------------------------------------
        all genotype comparisons except correct 0 -> 0

    equivalently:

        NRC = (N11 + N22) / (Ntotal - N00)

    Returns:
        nrc
        numerator   = number of correct non-reference calls
        denominator = number of NRC-evaluable comparisons
    """

    matrix = np.asarray(matrix)

    total = int(np.sum(matrix))
    concordant_ref = int(matrix[0, 0])

    correct_nonref = int(
        matrix[1, 1] +
        matrix[2, 2]
    )

    evaluable = total - concordant_ref

    if evaluable == 0:
        return np.nan, correct_nonref, evaluable

    nrc = correct_nonref / evaluable

    return nrc, correct_nonref, evaluable


# ============================================================
# Read truth VCF
# ============================================================

def read_truth_vcf(path):
    """
    Read truth VCF.

    Returns
    -------
    samples : list
        Sample IDs.

    data : dict
        snpid -> {
            "maf": MAF calculated from AC/AN,
            "gt":  list of observed genotypes coded 0/1/2
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

            # ------------------------------------------------
            # Parse INFO
            # ------------------------------------------------

            info = {}

            for x in cols[7].split(";"):
                if "=" in x:
                    k, v = x.split("=", 1)
                    info[k] = v

            # ------------------------------------------------
            # Calculate MAF from AC / AN
            # ------------------------------------------------

            try:
                ac = float(info["AC"].split(",")[0])
                an = float(info["AN"])

                af = ac / an
                maf = min(af, 1.0 - af)

            except Exception:
                maf = np.nan

            # ------------------------------------------------
            # GT field
            # ------------------------------------------------

            fmt = cols[8].split(":")

            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                raise RuntimeError(
                    f"GT field not found in truth VCF: {path}"
                )

            gt_values = []

            for s in cols[9:]:

                parts = s.split(":")

                if gt_idx < len(parts):
                    gt_values.append(
                        parse_gt(parts[gt_idx])
                    )
                else:
                    gt_values.append(None)

            data[snpid] = {
                "maf": maf,
                "gt": gt_values
            }

    return samples, data


# ============================================================
# Read imputed VCF
# ============================================================

def read_imputed_vcf(path):
    """
    Read IMPUTE5 VCF.

    Returns
    -------
    samples : list

    data : dict
        snpid -> {
            "info": INFO score,
            "ds":   imputed dosages,
            "gp":   genotype probabilities
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

            # ------------------------------------------------
            # Parse INFO score
            # ------------------------------------------------

            info_dict = {}

            for x in cols[7].split(";"):
                if "=" in x:
                    k, v = x.split("=", 1)
                    info_dict[k] = v

            try:
                info_score = float(
                    info_dict["INFO"].split(",")[0]
                )
            except Exception:
                info_score = np.nan

            # ------------------------------------------------
            # FORMAT
            # ------------------------------------------------

            fmt = cols[8].split(":")

            try:
                ds_idx = fmt.index("DS")
            except ValueError:
                ds_idx = None

            try:
                gp_idx = fmt.index("GP")
            except ValueError:
                gp_idx = None

            ds_values = []
            gp_values = []

            # ------------------------------------------------
            # Sample fields
            # ------------------------------------------------

            for s in cols[9:]:

                parts = s.split(":")

                # DS
                try:
                    ds = (
                        float(parts[ds_idx])
                        if ds_idx is not None
                        else np.nan
                    )
                except Exception:
                    ds = np.nan

                ds_values.append(ds)

                # GP
                try:

                    if gp_idx is None:
                        raise ValueError("GP absent")

                    gp = [
                        float(x)
                        for x in parts[gp_idx].split(",")
                    ]

                    if len(gp) != 3:
                        gp = [
                            np.nan,
                            np.nan,
                            np.nan
                        ]

                except Exception:

                    gp = [
                        np.nan,
                        np.nan,
                        np.nan
                    ]

                gp_values.append(gp)

            data[snpid] = {
                "info": info_score,
                "ds": ds_values,
                "gp": gp_values
            }

    return samples, data


# ============================================================
# Bucket for Rare / Common / All
# ============================================================

def new_bucket():

    return {

        # SNP counts
        "n_validation_snps": 0,
        "n_retained_snps": 0,

        # -----------------------------------------------
        # Existing per-SNP statistics
        # Keep these unchanged for compatibility
        # -----------------------------------------------

        "per_snp_r2_raw": [],
        "per_snp_conc_raw": [],

        "per_snp_r2_filtered": [],
        "per_snp_conc_filtered": [],

        # -----------------------------------------------
        # Existing pooled dosage R²
        # -----------------------------------------------

        "pooled_ds_raw": [],
        "pooled_gt_raw": [],

        "pooled_ds_filtered": [],
        "pooled_gt_filtered": [],

        # -----------------------------------------------
        # Existing pooled GP-based concordance
        # -----------------------------------------------

        "gp_correct_raw_sum": 0.0,
        "gp_correct_raw_n": 0,

        "gp_correct_filtered_sum": 0.0,
        "gp_correct_filtered_n": 0,

        # -----------------------------------------------
        # NEW: pooled hard-call confusion matrices for NRC
        #
        # rows    = TRUE genotype
        # columns = IMPUTED MAP genotype
        # -----------------------------------------------

        "nrc_matrix_raw": np.zeros(
            (3, 3),
            dtype=np.int64
        ),

        "nrc_matrix_filtered": np.zeros(
            (3, 3),
            dtype=np.int64
        )
    }


# ============================================================
# Analyze one SNP
# ============================================================

def analyze_snp(gt_values, ds_values, gp_values):
    """
    Analyze one SNP.

    Returns
    -------
    r2 :
        Per-SNP dosage R².

    concordance :
        Per-SNP GP-based concordance.

    valid_ds / valid_gt :
        Used later for pooled dosage R².

    gp_correct :
        Probability assigned to the observed genotype
        for each valid sample; used later for pooled
        GP-based concordance.

    nrc_matrix :
        3 x 3 hard-call confusion matrix for this SNP,
        used later for pooled NRC.
    """

    valid_ds = []
    valid_gt = []

    gp_correct = []

    # rows = true GT
    # cols = imputed MAP GT
    nrc_matrix = np.zeros(
        (3, 3),
        dtype=np.int64
    )

    for gt, ds, gp in zip(
        gt_values,
        ds_values,
        gp_values
    ):

        if gt is None:
            continue

        # ------------------------------------------------
        # Dosage R²
        # ------------------------------------------------

        if not np.isnan(ds):

            valid_gt.append(
                float(gt)
            )

            valid_ds.append(
                float(ds)
            )

        # ------------------------------------------------
        # GP-based concordance and NRC
        # ------------------------------------------------

        if (
            gt in (0, 1, 2)
            and valid_gp(gp)
        ):

            # Existing probabilistic concordance:
            # probability assigned to the observed genotype
            gp_correct.append(
                float(gp[gt])
            )

            # NEW:
            # Convert GP to hard genotype call for NRC
            pred_gt = int(
                np.argmax(gp)
            )

            nrc_matrix[
                gt,
                pred_gt
            ] += 1

    # ----------------------------------------------------
    # Per-SNP R²
    # ----------------------------------------------------

    r2 = pearson_r2(
        valid_ds,
        valid_gt
    )

    # ----------------------------------------------------
    # Per-SNP GP-based concordance
    # ----------------------------------------------------

    if len(gp_correct) > 0:

        concordance = float(
            np.mean(gp_correct)
        )

    else:

        concordance = np.nan

    return (
        r2,
        concordance,
        valid_ds,
        valid_gt,
        gp_correct,
        nrc_matrix
    )


# ============================================================
# Main analysis
# ============================================================

results = []


for panel in PANEL_SIZES:

    print(
        f"\n===== Reference panel size: {panel} ====="
    )

    buckets = {
        "Rare": new_bucket(),
        "Common": new_bucket(),
        "All": new_bucket()
    }

    # ========================================================
    # Ten site-based CV folds
    # ========================================================

    for fold in FOLDS:

        truth_vcf = (
            f"{SITE_DIR}/"
            f"evaluation_snps_fold{fold}_AC.vcf.gz"
        )

        imputed_vcf = (
            f"{BASE}/population-{panel}/"
            f"repeat-1/fold{fold}/"
            f"impute5_output_{panel}_1_fold{fold}_real.vcf.gz"
        )

        # ----------------------------------------------------
        # File checks
        # ----------------------------------------------------

        if not os.path.exists(truth_vcf):
            raise FileNotFoundError(
                truth_vcf
            )

        if not os.path.exists(imputed_vcf):
            raise FileNotFoundError(
                imputed_vcf
            )

        # ----------------------------------------------------
        # Read VCFs
        # ----------------------------------------------------

        truth_samples, truth = read_truth_vcf(
            truth_vcf
        )

        imp_samples, imputed = read_imputed_vcf(
            imputed_vcf
        )

        # Exact sample ID and order check
        if truth_samples != imp_samples:

            raise RuntimeError(
                f"Sample mismatch: "
                f"panel={panel}, fold={fold}"
            )

        # ----------------------------------------------------
        # Matched validation sites
        # ----------------------------------------------------

        common_snps = (
            set(truth)
            &
            set(imputed)
        )

        print(
            f"panel={panel}, "
            f"fold={fold}, "
            f"truth={len(truth)}, "
            f"imputed={len(imputed)}, "
            f"matched={len(common_snps)}"
        )

        # ====================================================
        # SNP loop
        # ====================================================

        for snpid in common_snps:

            maf = truth[snpid]["maf"]

            if np.isnan(maf):
                continue

            # ------------------------------------------------
            # Rare / Common classification
            # ------------------------------------------------

            if maf < 0.05:
                maf_class = "Rare"
            else:
                maf_class = "Common"

            classes = [
                maf_class,
                "All"
            ]

            # ------------------------------------------------
            # Data
            # ------------------------------------------------

            gt = truth[snpid]["gt"]

            ds = imputed[snpid]["ds"]

            gp = imputed[snpid]["gp"]

            info = imputed[snpid]["info"]

            (
                r2,
                concordance,
                valid_ds,
                valid_gt,
                gp_correct,
                nrc_matrix
            ) = analyze_snp(
                gt,
                ds,
                gp
            )

            # =================================================
            # Add SNP to relevant class(es)
            # =================================================

            for cls in classes:

                b = buckets[cls]

                b[
                    "n_validation_snps"
                ] += 1

                # =============================================
                # RAW results
                # =============================================

                if not np.isnan(r2):

                    b[
                        "per_snp_r2_raw"
                    ].append(r2)

                if not np.isnan(concordance):

                    b[
                        "per_snp_conc_raw"
                    ].append(concordance)

                # pooled dosage R²
                b[
                    "pooled_ds_raw"
                ].extend(valid_ds)

                b[
                    "pooled_gt_raw"
                ].extend(valid_gt)

                # pooled GP-based concordance
                b[
                    "gp_correct_raw_sum"
                ] += sum(gp_correct)

                b[
                    "gp_correct_raw_n"
                ] += len(gp_correct)

                # NEW:
                # pooled NRC confusion matrix
                b[
                    "nrc_matrix_raw"
                ] += nrc_matrix

                # =============================================
                # INFO > 0.8 results
                # =============================================

                if (
                    not np.isnan(info)
                    and
                    info > INFO_THRESHOLD
                ):

                    b[
                        "n_retained_snps"
                    ] += 1

                    if not np.isnan(r2):

                        b[
                            "per_snp_r2_filtered"
                        ].append(r2)

                    if not np.isnan(concordance):

                        b[
                            "per_snp_conc_filtered"
                        ].append(concordance)

                    # pooled dosage R²
                    b[
                        "pooled_ds_filtered"
                    ].extend(valid_ds)

                    b[
                        "pooled_gt_filtered"
                    ].extend(valid_gt)

                    # pooled GP-based concordance
                    b[
                        "gp_correct_filtered_sum"
                    ] += sum(gp_correct)

                    b[
                        "gp_correct_filtered_n"
                    ] += len(gp_correct)

                    # NEW:
                    # pooled NRC confusion matrix
                    b[
                        "nrc_matrix_filtered"
                    ] += nrc_matrix


    # ========================================================
    # Summarize panel
    # ========================================================

    for cls in [
        "Rare",
        "Common",
        "All"
    ]:

        b = buckets[cls]

        n_total = b[
            "n_validation_snps"
        ]

        n_ret = b[
            "n_retained_snps"
        ]

        # ----------------------------------------------------
        # Retention
        # ----------------------------------------------------

        retention = (
            n_ret / n_total
            if n_total > 0
            else np.nan
        )

        # ====================================================
        # RAW
        # ====================================================

        mean_r2_raw = (
            np.mean(
                b["per_snp_r2_raw"]
            )
            if b["per_snp_r2_raw"]
            else np.nan
        )

        mean_conc_raw = (
            np.mean(
                b["per_snp_conc_raw"]
            )
            if b["per_snp_conc_raw"]
            else np.nan
        )

        pooled_r2_raw = pearson_r2(
            b["pooled_ds_raw"],
            b["pooled_gt_raw"]
        )

        pooled_conc_raw = (
            b["gp_correct_raw_sum"]
            /
            b["gp_correct_raw_n"]
            if b["gp_correct_raw_n"] > 0
            else np.nan
        )

        (
            pooled_nrc_raw,
            nrc_num_raw,
            nrc_den_raw
        ) = nrc_from_matrix(
            b["nrc_matrix_raw"]
        )

        # ====================================================
        # INFO > 0.8
        # ====================================================

        mean_r2_filt = (
            np.mean(
                b["per_snp_r2_filtered"]
            )
            if b["per_snp_r2_filtered"]
            else np.nan
        )

        mean_conc_filt = (
            np.mean(
                b["per_snp_conc_filtered"]
            )
            if b["per_snp_conc_filtered"]
            else np.nan
        )

        pooled_r2_filt = pearson_r2(
            b["pooled_ds_filtered"],
            b["pooled_gt_filtered"]
        )

        pooled_conc_filt = (
            b["gp_correct_filtered_sum"]
            /
            b["gp_correct_filtered_n"]
            if b[
                "gp_correct_filtered_n"
            ] > 0
            else np.nan
        )

        (
            pooled_nrc_filt,
            nrc_num_filt,
            nrc_den_filt
        ) = nrc_from_matrix(
            b["nrc_matrix_filtered"]
        )

        # ----------------------------------------------------
        # Store results
        # ----------------------------------------------------

        results.append([
            panel,
            cls,
            n_total,
            n_ret,
            retention,

            # RAW
            mean_r2_raw,
            pooled_r2_raw,
            mean_conc_raw,
            pooled_conc_raw,
            pooled_nrc_raw,

            # INFO > 0.8
            mean_r2_filt,
            pooled_r2_filt,
            mean_conc_filt,
            pooled_conc_filt,
            pooled_nrc_filt,

            # NRC QC counts
            nrc_num_raw,
            nrc_den_raw,
            nrc_num_filt,
            nrc_den_filt
        ])

        # ----------------------------------------------------
        # Console summary
        # ----------------------------------------------------

        print(
            f"\nPanel={panel}, "
            f"Class={cls}"
        )

        print(
            f"  Validation SNPs: "
            f"{n_total}"
        )

        print(
            f"  Retained SNPs:   "
            f"{n_ret}"
        )

        print(
            f"  Retention:       "
            f"{retention:.6f}"
            if not np.isnan(retention)
            else
            "  Retention:       NA"
        )

        print(
            f"  RAW pooled R2:   "
            f"{pooled_r2_raw:.6f}"
            if not np.isnan(pooled_r2_raw)
            else
            "  RAW pooled R2:   NA"
        )

        print(
            f"  RAW pooled Conc: "
            f"{pooled_conc_raw:.6f}"
            if not np.isnan(pooled_conc_raw)
            else
            "  RAW pooled Conc: NA"
        )

        print(
            f"  RAW pooled NRC:  "
            f"{pooled_nrc_raw:.6f}"
            if not np.isnan(pooled_nrc_raw)
            else
            "  RAW pooled NRC:  NA"
        )

        print(
            f"  INFO>0.8 R2:     "
            f"{pooled_r2_filt:.6f}"
            if not np.isnan(pooled_r2_filt)
            else
            "  INFO>0.8 R2:     NA"
        )

        print(
            f"  INFO>0.8 Conc:   "
            f"{pooled_conc_filt:.6f}"
            if not np.isnan(pooled_conc_filt)
            else
            "  INFO>0.8 Conc:   NA"
        )

        print(
            f"  INFO>0.8 NRC:    "
            f"{pooled_nrc_filt:.6f}"
            if not np.isnan(pooled_nrc_filt)
            else
            "  INFO>0.8 NRC:    NA"
        )


# ============================================================
# Output table
# ============================================================

header = [
    "Reference_panel_size",
    "MAF_class",
    "N_validation_SNPs",
    "N_INFOgt0.8",
    "Retention_INFOgt0.8",

    # RAW
    "Mean_perSNP_R2_raw",
    "Pooled_R2_raw",
    "Mean_perSNP_Concordance_raw",
    "Pooled_Concordance_raw",
    "Pooled_NRC_raw",

    # INFO > 0.8
    "Mean_perSNP_R2_INFOgt0.8",
    "Pooled_R2_INFOgt0.8",
    "Mean_perSNP_Concordance_INFOgt0.8",
    "Pooled_Concordance_INFOgt0.8",
    "Pooled_NRC_INFOgt0.8",

    # NRC QC information
    "NRC_correct_nonref_raw",
    "NRC_evaluable_raw",
    "NRC_correct_nonref_INFOgt0.8",
    "NRC_evaluable_INFOgt0.8"
]


outfile = (
    "rare_common_panel_size_summary_with_NRC.tsv"
)


with open(
    outfile,
    "w"
) as out:

    out.write(
        "\t".join(header)
        +
        "\n"
    )

    for row in results:

        formatted = []

        for x in row:

            if isinstance(
                x,
                (float, np.floating)
            ):

                formatted.append(
                    "NA"
                    if np.isnan(x)
                    else f"{x:.6f}"
                )

            else:

                formatted.append(
                    str(x)
                )

        out.write(
            "\t".join(formatted)
            +
            "\n"
        )


print(
    "\n=========================================="
)

print(
    "Analysis complete"
)

print(
    f"Output: {outfile}"
)

print(
    "=========================================="
)
