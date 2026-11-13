#!/usr/bin/env python3

import gzip
import argparse
import math
import sys


# ============================================================
# Stable online Pearson correlation
# Used for pooled/global R²
# ============================================================

class OnlineCorrelation:
    def __init__(self):
        self.n = 0
        self.mean_x = 0.0
        self.mean_y = 0.0
        self.M2x = 0.0
        self.M2y = 0.0
        self.C = 0.0

    def add(self, x, y):
        self.n += 1

        dx = x - self.mean_x
        self.mean_x += dx / self.n

        dy = y - self.mean_y
        self.mean_y += dy / self.n

        self.M2x += dx * (x - self.mean_x)
        self.M2y += dy * (y - self.mean_y)

        self.C += dx * (y - self.mean_y)

    def r2(self):
        if self.n < 3:
            return None

        if self.M2x <= 0 or self.M2y <= 0:
            return None

        r = self.C / math.sqrt(self.M2x * self.M2y)
        return r * r


# ============================================================
# VCF utilities
# ============================================================

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r")


def get_samples(path):
    with open_vcf(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]

    raise RuntimeError(f"#CHROM header not found in {path}")


def first_data_line(handle):
    for line in handle:
        if not line.startswith("#"):
            return line
    return None


def normalize_chrom(chrom):
    chrom = chrom.replace("chr", "").replace("Chr", "")

    try:
        return int(chrom)
    except ValueError:
        return chrom


def parse_key(cols):
    return (normalize_chrom(cols[0]), int(cols[1]))


def get_format_index(format_string, field):
    fields = format_string.split(":")
    try:
        return fields.index(field)
    except ValueError:
        return None


def extract_sample_field(sample_string, field_index):
    if field_index is None:
        return None

    parts = sample_string.split(":")

    if field_index >= len(parts):
        return None

    value = parts[field_index]

    if value in (".", "./.", ".|.", ""):
        return None

    return value


# ============================================================
# Genotype / dosage / GP parsing
# ============================================================

def gt_to_dosage(gt):
    """
    Convert biallelic diploid GT to ALT dosage:
    0/0 or 0|0 -> 0
    0/1, 1/0, 0|1, 1|0 -> 1
    1/1 or 1|1 -> 2
    """

    if gt is None:
        return None

    gt = gt.replace("|", "/")

    if gt == "0/0":
        return 0.0
    elif gt in ("0/1", "1/0"):
        return 1.0
    elif gt == "1/1":
        return 2.0
    else:
        return None


def parse_float(value):
    if value is None:
        return None

    try:
        x = float(value)
    except ValueError:
        return None

    if not math.isfinite(x):
        return None

    return x


def parse_gp(value):
    if value is None:
        return None

    try:
        gp = [float(x) for x in value.split(",")]
    except ValueError:
        return None

    if len(gp) != 3:
        return None

    if not all(math.isfinite(x) for x in gp):
        return None

    return gp


# ============================================================
# Per-variant dosage R²
# ============================================================

def pearson_r2(xs, ys):
    n = len(xs)

    if n < 3:
        return None

    mean_x = sum(xs) / n
    mean_y = sum(ys) / n

    ssx = sum((x - mean_x) ** 2 for x in xs)
    ssy = sum((y - mean_y) ** 2 for y in ys)

    if ssx <= 0 or ssy <= 0:
        return None

    ssxy = sum(
        (x - mean_x) * (y - mean_y)
        for x, y in zip(xs, ys)
    )

    r = ssxy / math.sqrt(ssx * ssy)

    return r * r


# ============================================================
# Main evaluation
# ============================================================

def evaluate(true_vcf, imputed_vcf, output_file, progress_every=100000):

    print("=" * 72)
    print("Imputation accuracy evaluation")
    print("=" * 72)
    print(f"Truth:   {true_vcf}")
    print(f"Imputed: {imputed_vcf}")
    print("=" * 72)

    # --------------------------------------------------------
    # Samples
    # --------------------------------------------------------

    true_samples = get_samples(true_vcf)
    imp_samples = get_samples(imputed_vcf)

    print(f"Truth samples:   {len(true_samples)}")
    print(f"Imputed samples: {len(imp_samples)}")

    true_set = set(true_samples)
    imp_set = set(imp_samples)

    if true_set != imp_set:
        only_true = sorted(true_set - imp_set)
        only_imp = sorted(imp_set - true_set)

        print("ERROR: sample sets differ.")

        if only_true:
            print("Only in truth:", only_true[:10])

        if only_imp:
            print("Only in imputed:", only_imp[:10])

        sys.exit(1)

    # Match by sample name, not by column position
    true_index = {s: i for i, s in enumerate(true_samples)}
    imp_index = {s: i for i, s in enumerate(imp_samples)}

    sample_pairs = [
        (true_index[s], imp_index[s])
        for s in true_samples
    ]

    if true_samples == imp_samples:
        print("Sample order: identical")
    else:
        print("Sample order: different; matched safely by sample name")

    print(f"Samples analyzed: {len(sample_pairs)}")

    # --------------------------------------------------------
    # Statistics
    # --------------------------------------------------------

    pooled_corr = OnlineCorrelation()

    matched_positions = 0
    allele_mismatch = 0
    truth_only = 0
    imputed_only = 0

    valid_pairs = 0

    per_variant_r2_sum = 0.0
    per_variant_r2_n = 0
    per_variant_r2_uncalculable = 0

    concordance_sum = 0.0
    concordance_n = 0

    # --------------------------------------------------------
    # Open files
    # --------------------------------------------------------

    true_f = open_vcf(true_vcf)
    imp_f = open_vcf(imputed_vcf)

    true_line = first_data_line(true_f)
    imp_line = first_data_line(imp_f)

    print("\nStreaming VCF files...")

    while true_line is not None and imp_line is not None:

        true_cols = true_line.rstrip("\n").split("\t")
        imp_cols = imp_line.rstrip("\n").split("\t")

        true_key = parse_key(true_cols)
        imp_key = parse_key(imp_cols)

        # ----------------------------------------------------
        # Numeric CHROM/POS matching
        # ----------------------------------------------------

        if true_key < imp_key:
            truth_only += 1
            true_line = true_f.readline() or None
            continue

        if true_key > imp_key:
            imputed_only += 1
            imp_line = imp_f.readline() or None
            continue

        # Same CHROM + POS
        matched_positions += 1

        # ----------------------------------------------------
        # REF/ALT compatibility
        # ----------------------------------------------------

        true_ref = true_cols[3]
        true_alt = true_cols[4]

        imp_ref = imp_cols[3]
        imp_alt = imp_cols[4]

        if true_ref != imp_ref or true_alt != imp_alt:

            allele_mismatch += 1

            true_line = true_f.readline() or None
            imp_line = imp_f.readline() or None
            continue

        # ----------------------------------------------------
        # FORMAT fields
        # ----------------------------------------------------

        true_format = true_cols[8]
        imp_format = imp_cols[8]

        true_gt_idx = get_format_index(true_format, "GT")

        imp_ds_idx = get_format_index(imp_format, "DS")
        imp_gp_idx = get_format_index(imp_format, "GP")

        if true_gt_idx is None:
            raise RuntimeError(
                f"GT not found in truth FORMAT at "
                f"{true_cols[0]}:{true_cols[1]}"
            )

        if imp_ds_idx is None:
            raise RuntimeError(
                f"DS not found in imputed FORMAT at "
                f"{imp_cols[0]}:{imp_cols[1]}"
            )

        if imp_gp_idx is None:
            raise RuntimeError(
                f"GP not found in imputed FORMAT at "
                f"{imp_cols[0]}:{imp_cols[1]}"
            )

        # ----------------------------------------------------
        # Per-site data
        # ----------------------------------------------------

        site_ds = []
        site_gt = []

        concordance_numerator = 0.0
        concordance_denominator = 0.0

        for true_i, imp_i in sample_pairs:

            true_sample = true_cols[9 + true_i]
            imp_sample = imp_cols[9 + imp_i]

            gt_string = extract_sample_field(
                true_sample,
                true_gt_idx
            )

            ds_string = extract_sample_field(
                imp_sample,
                imp_ds_idx
            )

            gp_string = extract_sample_field(
                imp_sample,
                imp_gp_idx
            )

            true_gt = gt_to_dosage(gt_string)
            ds = parse_float(ds_string)
            gp = parse_gp(gp_string)

            # -----------------------------------------------
            # R²
            # -----------------------------------------------

            if true_gt is not None and ds is not None:

                pooled_corr.add(ds, true_gt)

                valid_pairs += 1

                site_ds.append(ds)
                site_gt.append(true_gt)

            # -----------------------------------------------
            # GP-weighted concordance
            # -----------------------------------------------

            if true_gt is not None and gp is not None:

                gt_index = int(true_gt)

                gp_total = sum(gp)

                if gp_total > 0:

                    concordance_numerator += gp[gt_index]
                    concordance_denominator += gp_total

        # ----------------------------------------------------
        # Per-variant dosage R²
        # ----------------------------------------------------

        site_r2 = pearson_r2(site_ds, site_gt)

        if site_r2 is not None:
            per_variant_r2_sum += site_r2
            per_variant_r2_n += 1
        else:
            per_variant_r2_uncalculable += 1

        # ----------------------------------------------------
        # Per-variant GP-weighted concordance
        # ----------------------------------------------------

        if concordance_denominator > 0:

            site_concordance = (
                concordance_numerator /
                concordance_denominator
            )

            concordance_sum += site_concordance
            concordance_n += 1

        # ----------------------------------------------------
        # Advance both
        # ----------------------------------------------------

        true_line = true_f.readline() or None
        imp_line = imp_f.readline() or None

        if (
            progress_every > 0
            and matched_positions % progress_every == 0
        ):
            print(
                f"Matched positions: {matched_positions:,} | "
                f"valid pairs: {valid_pairs:,}"
            )

    # Remaining records
    while true_line is not None:
        truth_only += 1
        true_line = true_f.readline() or None

    while imp_line is not None:
        imputed_only += 1
        imp_line = imp_f.readline() or None

    true_f.close()
    imp_f.close()

    # --------------------------------------------------------
    # Final metrics
    # --------------------------------------------------------

    global_r2 = pooled_corr.r2()

    if per_variant_r2_n > 0:
        mean_per_variant_r2 = (
            per_variant_r2_sum /
            per_variant_r2_n
        )
    else:
        mean_per_variant_r2 = None

    if concordance_n > 0:
        mean_concordance = (
            concordance_sum /
            concordance_n
        )
    else:
        mean_concordance = None

    # --------------------------------------------------------
    # Print
    # --------------------------------------------------------

    print("\n" + "=" * 72)
    print("RESULTS")
    print("=" * 72)

    print(f"Matched positions:              {matched_positions:,}")
    print(f"REF/ALT mismatches:             {allele_mismatch:,}")
    print(f"Truth-only positions:           {truth_only:,}")
    print(f"Imputed-only positions:         {imputed_only:,}")

    print(f"Valid DS-GT pairs:              {valid_pairs:,}")

    if global_r2 is not None:
        print(f"Pooled empirical R2:            {global_r2:.10f}")
    else:
        print("Pooled empirical R2:            NA")

    print(f"Variants with calculable R2:    {per_variant_r2_n:,}")
    print(
        f"Variants without calculable R2: "
        f"{per_variant_r2_uncalculable:,}"
    )

    if mean_per_variant_r2 is not None:
        print(
            f"Mean per-variant dosage R2:      "
            f"{mean_per_variant_r2:.10f}"
        )
    else:
        print("Mean per-variant dosage R2:      NA")

    print(
        f"Variants with concordance:      "
        f"{concordance_n:,}"
    )

    if mean_concordance is not None:
        print(
            f"Mean GP-weighted concordance:    "
            f"{mean_concordance:.10f}"
        )
    else:
        print("Mean GP-weighted concordance:    NA")

    print("=" * 72)

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    with open(output_file, "w") as out:

        out.write("Metric\tValue\n")

        out.write(
            f"Matched_positions\t"
            f"{matched_positions}\n"
        )

        out.write(
            f"Allele_mismatch\t"
            f"{allele_mismatch}\n"
        )

        out.write(
            f"Truth_only_positions\t"
            f"{truth_only}\n"
        )

        out.write(
            f"Imputed_only_positions\t"
            f"{imputed_only}\n"
        )

        out.write(
            f"Valid_DS_GT_pairs\t"
            f"{valid_pairs}\n"
        )

        out.write(
            f"Pooled_empirical_R2\t"
            f"{global_r2 if global_r2 is not None else 'NA'}\n"
        )

        out.write(
            f"Variants_with_calculable_R2\t"
            f"{per_variant_r2_n}\n"
        )

        out.write(
            f"Variants_without_calculable_R2\t"
            f"{per_variant_r2_uncalculable}\n"
        )

        out.write(
            f"Mean_per_variant_dosage_R2\t"
            f"{mean_per_variant_r2 if mean_per_variant_r2 is not None else 'NA'}\n"
        )

        out.write(
            f"Variants_with_concordance\t"
            f"{concordance_n}\n"
        )

        out.write(
            f"Mean_GP_weighted_concordance\t"
            f"{mean_concordance if mean_concordance is not None else 'NA'}\n"
        )

    print(f"\nResults written to: {output_file}")


# ============================================================
# Command line
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Evaluate imputation accuracy using pooled empirical R2, "
            "mean per-variant dosage R2, and GP-weighted concordance."
        )
    )

    parser.add_argument(
        "-t",
        "--truth",
        required=True,
        help="Truth VCF/VCF.gz"
    )

    parser.add_argument(
        "-i",
        "--imputed",
        required=True,
        help="IMPUTE5 output VCF/VCF.gz"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output summary TSV"
    )

    parser.add_argument(
        "--progress",
        type=int,
        default=100000,
        help="Print progress every N matched variants"
    )

    args = parser.parse_args()

    evaluate(
        args.truth,
        args.imputed,
        args.output,
        args.progress
    )


if __name__ == "__main__":
    main()
