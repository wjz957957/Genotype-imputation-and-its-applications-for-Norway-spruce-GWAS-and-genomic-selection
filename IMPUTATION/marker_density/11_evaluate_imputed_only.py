#!/usr/bin/env python3

import gzip
import argparse
import math
import sys


# ============================================================
# Helpers
# ============================================================

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom(chrom):
    chrom = chrom.lower()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    try:
        return int(chrom)
    except ValueError:
        return chrom


def get_samples(path):
    with open_vcf(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise RuntimeError(f"#CHROM header not found: {path}")


def first_data_line(handle):
    for line in handle:
        if not line.startswith("#"):
            return line
    return None


def site_key(cols):
    return (normalize_chrom(cols[0]), int(cols[1]))


def get_format_index(fmt, field):
    fields = fmt.split(":")
    try:
        return fields.index(field)
    except ValueError:
        return None


def get_field(sample_string, idx):
    if idx is None:
        return None

    parts = sample_string.split(":")
    if idx >= len(parts):
        return None

    value = parts[idx]

    if value in ("", ".", "./.", ".|."):
        return None

    return value


def gt_to_dosage(gt):
    """Convert biallelic diploid GT to ALT dosage 0/1/2."""
    if gt is None:
        return None

    gt = gt.replace("|", "/")

    if gt == "0/0":
        return 0
    if gt in ("0/1", "1/0"):
        return 1
    if gt == "1/1":
        return 2

    return None


def parse_float(value):
    if value is None:
        return None
    try:
        x = float(value)
    except ValueError:
        return None

    return x if math.isfinite(x) else None


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


def pearson_r2(xs, ys):
    """Squared Pearson correlation for one variant."""
    n = len(xs)

    if n < 3:
        return None

    mx = sum(xs) / n
    my = sum(ys) / n

    ssx = sum((x - mx) ** 2 for x in xs)
    ssy = sum((y - my) ** 2 for y in ys)

    if ssx <= 0 or ssy <= 0:
        return None

    ssxy = sum(
        (x - mx) * (y - my)
        for x, y in zip(xs, ys)
    )

    r = ssxy / math.sqrt(ssx * ssy)

    return r * r


# ============================================================
# Stable pooled correlation
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
# Read input scaffold sites
# ============================================================

def load_excluded_sites(path):

    sites = set()

    with open(path) as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 2:
                continue

            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])

            sites.add((chrom, pos))

    return sites


# ============================================================
# Evaluation
# ============================================================

def evaluate(truth_vcf, imputed_vcf, exclude_file,
             output_file, progress_every):

    print("=" * 78)
    print("IMPUTED-ONLY ACCURACY EVALUATION")
    print("=" * 78)

    print(f"Truth VCF:       {truth_vcf}")
    print(f"Imputed VCF:     {imputed_vcf}")
    print(f"Scaffold sites:  {exclude_file}")

    # --------------------------------------------------------
    # Read scaffold SNPs
    # --------------------------------------------------------

    exclude_sites = load_excluded_sites(exclude_file)

    print(f"\nInput scaffold SNPs loaded: {len(exclude_sites):,}")

    # --------------------------------------------------------
    # Match samples by NAME
    # --------------------------------------------------------

    true_samples = get_samples(truth_vcf)
    imp_samples = get_samples(imputed_vcf)

    if set(true_samples) != set(imp_samples):

        print("ERROR: truth and imputed sample sets differ.")
        sys.exit(1)

    true_idx = {s: i for i, s in enumerate(true_samples)}
    imp_idx = {s: i for i, s in enumerate(imp_samples)}

    sample_pairs = [
        (true_idx[s], imp_idx[s])
        for s in true_samples
    ]

    print(f"Truth samples:   {len(true_samples)}")
    print(f"Imputed samples: {len(imp_samples)}")

    if true_samples == imp_samples:
        print("Sample order:    identical")
    else:
        print("Sample order:    different; matched by sample name")

    # --------------------------------------------------------
    # Counters
    # --------------------------------------------------------

    matched_positions = 0
    allele_mismatch = 0

    scaffold_excluded = 0
    evaluated_variants = 0

    truth_only = 0
    imputed_only = 0

    valid_ds_gt_pairs = 0

    # Pooled empirical R2
    pooled_corr = OnlineCorrelation()

    # Per-variant dosage R2
    variant_r2_sum = 0.0
    variant_r2_n = 0
    variant_r2_missing = 0

    # Legacy GP-weighted concordance
    gp_conc_sum = 0.0
    gp_conc_n = 0

    # Hard genotype concordance
    hard_match = 0
    hard_total = 0

    # NRC counts
    #
    # m0 = correctly called true 0/0
    # m1 = correctly called true 0/1
    # m2 = correctly called true 1/1
    #
    # e0 = incorrectly called true 0/0
    # e1 = incorrectly called true 0/1
    # e2 = incorrectly called true 1/1
    #
    m0 = m1 = m2 = 0
    e0 = e1 = e2 = 0

    # --------------------------------------------------------
    # Stream VCFs
    # --------------------------------------------------------

    true_f = open_vcf(truth_vcf)
    imp_f = open_vcf(imputed_vcf)

    true_line = first_data_line(true_f)
    imp_line = first_data_line(imp_f)

    while true_line is not None and imp_line is not None:

        true_cols = true_line.rstrip("\n").split("\t")
        imp_cols = imp_line.rstrip("\n").split("\t")

        tk = site_key(true_cols)
        ik = site_key(imp_cols)

        # Numeric chromosome/position comparison
        if tk < ik:

            truth_only += 1
            true_line = true_f.readline() or None
            continue

        if tk > ik:

            imputed_only += 1
            imp_line = imp_f.readline() or None
            continue

        matched_positions += 1

        # ----------------------------------------------------
        # Check alleles
        # ----------------------------------------------------

        if (
            true_cols[3] != imp_cols[3]
            or true_cols[4] != imp_cols[4]
        ):

            allele_mismatch += 1

            true_line = true_f.readline() or None
            imp_line = imp_f.readline() or None

            continue

        # ----------------------------------------------------
        # CRITICAL STEP:
        # remove observed scaffold/input SNPs
        # ----------------------------------------------------

        if tk in exclude_sites:

            scaffold_excluded += 1

            true_line = true_f.readline() or None
            imp_line = imp_f.readline() or None

            continue

        evaluated_variants += 1

        # ----------------------------------------------------
        # FORMAT
        # ----------------------------------------------------

        true_gt_idx = get_format_index(
            true_cols[8], "GT"
        )

        imp_gt_idx = get_format_index(
            imp_cols[8], "GT"
        )

        imp_ds_idx = get_format_index(
            imp_cols[8], "DS"
        )

        imp_gp_idx = get_format_index(
            imp_cols[8], "GP"
        )

        if true_gt_idx is None:
            raise RuntimeError(
                f"Truth GT missing at {true_cols[0]}:{true_cols[1]}"
            )

        if imp_gt_idx is None:
            raise RuntimeError(
                f"Imputed GT missing at {imp_cols[0]}:{imp_cols[1]}"
            )

        if imp_ds_idx is None:
            raise RuntimeError(
                f"Imputed DS missing at {imp_cols[0]}:{imp_cols[1]}"
            )

        if imp_gp_idx is None:
            raise RuntimeError(
                f"Imputed GP missing at {imp_cols[0]}:{imp_cols[1]}"
            )

        site_ds = []
        site_truth = []

        # for legacy GP-weighted concordance
        gp_numerator = 0.0
        gp_denominator = 0.0

        # ----------------------------------------------------
        # Samples
        # ----------------------------------------------------

        for ti, ii in sample_pairs:

            true_sample = true_cols[9 + ti]
            imp_sample = imp_cols[9 + ii]

            true_gt_string = get_field(
                true_sample, true_gt_idx
            )

            imp_gt_string = get_field(
                imp_sample, imp_gt_idx
            )

            ds_string = get_field(
                imp_sample, imp_ds_idx
            )

            gp_string = get_field(
                imp_sample, imp_gp_idx
            )

            truth_gt = gt_to_dosage(true_gt_string)
            imp_gt = gt_to_dosage(imp_gt_string)

            ds = parse_float(ds_string)
            gp = parse_gp(gp_string)

            # ================================================
            # DS-based R2
            # ================================================

            if truth_gt is not None and ds is not None:

                pooled_corr.add(ds, truth_gt)

                valid_ds_gt_pairs += 1

                site_ds.append(ds)
                site_truth.append(float(truth_gt))

            # ================================================
            # GP-weighted concordance
            # ================================================

            if truth_gt is not None and gp is not None:

                gp_total = sum(gp)

                if gp_total > 0:

                    gp_numerator += gp[int(truth_gt)]
                    gp_denominator += gp_total

            # ================================================
            # Hard-call concordance + NRC
            # ================================================

            if truth_gt is not None and imp_gt is not None:

                hard_total += 1

                if truth_gt == imp_gt:

                    hard_match += 1

                    if truth_gt == 0:
                        m0 += 1
                    elif truth_gt == 1:
                        m1 += 1
                    elif truth_gt == 2:
                        m2 += 1

                else:

                    if truth_gt == 0:
                        e0 += 1
                    elif truth_gt == 1:
                        e1 += 1
                    elif truth_gt == 2:
                        e2 += 1

        # ----------------------------------------------------
        # Variant-wise dosage R2
        # ----------------------------------------------------

        r2 = pearson_r2(site_ds, site_truth)

        if r2 is None:

            variant_r2_missing += 1

        else:

            variant_r2_sum += r2
            variant_r2_n += 1

        # ----------------------------------------------------
        # Variant GP concordance
        # ----------------------------------------------------

        if gp_denominator > 0:

            gp_conc_sum += (
                gp_numerator / gp_denominator
            )

            gp_conc_n += 1

        # ----------------------------------------------------
        # Move both streams
        # ----------------------------------------------------

        true_line = true_f.readline() or None
        imp_line = imp_f.readline() or None

        if (
            progress_every > 0
            and matched_positions % progress_every == 0
        ):

            print(
                f"Matched={matched_positions:,} | "
                f"Excluded scaffold={scaffold_excluded:,} | "
                f"Evaluated={evaluated_variants:,}"
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

    # ========================================================
    # Final metrics
    # ========================================================

    pooled_r2 = pooled_corr.r2()

    mean_variant_r2 = (
        variant_r2_sum / variant_r2_n
        if variant_r2_n > 0
        else None
    )

    mean_gp_concordance = (
        gp_conc_sum / gp_conc_n
        if gp_conc_n > 0
        else None
    )

    hard_concordance = (
        hard_match / hard_total
        if hard_total > 0
        else None
    )

    #
    # NRC:
    #
    # Correct non-reference calls:
    #     m1 + m2
    #
    # Denominator:
    #     all genotype errors + correct non-reference calls
    #
    nrc_denominator = (
        m1 + m2 + e0 + e1 + e2
    )

    nrc = (
        (m1 + m2) / nrc_denominator
        if nrc_denominator > 0
        else None
    )

    # ========================================================
    # Print results
    # ========================================================

    print("\n" + "=" * 78)
    print("IMPUTED-ONLY RESULTS")
    print("=" * 78)

    print(f"Matched positions:                    {matched_positions:,}")
    print(f"Allele mismatches:                    {allele_mismatch:,}")
    print(f"Input scaffold SNPs excluded:         {scaffold_excluded:,}")
    print(f"Variants evaluated:                   {evaluated_variants:,}")

    print(f"Truth-only positions:                 {truth_only:,}")
    print(f"Imputed-only positions:               {imputed_only:,}")

    print(f"Valid DS-GT pairs:                    {valid_ds_gt_pairs:,}")

    print(
        "Pooled empirical R2:                 "
        + (f"{pooled_r2:.10f}" if pooled_r2 is not None else "NA")
    )

    print(f"Variants with calculable dosage R2:   {variant_r2_n:,}")
    print(f"Variants without calculable R2:       {variant_r2_missing:,}")

    print(
        "Mean per-variant dosage R2:          "
        + (
            f"{mean_variant_r2:.10f}"
            if mean_variant_r2 is not None
            else "NA"
        )
    )

    print(
        "Mean GP-weighted concordance:         "
        + (
            f"{mean_gp_concordance:.10f}"
            if mean_gp_concordance is not None
            else "NA"
        )
    )

    print(
        "Hard genotype concordance:           "
        + (
            f"{hard_concordance:.10f}"
            if hard_concordance is not None
            else "NA"
        )
    )

    print(
        "Non-reference concordance (NRC):     "
        + (
            f"{nrc:.10f}"
            if nrc is not None
            else "NA"
        )
    )

    print()
    print("NRC components:")
    print(f"  Correct 0/0 (m0; NOT in NRC numerator): {m0:,}")
    print(f"  Correct 0/1 (m1):                       {m1:,}")
    print(f"  Correct 1/1 (m2):                       {m2:,}")
    print(f"  Errors for true 0/0 (e0):               {e0:,}")
    print(f"  Errors for true 0/1 (e1):               {e1:,}")
    print(f"  Errors for true 1/1 (e2):               {e2:,}")

    print("=" * 78)

    # ========================================================
    # Save output
    # ========================================================

    results = [
        ("Matched_positions", matched_positions),
        ("Allele_mismatches", allele_mismatch),
        ("Input_scaffold_SNPs_total", len(exclude_sites)),
        ("Input_scaffold_SNPs_excluded", scaffold_excluded),
        ("Imputed_variants_evaluated", evaluated_variants),
        ("Truth_only_positions", truth_only),
        ("Imputed_only_positions", imputed_only),
        ("Valid_DS_GT_pairs", valid_ds_gt_pairs),

        ("Pooled_empirical_R2",
         pooled_r2 if pooled_r2 is not None else "NA"),

        ("Variants_with_calculable_dosage_R2",
         variant_r2_n),

        ("Variants_without_calculable_dosage_R2",
         variant_r2_missing),

        ("Mean_per_variant_dosage_R2",
         mean_variant_r2 if mean_variant_r2 is not None else "NA"),

        ("Mean_GP_weighted_concordance",
         mean_gp_concordance
         if mean_gp_concordance is not None
         else "NA"),

        ("Hard_genotype_concordance",
         hard_concordance
         if hard_concordance is not None
         else "NA"),

        ("Non_reference_concordance",
         nrc if nrc is not None else "NA"),

        ("NRC_correct_00_m0", m0),
        ("NRC_correct_01_m1", m1),
        ("NRC_correct_11_m2", m2),

        ("NRC_error_true_00_e0", e0),
        ("NRC_error_true_01_e1", e1),
        ("NRC_error_true_11_e2", e2),
    ]

    with open(output_file, "w") as out:

        out.write("Metric\tValue\n")

        for metric, value in results:
            out.write(f"{metric}\t{value}\n")

    print(f"\nResults saved to: {output_file}")


# ============================================================
# CLI
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Evaluate imputation accuracy while excluding "
            "observed target/scaffold SNPs."
        )
    )

    parser.add_argument(
        "-t", "--truth",
        required=True,
        help="Complete WGS truth VCF.gz"
    )

    parser.add_argument(
        "-i", "--imputed",
        required=True,
        help="IMPUTE5 output VCF.gz"
    )

    parser.add_argument(
        "-x", "--exclude-sites",
        required=True,
        help="CHROM POS file containing observed target SNPs"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output summary TSV"
    )

    parser.add_argument(
        "--progress",
        type=int,
        default=100000,
        help="Progress interval"
    )

    args = parser.parse_args()

    evaluate(
        args.truth,
        args.imputed,
        args.exclude_sites,
        args.output,
        args.progress
    )


if __name__ == "__main__":
    main()
