#!/usr/bin/env python3

"""
MAF-stratified Non-Reference Concordance (NRC) calculator.

Design:
1. Match truth and imputed VCF sites.
2. Retain only sites with the IMP flag in the imputed VCF.
3. Obtain MAF from the truth VCF INFO field.
4. Convert imputed GP to hard genotype calls using argmax(GP).
5. Calculate NRC within each MAF bin:

       NRC = (N11 + N22) / (Ntotal - N00)

   where rows are true genotypes and columns are imputed hard genotypes.

6. Calculate both:
   - Raw IMP-only NRC
   - IMP-only + INFO >= 0.8 NRC

The MAF bins are identical to those used for the existing Figure 2B analysis.
"""

import gzip
import argparse
import sys
import numpy as np


# ============================================================
# Basic helper functions
# ============================================================

def normalize_chrom(chrom):
    """Remove a leading 'chr' or 'Chr'."""
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def chrom_sort_key(chrom):
    """
    Natural chromosome sorting key.
    Numeric chromosomes are ordered numerically.
    """
    chrom = normalize_chrom(chrom)

    try:
        return (0, int(chrom))
    except ValueError:
        return (1, chrom)


def site_key(cols):
    """
    Return a sortable genomic key:
        ((chromosome sort key), position)
    """
    chrom = cols[0]
    pos = int(cols[1])

    return (
        chrom_sort_key(chrom),
        pos
    )


def parse_info(info_str):
    """
    Parse VCF INFO field.

    Returns
    -------
    values : dict
        KEY=VALUE fields.
    flags : set
        Flag-only fields such as IMP.
    """

    values = {}
    flags = set()

    if not info_str or info_str == ".":
        return values, flags

    for field in info_str.split(";"):

        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
        else:
            flags.add(field)

    return values, flags


def extract_maf(info_str):
    """
    Extract target-population MAF from truth VCF INFO.

    Priority:
      1. MAF
      2. AF converted to MAF
    """

    values, _ = parse_info(info_str)

    if "MAF" in values:
        try:
            maf = float(values["MAF"].split(",")[0])

            if 0 <= maf <= 0.5:
                return maf
        except Exception:
            pass

    if "AF" in values:
        try:
            af = float(values["AF"].split(",")[0])
            maf = min(af, 1.0 - af)

            if 0 <= maf <= 0.5:
                return maf
        except Exception:
            pass

    return None


def has_imp_tag(info_str):
    """Return True when the imputed VCF INFO contains the IMP flag."""

    _, flags = parse_info(info_str)

    return "IMP" in flags


def extract_info_score(info_str):
    """Extract IMPUTE5 INFO score."""

    values, _ = parse_info(info_str)

    try:
        return float(values["INFO"].split(",")[0])
    except Exception:
        return np.nan


# ============================================================
# GT / GP parsing
# ============================================================

def parse_true_gt(sample_field, format_fields):
    """
    Parse truth GT and convert to 0/1/2 alternate-allele count.

    Supports both:
      0/0, 0/1, 1/0, 1/1
      0|0, 0|1, 1|0, 1|1

    Missing and non-biallelic genotypes are ignored.
    """

    try:
        gt_idx = format_fields.index("GT")
    except ValueError:
        return None

    parts = sample_field.split(":")

    if gt_idx >= len(parts):
        return None

    gt = parts[gt_idx]

    if gt in (".", "./.", ".|."):
        return None

    if "|" in gt:
        alleles = gt.split("|")
    elif "/" in gt:
        alleles = gt.split("/")
    else:
        return None

    if len(alleles) != 2:
        return None

    if "." in alleles:
        return None

    try:
        a1 = int(alleles[0])
        a2 = int(alleles[1])
    except ValueError:
        return None

    # This analysis assumes biallelic SNPs
    if a1 not in (0, 1) or a2 not in (0, 1):
        return None

    return a1 + a2


def parse_gp(sample_field, format_fields):
    """
    Extract GP from the imputed VCF according to FORMAT.

    Returns
    -------
    list of three floats or None
    """

    try:
        gp_idx = format_fields.index("GP")
    except ValueError:
        return None

    parts = sample_field.split(":")

    if gp_idx >= len(parts):
        return None

    gp_text = parts[gp_idx]

    if gp_text in (".", ""):
        return None

    try:
        gp = [float(x) for x in gp_text.split(",")]
    except Exception:
        return None

    if len(gp) != 3:
        return None

    if not all(np.isfinite(x) for x in gp):
        return None

    return gp


# ============================================================
# NRC accumulator
# ============================================================

class NRCAccumulator:

    def __init__(self):

        # rows    = true genotype
        # columns = imputed hard genotype
        self.matrix = np.zeros(
            (3, 3),
            dtype=np.int64
        )

        self.n_snps = 0
        self.n_observations = 0

    def add_snp(self, true_samples, imputed_samples,
                true_format, imputed_format):

        """
        Add all valid sample comparisons for one SNP.

        Returns number of valid genotype comparisons at this SNP.
        """

        valid = 0

        for true_sample, imp_sample in zip(
            true_samples,
            imputed_samples
        ):

            true_gt = parse_true_gt(
                true_sample,
                true_format
            )

            if true_gt is None:
                continue

            gp = parse_gp(
                imp_sample,
                imputed_format
            )

            if gp is None:
                continue

            # Hard genotype call for NRC
            pred_gt = int(np.argmax(gp))

            self.matrix[
                true_gt,
                pred_gt
            ] += 1

            valid += 1

        if valid > 0:
            self.n_snps += 1
            self.n_observations += valid

        return valid

    def get_nrc(self):

        """
        NRC =
          correctly imputed non-reference genotypes
          -----------------------------------------
          all comparisons except correct 0 -> 0
        """

        total = int(
            np.sum(self.matrix)
        )

        n00 = int(
            self.matrix[0, 0]
        )

        correct_nonref = int(
            self.matrix[1, 1]
            +
            self.matrix[2, 2]
        )

        evaluable = total - n00

        if evaluable == 0:
            nrc = np.nan
        else:
            nrc = (
                correct_nonref
                /
                evaluable
            )

        return {
            "nrc": nrc,
            "correct_nonref": correct_nonref,
            "evaluable": evaluable,
            "total_observations": total,
            "matrix": self.matrix
        }


# ============================================================
# MAF bins
# ============================================================

def get_default_maf_bins():

    return [
        (0.00, 0.01),
        (0.01, 0.02),
        (0.02, 0.03),
        (0.03, 0.04),
        (0.04, 0.05),
        (0.05, 0.06),
        (0.06, 0.07),
        (0.07, 0.08),
        (0.08, 0.09),
        (0.09, 0.10),
        (0.10, 0.15),
        (0.15, 0.20),
        (0.20, 0.25),
        (0.25, 0.30),
        (0.30, 0.35),
        (0.35, 0.40),
        (0.40, 0.45),
        (0.45, 0.51)
    ]


def maf_bin_key(maf, maf_bins):

    for low, high in maf_bins:

        if low <= maf < high:
            return f"{low:.3f}-{high:.3f}"

    return None


# ============================================================
# Open VCF and obtain samples + first record
# ============================================================

def open_vcf(path):

    f = gzip.open(
        path,
        "rt"
    )

    samples = None
    first_data_line = None

    for line in f:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):

            cols = line.rstrip().split("\t")
            samples = cols[9:]
            continue

        if not line.startswith("#"):

            first_data_line = line
            break

    if samples is None:
        f.close()
        raise RuntimeError(
            f"#CHROM header not found in {path}"
        )

    return (
        f,
        samples,
        first_data_line
    )


# ============================================================
# Main analysis
# ============================================================

def calculate_nrc(
    true_vcf,
    imputed_vcf,
    maf_bins,
    info_threshold=0.8,
    verbose=True
):

    print("=" * 72)
    print("MAF-Stratified NRC Calculator")
    print("IMP-only, raw + INFO-filtered modes")
    print("=" * 72)

    print(
        f"True VCF:       {true_vcf}"
    )

    print(
        f"Imputed VCF:    {imputed_vcf}"
    )

    print(
        f"INFO threshold: >= {info_threshold}"
    )

    print(
        "Only imputed sites carrying the IMP flag are analyzed."
    )

    print("=" * 72)

    # --------------------------------------------------------
    # Open VCFs
    # --------------------------------------------------------

    (
        true_f,
        true_samples,
        true_line
    ) = open_vcf(true_vcf)

    (
        imp_f,
        imp_samples,
        imp_line
    ) = open_vcf(imputed_vcf)

    print(
        f"\nTruth samples:   {len(true_samples)}"
    )

    print(
        f"Imputed samples: {len(imp_samples)}"
    )

    # Exact sample ID + order check
    if true_samples != imp_samples:

        true_f.close()
        imp_f.close()

        print(
            "\nERROR: Sample IDs or sample order differ "
            "between the two VCF files."
        )

        sys.exit(1)

    n_samples = len(
        true_samples
    )

    print(
        f"Sample order confirmed for {n_samples} individuals."
    )

    # --------------------------------------------------------
    # Initialize accumulators
    # --------------------------------------------------------

    raw = {}

    filtered = {}

    for low, high in maf_bins:

        key = f"{low:.3f}-{high:.3f}"

        raw[key] = NRCAccumulator()
        filtered[key] = NRCAccumulator()

    # --------------------------------------------------------
    # QC counters
    # --------------------------------------------------------

    matched = 0
    unmatched_steps = 0
    no_imp = 0
    no_maf = 0
    no_info = 0

    raw_valid_snps = 0
    filtered_valid_snps = 0

    # --------------------------------------------------------
    # Streaming site matching
    # --------------------------------------------------------

    while (
        true_line is not None
        and
        imp_line is not None
    ):

        true_cols = (
            true_line.rstrip().split("\t")
        )

        imp_cols = (
            imp_line.rstrip().split("\t")
        )

        true_key = site_key(
            true_cols
        )

        imp_key = site_key(
            imp_cols
        )

        if true_key < imp_key:

            true_line = true_f.readline()

            if true_line == "":
                true_line = None

            unmatched_steps += 1
            continue

        if true_key > imp_key:

            imp_line = imp_f.readline()

            if imp_line == "":
                imp_line = None

            unmatched_steps += 1
            continue

        # ----------------------------------------------------
        # Same genomic position
        # ----------------------------------------------------

        matched += 1

        imp_info = (
            imp_cols[7]
            if len(imp_cols) > 7
            else ""
        )

        # IMP-only
        if not has_imp_tag(
            imp_info
        ):

            no_imp += 1

            true_line = true_f.readline()
            imp_line = imp_f.readline()

            if true_line == "":
                true_line = None

            if imp_line == "":
                imp_line = None

            continue

        # ----------------------------------------------------
        # MAF from truth VCF
        # ----------------------------------------------------

        truth_info = (
            true_cols[7]
            if len(true_cols) > 7
            else ""
        )

        maf = extract_maf(
            truth_info
        )

        if maf is None or not np.isfinite(maf):

            no_maf += 1

            true_line = true_f.readline()
            imp_line = imp_f.readline()

            if true_line == "":
                true_line = None

            if imp_line == "":
                imp_line = None

            continue

        bin_key = maf_bin_key(
            maf,
            maf_bins
        )

        if bin_key is None:

            true_line = true_f.readline()
            imp_line = imp_f.readline()

            if true_line == "":
                true_line = None

            if imp_line == "":
                imp_line = None

            continue

        # ----------------------------------------------------
        # FORMAT fields
        # ----------------------------------------------------

        true_format = (
            true_cols[8].split(":")
        )

        imp_format = (
            imp_cols[8].split(":")
        )

        # ----------------------------------------------------
        # RAW NRC
        # ----------------------------------------------------

        n_raw = raw[
            bin_key
        ].add_snp(
            true_cols[9:],
            imp_cols[9:],
            true_format,
            imp_format
        )

        if n_raw > 0:
            raw_valid_snps += 1

        # ----------------------------------------------------
        # INFO >= threshold NRC
        # ----------------------------------------------------

        info_score = extract_info_score(
            imp_info
        )

        if not np.isfinite(info_score):

            no_info += 1

        elif info_score >= info_threshold:

            n_filtered = filtered[
                bin_key
            ].add_snp(
                true_cols[9:],
                imp_cols[9:],
                true_format,
                imp_format
            )

            if n_filtered > 0:
                filtered_valid_snps += 1

        # ----------------------------------------------------
        # Advance both files
        # ----------------------------------------------------

        true_line = true_f.readline()
        imp_line = imp_f.readline()

        if true_line == "":
            true_line = None

        if imp_line == "":
            imp_line = None

        if (
            verbose
            and
            matched % 100000 == 0
        ):

            print(
                f"Matched={matched:,}; "
                f"raw SNPs={raw_valid_snps:,}; "
                f"filtered SNPs={filtered_valid_snps:,}"
            )

    true_f.close()
    imp_f.close()

    print("\nProcessing complete.")

    print(
        f"Matched genomic positions: {matched:,}"
    )

    print(
        f"Unmatched streaming steps: {unmatched_steps:,}"
    )

    print(
        f"Matched sites without IMP: {no_imp:,}"
    )

    print(
        f"IMP sites without MAF:     {no_maf:,}"
    )

    print(
        f"IMP sites without INFO:    {no_info:,}"
    )

    # --------------------------------------------------------
    # Collect results
    # --------------------------------------------------------

    results = []

    for low, high in maf_bins:

        key = f"{low:.3f}-{high:.3f}"

        raw_stats = (
            raw[key].get_nrc()
        )

        filt_stats = (
            filtered[key].get_nrc()
        )

        results.append({
            "MAF_Bin": key,

            "Raw_SNPs":
                raw[key].n_snps,

            "Raw_Observations":
                raw_stats["total_observations"],

            "Raw_Correct_NonRef":
                raw_stats["correct_nonref"],

            "Raw_NRC_Evaluable":
                raw_stats["evaluable"],

            "Raw_NRC":
                raw_stats["nrc"],

            "Filtered_SNPs":
                filtered[key].n_snps,

            "Filtered_Observations":
                filt_stats["total_observations"],

            "Filtered_Correct_NonRef":
                filt_stats["correct_nonref"],

            "Filtered_NRC_Evaluable":
                filt_stats["evaluable"],

            "Filtered_NRC":
                filt_stats["nrc"],

            "Raw_Matrix":
                raw_stats["matrix"],

            "Filtered_Matrix":
                filt_stats["matrix"]
        })

    return results


# ============================================================
# Save output
# ============================================================

def save_summary(
    results,
    output
):

    with open(
        output,
        "w"
    ) as f:

        f.write(
            "# MAF-Stratified NRC Results\n"
        )

        f.write(
            "# Raw = IMP-only\n"
        )

        f.write(
            "# Filtered = IMP-only + INFO>=0.8\n"
        )

        f.write(
            "# NRC=(N11+N22)/(Ntotal-N00)\n"
        )

        f.write(
            "MAF_Bin\t"
            "Raw_SNPs\t"
            "Raw_Observations\t"
            "Raw_Correct_NonRef\t"
            "Raw_NRC_Evaluable\t"
            "Raw_NRC\t"
            "Filtered_SNPs\t"
            "Filtered_Observations\t"
            "Filtered_Correct_NonRef\t"
            "Filtered_NRC_Evaluable\t"
            "Filtered_NRC\n"
        )

        for res in results:

            raw_nrc = (
                "NA"
                if np.isnan(res["Raw_NRC"])
                else f'{res["Raw_NRC"]:.10f}'
            )

            filt_nrc = (
                "NA"
                if np.isnan(res["Filtered_NRC"])
                else f'{res["Filtered_NRC"]:.10f}'
            )

            f.write(
                f'{res["MAF_Bin"]}\t'
                f'{res["Raw_SNPs"]}\t'
                f'{res["Raw_Observations"]}\t'
                f'{res["Raw_Correct_NonRef"]}\t'
                f'{res["Raw_NRC_Evaluable"]}\t'
                f'{raw_nrc}\t'
                f'{res["Filtered_SNPs"]}\t'
                f'{res["Filtered_Observations"]}\t'
                f'{res["Filtered_Correct_NonRef"]}\t'
                f'{res["Filtered_NRC_Evaluable"]}\t'
                f'{filt_nrc}\n'
            )


def save_matrices(
    results,
    output
):

    with open(
        output,
        "w"
    ) as f:

        f.write(
            "# NRC hard-call confusion matrices\n"
        )

        f.write(
            "# Rows=true GT; columns=imputed MAP GT\n\n"
        )

        for res in results:

            f.write(
                f'## MAF {res["MAF_Bin"]}\n'
            )

            f.write(
                "RAW\n"
            )

            for row in res[
                "Raw_Matrix"
            ]:

                f.write(
                    "\t".join(
                        str(int(x))
                        for x in row
                    )
                    + "\n"
                )

            f.write(
                "FILTERED_INFO_GE_0.8\n"
            )

            for row in res[
                "Filtered_Matrix"
            ]:

                f.write(
                    "\t".join(
                        str(int(x))
                        for x in row
                    )
                    + "\n"
                )

            f.write("\n")


# ============================================================
# Command line
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate MAF-stratified NRC for "
            "IMPUTE5 imputed sites."
        )
    )

    parser.add_argument(
        "-t",
        "--true-vcf",
        required=True,
        help="Truth VCF containing MAF or AF"
    )

    parser.add_argument(
        "-i",
        "--imputed-vcf",
        required=True,
        help="Imputed VCF containing GP, INFO and IMP"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="results_MAF-NRC-all.txt",
        help="Summary output file"
    )

    parser.add_argument(
        "--matrix-output",
        default="results_MAF-NRC-matrices.txt",
        help="Hard-call confusion-matrix output"
    )

    parser.add_argument(
        "--info-threshold",
        type=float,
        default=0.8,
        help="INFO threshold; default 0.8"
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress progress messages"
    )

    args = parser.parse_args()

    maf_bins = (
        get_default_maf_bins()
    )

    results = calculate_nrc(
        args.true_vcf,
        args.imputed_vcf,
        maf_bins,
        info_threshold=args.info_threshold,
        verbose=not args.quiet
    )

    save_summary(
        results,
        args.output
    )

    save_matrices(
        results,
        args.matrix_output
    )

    print("\n==========================================")
    print("NRC analysis complete")
    print(f"Summary:  {args.output}")
    print(f"Matrices: {args.matrix_output}")
    print("==========================================")


if __name__ == "__main__":
    main()
