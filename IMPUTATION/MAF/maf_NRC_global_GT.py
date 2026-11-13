#!/usr/bin/env python3

"""
Fast MAF-stratified NRC calculator for IMPUTE5 output.

Key features
------------
1. Only variants carrying the IMP flag in the imputed VCF are analyzed.
2. MAF is obtained from the truth VCF INFO field.
3. Imputed hard genotype calls are read directly from the IMPUTE5 GT field.
4. Raw and INFO >= threshold NRC are calculated simultaneously.
5. Supports indexed genomic regions for parallel chromosome analysis.
6. NRC is calculated from pooled genotype comparisons:

       NRC = (N11 + N22) / (Ntotal - N00)

   where rows are true genotypes and columns are imputed hard calls.

Genotype coding
---------------
0 = homozygous reference
1 = heterozygous
2 = homozygous alternate
"""

import argparse
import bisect
import math
import subprocess
import sys


# ============================================================
# MAF bins
#
# Same scheme used in Figure 2B:
#   0.01-wide bins from 0.00 to 0.10
#   0.05-wide bins from 0.10 onward
#
# Final upper boundary is set to 0.51 so that MAF = 0.50
# is included in the last interval.
# ============================================================

MAF_BINS = [
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
    (0.45, 0.51),
]

BIN_UPPERS = [high for _, high in MAF_BINS]

BIN_LABELS = [
    f"{low:.3f}-{high:.3f}"
    for low, high in MAF_BINS
]


# ============================================================
# Genotype coding
# ============================================================

GT_TO_DOSAGE = {
    "0/0": 0,
    "0|0": 0,

    "0/1": 1,
    "1/0": 1,
    "0|1": 1,
    "1|0": 1,

    "1/1": 2,
    "1|1": 2,
}


# ============================================================
# Sample IDs
# ============================================================

def get_samples(vcf_path):
    """
    Read sample IDs from an indexed VCF/BCF using bcftools.
    """

    cmd = [
        "bcftools",
        "query",
        "-l",
        vcf_path,
    ]

    try:
        result = subprocess.check_output(
            cmd,
            text=True,
        )

    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Failed to read sample IDs from {vcf_path}"
        ) from exc

    return [
        x
        for x in result.splitlines()
        if x
    ]


# ============================================================
# Open an indexed genomic region
# ============================================================

def open_region(vcf_path, region=None):
    """
    Open VCF records through bcftools view.

    Using -r requires an index and allows each SLURM array task
    to read only its own genomic interval.
    """

    cmd = [
        "bcftools",
        "view",
        "-H",
    ]

    if region is not None:
        cmd.extend([
            "-r",
            region,
        ])

    cmd.append(vcf_path)

    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=None,
        text=True,
        bufsize=1024 * 1024,
    )


# ============================================================
# INFO parsing
# ============================================================

def info_value(info_str, key):
    """
    Retrieve KEY=value from the VCF INFO string.

    Returns the value as text, or None when absent.
    """

    needle = key + "="
    start = 0

    while True:

        idx = info_str.find(
            needle,
            start,
        )

        if idx < 0:
            return None

        # Make sure the match starts at an INFO-field boundary
        if idx == 0 or info_str[idx - 1] == ";":

            value_start = idx + len(needle)

            value_end = info_str.find(
                ";",
                value_start,
            )

            if value_end < 0:
                value_end = len(info_str)

            return info_str[
                value_start:value_end
            ]

        start = idx + 1


def has_imp_tag(info_str):
    """
    Test whether INFO contains the flag IMP.
    """

    return (
        info_str == "IMP"
        or info_str.startswith("IMP;")
        or info_str.endswith(";IMP")
        or ";IMP;" in info_str
    )


def extract_maf(info_str):
    """
    Extract MAF from truth VCF INFO.

    Priority:
      1. MAF=
      2. AF= converted to minor allele frequency
    """

    value = info_value(
        info_str,
        "MAF",
    )

    if value is not None:

        try:
            maf = float(
                value.split(",", 1)[0]
            )

            if 0.0 <= maf <= 0.5:
                return maf

        except ValueError:
            pass

    value = info_value(
        info_str,
        "AF",
    )

    if value is not None:

        try:
            af = float(
                value.split(",", 1)[0]
            )

            maf = min(
                af,
                1.0 - af,
            )

            if 0.0 <= maf <= 0.5:
                return maf

        except ValueError:
            pass

    return None


def extract_info_score(info_str):
    """
    Extract IMPUTE5 INFO score from INFO.
    """

    value = info_value(
        info_str,
        "INFO",
    )

    if value is None:
        return None

    try:
        return float(
            value.split(",", 1)[0]
        )

    except ValueError:
        return None


# ============================================================
# Genotype parsing
# ============================================================

def gt_index(format_str):
    """
    Find GT position in FORMAT.

    Fast path handles GT:DS:GP, which is the observed
    IMPUTE5 FORMAT in this dataset.
    """

    if (
        format_str == "GT"
        or format_str.startswith("GT:")
    ):
        return 0

    fields = format_str.split(":")

    try:
        return fields.index("GT")

    except ValueError:
        return None


def extract_gt_dosage(sample_field, index):
    """
    Convert sample GT into genotype dosage 0/1/2.

    Both phased and unphased diploid genotypes are supported.
    """

    if index is None:
        return None

    # Fast path: GT is first FORMAT field
    if index == 0:

        gt = sample_field.partition(":")[0]

    else:

        parts = sample_field.split(":")

        if index >= len(parts):
            return None

        gt = parts[index]

    return GT_TO_DOSAGE.get(gt)


# ============================================================
# MAF bin assignment
# ============================================================

def get_bin_index(maf):
    """
    Return index of the MAF interval containing maf.
    """

    idx = bisect.bisect_right(
        BIN_UPPERS,
        maf,
    )

    if idx >= len(MAF_BINS):
        return None

    low, high = MAF_BINS[idx]

    if low <= maf < high:
        return idx

    return None


# ============================================================
# NRC calculations
#
# Matrix is stored as a flattened list:
#
#              imputed
#             0    1    2
#
# truth  0   N00  N01  N02
#        1   N10  N11  N12
#        2   N20  N21  N22
#
# flattened:
#
# [N00, N01, N02,
#  N10, N11, N12,
#  N20, N21, N22]
#
# NRC = (N11 + N22) / (Ntotal - N00)
# ============================================================

def nrc_stats(matrix):
    """
    Calculate NRC numerator, denominator and NRC.
    """

    total = sum(matrix)

    n00 = matrix[0]

    correct_nonref = (
        matrix[4]
        +
        matrix[8]
    )

    evaluable = (
        total
        -
        n00
    )

    if evaluable > 0:

        nrc = (
            correct_nonref
            /
            evaluable
        )

    else:

        nrc = math.nan

    return (
        total,
        correct_nonref,
        evaluable,
        nrc,
    )


def format_float(value):

    if math.isnan(value):
        return "NA"

    return f"{value:.10f}"


# ============================================================
# Main analysis
# ============================================================

def run_analysis(
    true_vcf,
    imputed_vcf,
    region,
    info_threshold,
    verbose=True,
):

    # --------------------------------------------------------
    # Sample-order QC
    # --------------------------------------------------------

    true_samples = get_samples(
        true_vcf
    )

    imputed_samples = get_samples(
        imputed_vcf
    )

    if true_samples != imputed_samples:

        print(
            "ERROR: sample IDs or sample order differ "
            "between truth and imputed VCFs.",
            file=sys.stderr,
        )

        sys.exit(2)

    print(
        f"Samples: {len(true_samples)}",
        flush=True,
    )

    print(
        f"Region: "
        f"{region if region else 'whole file'}",
        flush=True,
    )

    print(
        f"INFO threshold: >= {info_threshold}",
        flush=True,
    )

    print(
        "Hard calls: imputed GT field",
        flush=True,
    )

    print(
        "Site set: IMP-tagged imputed sites only",
        flush=True,
    )

    # --------------------------------------------------------
    # Initialize matrices
    # --------------------------------------------------------

    raw_matrices = [
        [0] * 9
        for _ in MAF_BINS
    ]

    filtered_matrices = [
        [0] * 9
        for _ in MAF_BINS
    ]

    raw_snps = [
        0
        for _ in MAF_BINS
    ]

    filtered_snps = [
        0
        for _ in MAF_BINS
    ]

    raw_observations = [
        0
        for _ in MAF_BINS
    ]

    filtered_observations = [
        0
        for _ in MAF_BINS
    ]

    # --------------------------------------------------------
    # QC counters
    # --------------------------------------------------------

    matched = 0

    unmatched_truth = 0

    unmatched_imputed = 0

    no_imp = 0

    no_maf = 0

    no_info = 0

    allele_mismatch_same_pos = 0

    # --------------------------------------------------------
    # Open indexed streams
    # --------------------------------------------------------

    true_proc = open_region(
        true_vcf,
        region,
    )

    imputed_proc = open_region(
        imputed_vcf,
        region,
    )

    true_line = (
        true_proc.stdout.readline()
    )

    imputed_line = (
        imputed_proc.stdout.readline()
    )

    # --------------------------------------------------------
    # Merge streams by chromosome and position
    # --------------------------------------------------------

    while true_line and imputed_line:

        true_cols = (
            true_line
            .rstrip("\n")
            .split("\t")
        )

        imputed_cols = (
            imputed_line
            .rstrip("\n")
            .split("\t")
        )

        true_chrom = true_cols[0]
        imputed_chrom = imputed_cols[0]

        true_pos = int(
            true_cols[1]
        )

        imputed_pos = int(
            imputed_cols[1]
        )

        true_key = (
            true_chrom,
            true_pos,
        )

        imputed_key = (
            imputed_chrom,
            imputed_pos,
        )

        # ----------------------------------------------------
        # Truth site occurs first
        # ----------------------------------------------------

        if true_key < imputed_key:

            unmatched_truth += 1

            true_line = (
                true_proc.stdout.readline()
            )

            continue

        # ----------------------------------------------------
        # Imputed site occurs first
        # ----------------------------------------------------

        if true_key > imputed_key:

            unmatched_imputed += 1

            imputed_line = (
                imputed_proc.stdout.readline()
            )

            continue

        # ----------------------------------------------------
        # Same CHROM:POS
        # ----------------------------------------------------

        matched += 1

        # REF/ALT mismatch is counted for QC.
        # Historical Figure 2B matching was CHROM + POS.
        if (
            len(true_cols) > 4
            and
            len(imputed_cols) > 4
            and
            (
                true_cols[3] != imputed_cols[3]
                or
                true_cols[4] != imputed_cols[4]
            )
        ):

            allele_mismatch_same_pos += 1

        # ----------------------------------------------------
        # Require IMP flag
        # ----------------------------------------------------

        imputed_info = (
            imputed_cols[7]
            if len(imputed_cols) > 7
            else ""
        )

        if not has_imp_tag(
            imputed_info
        ):

            no_imp += 1

            true_line = (
                true_proc.stdout.readline()
            )

            imputed_line = (
                imputed_proc.stdout.readline()
            )

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

        if maf is None:

            no_maf += 1

            true_line = (
                true_proc.stdout.readline()
            )

            imputed_line = (
                imputed_proc.stdout.readline()
            )

            continue

        bin_idx = get_bin_index(
            maf
        )

        if bin_idx is None:

            true_line = (
                true_proc.stdout.readline()
            )

            imputed_line = (
                imputed_proc.stdout.readline()
            )

            continue

        # ----------------------------------------------------
        # Locate GT in FORMAT
        # ----------------------------------------------------

        true_gt_idx = gt_index(
            true_cols[8]
        )

        imputed_gt_idx = gt_index(
            imputed_cols[8]
        )

        if (
            true_gt_idx is None
            or
            imputed_gt_idx is None
        ):

            true_line = (
                true_proc.stdout.readline()
            )

            imputed_line = (
                imputed_proc.stdout.readline()
            )

            continue

        # ----------------------------------------------------
        # INFO filtering
        # ----------------------------------------------------

        info_score = extract_info_score(
            imputed_info
        )

        passes_info = (
            info_score is not None
            and
            info_score >= info_threshold
        )

        if info_score is None:
            no_info += 1

        # ----------------------------------------------------
        # References to current bin matrices
        # ----------------------------------------------------

        raw_matrix = (
            raw_matrices[bin_idx]
        )

        filtered_matrix = (
            filtered_matrices[bin_idx]
        )

        valid_observations = 0

        # ----------------------------------------------------
        # Process every individual exactly once
        # ----------------------------------------------------

        for (
            truth_sample,
            imputed_sample,
        ) in zip(
            true_cols[9:],
            imputed_cols[9:],
        ):

            true_gt = extract_gt_dosage(
                truth_sample,
                true_gt_idx,
            )

            if true_gt is None:
                continue

            # Use IMPUTE5 best-guess GT directly.
            pred_gt = extract_gt_dosage(
                imputed_sample,
                imputed_gt_idx,
            )

            if pred_gt is None:
                continue

            # Flattened confusion-matrix cell
            cell = (
                true_gt * 3
                +
                pred_gt
            )

            # Raw NRC matrix
            raw_matrix[
                cell
            ] += 1

            # INFO-filtered NRC matrix
            if passes_info:

                filtered_matrix[
                    cell
                ] += 1

            valid_observations += 1

        # ----------------------------------------------------
        # SNP and observation counts
        # ----------------------------------------------------

        if valid_observations > 0:

            raw_snps[
                bin_idx
            ] += 1

            raw_observations[
                bin_idx
            ] += valid_observations

            if passes_info:

                filtered_snps[
                    bin_idx
                ] += 1

                filtered_observations[
                    bin_idx
                ] += valid_observations

        # ----------------------------------------------------
        # Advance both streams
        # ----------------------------------------------------

        true_line = (
            true_proc.stdout.readline()
        )

        imputed_line = (
            imputed_proc.stdout.readline()
        )

        # ----------------------------------------------------
        # Progress reporting
        # ----------------------------------------------------

        if (
            verbose
            and
            matched % 100000 == 0
        ):

            print(
                f"matched={matched:,} "
                f"raw_snps={sum(raw_snps):,} "
                f"filtered_snps={sum(filtered_snps):,}",
                flush=True,
            )

    # ========================================================
    # IMPORTANT:
    #
    # Drain remaining bcftools output before closing pipes.
    #
    # One VCF can reach EOF slightly before the other because
    # not every truth record must occur in the imputed VCF.
    #
    # Closing a still-producing bcftools stdout pipe directly
    # can cause SIGPIPE and return code -13.
    #
    # Draining both streams allows bcftools to terminate
    # normally.
    # ========================================================

    if true_proc.stdout is not None:

        for _ in true_proc.stdout:
            pass

        true_proc.stdout.close()

    if imputed_proc.stdout is not None:

        for _ in imputed_proc.stdout:
            pass

        imputed_proc.stdout.close()

    # --------------------------------------------------------
    # Wait for bcftools processes
    # --------------------------------------------------------

    true_return_code = (
        true_proc.wait()
    )

    imputed_return_code = (
        imputed_proc.wait()
    )

    # --------------------------------------------------------
    # Confirm successful bcftools termination
    # --------------------------------------------------------

    if (
        true_return_code != 0
        or
        imputed_return_code != 0
    ):

        raise RuntimeError(
            "bcftools failed: "
            f"truth rc={true_return_code}, "
            f"imputed rc={imputed_return_code}"
        )

    # --------------------------------------------------------
    # QC summary
    # --------------------------------------------------------

    qc = {
        "matched":
            matched,

        "unmatched_truth":
            unmatched_truth,

        "unmatched_imputed":
            unmatched_imputed,

        "no_imp":
            no_imp,

        "no_maf":
            no_maf,

        "no_info":
            no_info,

        "allele_mismatch_same_pos":
            allele_mismatch_same_pos,
    }

    return (
        raw_matrices,
        filtered_matrices,
        raw_snps,
        filtered_snps,
        raw_observations,
        filtered_observations,
        qc,
    )


# ============================================================
# Output
# ============================================================

def write_outputs(
    output,
    matrix_output,
    raw_matrices,
    filtered_matrices,
    raw_snps,
    filtered_snps,
    raw_observations,
    filtered_observations,
    region,
    info_threshold,
    qc,
):

    # --------------------------------------------------------
    # NRC summary
    # --------------------------------------------------------

    with open(
        output,
        "w",
    ) as out:

        out.write(
            "# MAF-Stratified NRC Results\n"
        )

        out.write(
            "# Hard call source: imputed GT field\n"
        )

        out.write(
            "# Site set: IMP-tagged imputed sites only\n"
        )

        out.write(
            f"# Region: "
            f"{region if region else 'whole file'}\n"
        )

        out.write(
            f"# Filtered definition: "
            f"INFO>={info_threshold}\n"
        )

        out.write(
            "# NRC=(N11+N22)/(Ntotal-N00)\n"
        )

        out.write(
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

        for (
            bin_idx,
            label,
        ) in enumerate(
            BIN_LABELS
        ):

            (
                raw_total,
                raw_correct_nonref,
                raw_evaluable,
                raw_nrc,
            ) = nrc_stats(
                raw_matrices[
                    bin_idx
                ]
            )

            (
                filtered_total,
                filtered_correct_nonref,
                filtered_evaluable,
                filtered_nrc,
            ) = nrc_stats(
                filtered_matrices[
                    bin_idx
                ]
            )

            # ------------------------------------------------
            # Internal consistency check
            # ------------------------------------------------

            if (
                raw_total
                !=
                raw_observations[
                    bin_idx
                ]
            ):

                raise RuntimeError(
                    "Raw observation-count mismatch "
                    f"in MAF bin {label}: "
                    f"matrix={raw_total}, "
                    f"counter="
                    f"{raw_observations[bin_idx]}"
                )

            if (
                filtered_total
                !=
                filtered_observations[
                    bin_idx
                ]
            ):

                raise RuntimeError(
                    "Filtered observation-count mismatch "
                    f"in MAF bin {label}: "
                    f"matrix={filtered_total}, "
                    f"counter="
                    f"{filtered_observations[bin_idx]}"
                )

            # ------------------------------------------------
            # Write result row
            # ------------------------------------------------

            out.write(
                f"{label}\t"
                f"{raw_snps[bin_idx]}\t"
                f"{raw_observations[bin_idx]}\t"
                f"{raw_correct_nonref}\t"
                f"{raw_evaluable}\t"
                f"{format_float(raw_nrc)}\t"
                f"{filtered_snps[bin_idx]}\t"
                f"{filtered_observations[bin_idx]}\t"
                f"{filtered_correct_nonref}\t"
                f"{filtered_evaluable}\t"
                f"{format_float(filtered_nrc)}\n"
            )

    # --------------------------------------------------------
    # Full 3x3 confusion matrices
    # --------------------------------------------------------

    with open(
        matrix_output,
        "w",
    ) as out:

        out.write(
            "MAF_Bin\t"
            "Mode\t"
            "N00\t"
            "N01\t"
            "N02\t"
            "N10\t"
            "N11\t"
            "N12\t"
            "N20\t"
            "N21\t"
            "N22\n"
        )

        for (
            bin_idx,
            label,
        ) in enumerate(
            BIN_LABELS
        ):

            out.write(
                label
                +
                "\tRAW\t"
                +
                "\t".join(
                    str(x)
                    for x in raw_matrices[
                        bin_idx
                    ]
                )
                +
                "\n"
            )

            out.write(
                label
                +
                "\tFILTERED\t"
                +
                "\t".join(
                    str(x)
                    for x in filtered_matrices[
                        bin_idx
                    ]
                )
                +
                "\n"
            )

    # --------------------------------------------------------
    # Console QC report
    # --------------------------------------------------------

    print(
        "\nQC:",
        flush=True,
    )

    for key, value in qc.items():

        print(
            f"  {key}: {value:,}",
            flush=True,
        )

    print(
        f"\nSummary written to: {output}",
        flush=True,
    )

    print(
        f"Matrices written to: {matrix_output}",
        flush=True,
    )


# ============================================================
# Command-line interface
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Fast MAF-stratified NRC calculation "
            "using IMPUTE5 GT hard calls."
        )
    )

    parser.add_argument(
        "-t",
        "--true-vcf",
        required=True,
        help=(
            "Truth VCF containing target-population "
            "MAF or AF in INFO."
        ),
    )

    parser.add_argument(
        "-i",
        "--imputed-vcf",
        required=True,
        help=(
            "IMPUTE5 output VCF containing "
            "GT, INFO and IMP."
        ),
    )

    parser.add_argument(
        "-r",
        "--region",
        default=None,
        help=(
            "Indexed genomic region, "
            "for example 11:1-119321753."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help=(
            "Output NRC summary TSV."
        ),
    )

    parser.add_argument(
        "--matrix-output",
        required=True,
        help=(
            "Output file containing pooled "
            "3x3 genotype confusion matrices."
        ),
    )

    parser.add_argument(
        "--info-threshold",
        type=float,
        default=0.8,
        help=(
            "INFO threshold for filtered NRC "
            "(default: 0.8; criterion is >= threshold)."
        ),
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help=(
            "Suppress progress messages."
        ),
    )

    args = parser.parse_args()

    (
        raw_matrices,
        filtered_matrices,
        raw_snps,
        filtered_snps,
        raw_observations,
        filtered_observations,
        qc,
    ) = run_analysis(
        true_vcf=args.true_vcf,
        imputed_vcf=args.imputed_vcf,
        region=args.region,
        info_threshold=args.info_threshold,
        verbose=not args.quiet,
    )

    write_outputs(
        output=args.output,
        matrix_output=args.matrix_output,
        raw_matrices=raw_matrices,
        filtered_matrices=filtered_matrices,
        raw_snps=raw_snps,
        filtered_snps=filtered_snps,
        raw_observations=raw_observations,
        filtered_observations=filtered_observations,
        region=args.region,
        info_threshold=args.info_threshold,
        qc=qc,
    )


if __name__ == "__main__":
    main()
