#!/usr/bin/env python3

"""
Unified MAF-stratified imputation accuracy analysis.

For each MAF bin, simultaneously calculates:

1. Dosage R²
   Pearson correlation squared between:
       truth GT dosage (0/1/2)
       imputed DS

2. Overall concordance
   GP-based probability-weighted concordance:
       sum(GP assigned to the true genotype)
       --------------------------------------
       sum(all GP probabilities)

3. Non-reference concordance (NRC)
   Based on IMPUTE5 best-guess GT:
       NRC = (N11 + N22) / (Ntotal - N00)

Only sites carrying the IMP flag are analyzed.

Both:
    RAW
and:
    INFO >= 0.8
are calculated simultaneously.

IMPORTANT:
Only sample-level observations with valid:
    truth GT + imputed DS + imputed GP + imputed GT
are used.

Thus R², concordance and NRC are calculated from the same
SNP-by-individual observations.

Supports indexed genomic regions for SLURM-array parallelization.
"""

import argparse
import bisect
import math
import subprocess
import sys


# ============================================================
# MAF bins
# Same binning used for Figure 2B
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

BIN_UPPERS = [
    high
    for _, high in MAF_BINS
]

BIN_LABELS = [
    f"{low:.3f}-{high:.3f}"
    for low, high in MAF_BINS
]


# ============================================================
# GT coding
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
# Basic VCF helpers
# ============================================================

def get_samples(vcf_path):

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


def open_region(vcf_path, region=None):

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

    cmd.append(
        vcf_path
    )

    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=None,
        text=True,
        bufsize=1024 * 1024,
    )


def normalize_chrom(chrom):

    chrom = chrom.strip()

    if chrom.lower().startswith("chr"):
        return chrom[3:]

    return chrom


# ============================================================
# INFO parsing
# ============================================================

def info_value(info_str, key):

    needle = key + "="

    start = 0

    while True:

        idx = info_str.find(
            needle,
            start,
        )

        if idx < 0:
            return None

        if (
            idx == 0
            or
            info_str[idx - 1] == ";"
        ):

            value_start = (
                idx
                +
                len(needle)
            )

            value_end = info_str.find(
                ";",
                value_start,
            )

            if value_end < 0:
                value_end = len(
                    info_str
                )

            return info_str[
                value_start:value_end
            ]

        start = idx + 1


def has_imp_tag(info_str):

    return (
        info_str == "IMP"
        or
        info_str.startswith("IMP;")
        or
        info_str.endswith(";IMP")
        or
        ";IMP;" in info_str
    )


def extract_maf(info_str):

    # Prefer explicit MAF
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

    # Fall back to AF
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
# FORMAT handling
# ============================================================

def get_format_indices(format_str):

    fields = format_str.split(":")

    index = {
        key: i
        for i, key in enumerate(fields)
    }

    return index


def parse_gt_from_parts(parts, index):

    if index is None:
        return None

    if index >= len(parts):
        return None

    return GT_TO_DOSAGE.get(
        parts[index]
    )


def parse_ds_from_parts(parts, index):

    if index is None:
        return None

    if index >= len(parts):
        return None

    value = parts[index]

    if (
        value == "."
        or
        value == ""
    ):
        return None

    try:

        ds = float(value)

    except ValueError:

        return None

    if not math.isfinite(ds):
        return None

    return ds


def parse_gp_from_parts(parts, index):

    if index is None:
        return None

    if index >= len(parts):
        return None

    value = parts[index]

    if (
        value == "."
        or
        value == ""
    ):
        return None

    fields = value.split(",")

    if len(fields) != 3:
        return None

    try:

        gp0 = float(fields[0])
        gp1 = float(fields[1])
        gp2 = float(fields[2])

    except ValueError:

        return None

    if not (
        math.isfinite(gp0)
        and
        math.isfinite(gp1)
        and
        math.isfinite(gp2)
    ):
        return None

    return (
        gp0,
        gp1,
        gp2,
    )


# ============================================================
# MAF bin assignment
# ============================================================

def get_bin_index(maf):

    idx = bisect.bisect_right(
        BIN_UPPERS,
        maf,
    )

    if idx >= len(MAF_BINS):
        return None

    low, high = (
        MAF_BINS[idx]
    )

    if low <= maf < high:
        return idx

    return None


# ============================================================
# Sufficient-statistics accumulator
# ============================================================

class AccuracyAccumulator:

    def __init__(self):

        self.snp_count = 0
        self.n = 0

        # R²
        self.sum_ds = 0.0
        self.sum_gt = 0.0
        self.sum_ds_gt = 0.0
        self.sum_ds2 = 0.0
        self.sum_gt2 = 0.0

        # GP-based concordance
        self.conc_numerator = 0.0
        self.conc_denominator = 0.0

        # NRC matrix
        #
        # N00 N01 N02
        # N10 N11 N12
        # N20 N21 N22
        self.nrc_matrix = [
            0
            for _ in range(9)
        ]

    def add_site(
        self,
        site_n,
        sum_ds,
        sum_gt,
        sum_ds_gt,
        sum_ds2,
        sum_gt2,
        conc_num,
        conc_den,
        nrc_matrix,
    ):

        if site_n <= 0:
            return

        self.snp_count += 1
        self.n += site_n

        self.sum_ds += sum_ds
        self.sum_gt += sum_gt
        self.sum_ds_gt += sum_ds_gt
        self.sum_ds2 += sum_ds2
        self.sum_gt2 += sum_gt2

        self.conc_numerator += conc_num
        self.conc_denominator += conc_den

        for j in range(9):

            self.nrc_matrix[j] += (
                nrc_matrix[j]
            )


# ============================================================
# Metric calculations
# ============================================================

def calculate_r2(acc):

    n = acc.n

    if n < 3:
        return math.nan

    numerator = (
        n * acc.sum_ds_gt
        -
        acc.sum_ds * acc.sum_gt
    )

    denominator_ds = (
        n * acc.sum_ds2
        -
        acc.sum_ds * acc.sum_ds
    )

    denominator_gt = (
        n * acc.sum_gt2
        -
        acc.sum_gt * acc.sum_gt
    )

    if (
        denominator_ds <= 0
        or
        denominator_gt <= 0
    ):
        return math.nan

    denominator = math.sqrt(
        denominator_ds
        *
        denominator_gt
    )

    if denominator == 0:
        return math.nan

    r = (
        numerator
        /
        denominator
    )

    return r * r


def calculate_concordance(acc):

    if acc.conc_denominator <= 0:
        return math.nan

    return (
        acc.conc_numerator
        /
        acc.conc_denominator
    )


def calculate_nrc(acc):

    matrix = (
        acc.nrc_matrix
    )

    total = sum(
        matrix
    )

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

    if evaluable <= 0:
        return math.nan

    return (
        correct_nonref
        /
        evaluable
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
    # Verify exact sample order
    # --------------------------------------------------------

    true_samples = get_samples(
        true_vcf
    )

    imputed_samples = get_samples(
        imputed_vcf
    )

    if true_samples != imputed_samples:

        raise RuntimeError(
            "Sample IDs/order differ between truth "
            "and imputed VCFs."
        )

    print(
        f"Samples: {len(true_samples)}",
        flush=True,
    )

    print(
        f"Region: {region if region else 'whole file'}",
        flush=True,
    )

    print(
        f"INFO threshold: >= {info_threshold}",
        flush=True,
    )

    print(
        "Site set: IMP-tagged sites only",
        flush=True,
    )

    print(
        "Metrics: dosage R2 + GP concordance + GT NRC",
        flush=True,
    )

    print(
        "All metrics use the same valid sample-level observations.",
        flush=True,
    )

    # --------------------------------------------------------
    # Accumulators
    # --------------------------------------------------------

    raw = [
        AccuracyAccumulator()
        for _ in MAF_BINS
    ]

    filtered = [
        AccuracyAccumulator()
        for _ in MAF_BINS
    ]

    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------

    matched = 0

    unmatched_truth = 0

    unmatched_imputed = 0

    no_imp = 0

    no_maf = 0

    no_info = 0

    allele_mismatch_same_pos = 0

    no_common_valid_observation = 0

    # --------------------------------------------------------
    # Indexed regional streams
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
    # Merge by genomic position
    # --------------------------------------------------------

    while (
        true_line
        and
        imputed_line
    ):

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

        true_chrom = normalize_chrom(
            true_cols[0]
        )

        imputed_chrom = normalize_chrom(
            imputed_cols[0]
        )

        true_pos = int(
            true_cols[1]
        )

        imputed_pos = int(
            imputed_cols[1]
        )

        # In regional analysis this should normally be same chr
        if true_chrom != imputed_chrom:

            true_key = (
                true_chrom,
                true_pos,
            )

            imputed_key = (
                imputed_chrom,
                imputed_pos,
            )

            if true_key < imputed_key:

                unmatched_truth += 1

                true_line = (
                    true_proc.stdout.readline()
                )

                continue

            else:

                unmatched_imputed += 1

                imputed_line = (
                    imputed_proc.stdout.readline()
                )

                continue

        # Same chromosome; compare numerical POS
        if true_pos < imputed_pos:

            unmatched_truth += 1

            true_line = (
                true_proc.stdout.readline()
            )

            continue

        if true_pos > imputed_pos:

            unmatched_imputed += 1

            imputed_line = (
                imputed_proc.stdout.readline()
            )

            continue

        # ----------------------------------------------------
        # Same chromosome and position
        # ----------------------------------------------------

        matched += 1

        # REF/ALT QC
        if (
            true_cols[3] != imputed_cols[3]
            or
            true_cols[4] != imputed_cols[4]
        ):

            allele_mismatch_same_pos += 1

        # ----------------------------------------------------
        # IMP-only
        # ----------------------------------------------------

        imputed_info = (
            imputed_cols[7]
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
        # MAF from truth
        # ----------------------------------------------------

        maf = extract_maf(
            true_cols[7]
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
        # FORMAT indices
        # ----------------------------------------------------

        true_format = get_format_indices(
            true_cols[8]
        )

        imp_format = get_format_indices(
            imputed_cols[8]
        )

        true_gt_idx = (
            true_format.get("GT")
        )

        imp_gt_idx = (
            imp_format.get("GT")
        )

        ds_idx = (
            imp_format.get("DS")
        )

        gp_idx = (
            imp_format.get("GP")
        )

        if (
            true_gt_idx is None
            or
            imp_gt_idx is None
            or
            ds_idx is None
            or
            gp_idx is None
        ):

            no_common_valid_observation += 1

            true_line = (
                true_proc.stdout.readline()
            )

            imputed_line = (
                imputed_proc.stdout.readline()
            )

            continue

        # ----------------------------------------------------
        # INFO filter
        # ----------------------------------------------------

        info_score = extract_info_score(
            imputed_info
        )

        if info_score is None:

            no_info += 1

            passes_info = False

        else:

            passes_info = (
                info_score
                >=
                info_threshold
            )

        # ----------------------------------------------------
        # Local per-site statistics
        # ----------------------------------------------------

        site_n = 0

        site_sum_ds = 0.0

        site_sum_gt = 0.0

        site_sum_ds_gt = 0.0

        site_sum_ds2 = 0.0

        site_sum_gt2 = 0.0

        site_conc_num = 0.0

        site_conc_den = 0.0

        site_nrc_matrix = [
            0
            for _ in range(9)
        ]

        # ----------------------------------------------------
        # Same observation set for all three metrics
        # ----------------------------------------------------

        for (
            truth_sample,
            imputed_sample,
        ) in zip(
            true_cols[9:],
            imputed_cols[9:],
        ):

            truth_parts = (
                truth_sample.split(":")
            )

            imputed_parts = (
                imputed_sample.split(":")
            )

            truth_gt = parse_gt_from_parts(
                truth_parts,
                true_gt_idx,
            )

            if truth_gt is None:
                continue

            imputed_gt = parse_gt_from_parts(
                imputed_parts,
                imp_gt_idx,
            )

            if imputed_gt is None:
                continue

            ds = parse_ds_from_parts(
                imputed_parts,
                ds_idx,
            )

            if ds is None:
                continue

            gp = parse_gp_from_parts(
                imputed_parts,
                gp_idx,
            )

            if gp is None:
                continue

            # ------------------------------------------------
            # R² sufficient statistics
            # ------------------------------------------------

            gt_float = float(
                truth_gt
            )

            site_sum_ds += ds

            site_sum_gt += (
                gt_float
            )

            site_sum_ds_gt += (
                ds
                *
                gt_float
            )

            site_sum_ds2 += (
                ds
                *
                ds
            )

            site_sum_gt2 += (
                gt_float
                *
                gt_float
            )

            # ------------------------------------------------
            # GP-based concordance
            #
            # Equivalent to the old probability-weighted
            # confusion-matrix diagonal / total.
            # ------------------------------------------------

            gp_sum = (
                gp[0]
                +
                gp[1]
                +
                gp[2]
            )

            if truth_gt == 0:

                site_conc_num += gp[0]

            elif truth_gt == 1:

                site_conc_num += gp[1]

            else:

                site_conc_num += gp[2]

            site_conc_den += gp_sum

            # ------------------------------------------------
            # GT-based NRC
            # ------------------------------------------------

            cell = (
                truth_gt * 3
                +
                imputed_gt
            )

            site_nrc_matrix[
                cell
            ] += 1

            site_n += 1

        # ----------------------------------------------------
        # Add site only when at least one common-valid
        # observation exists
        # ----------------------------------------------------

        if site_n > 0:

            raw[
                bin_idx
            ].add_site(
                site_n,
                site_sum_ds,
                site_sum_gt,
                site_sum_ds_gt,
                site_sum_ds2,
                site_sum_gt2,
                site_conc_num,
                site_conc_den,
                site_nrc_matrix,
            )

            if passes_info:

                filtered[
                    bin_idx
                ].add_site(
                    site_n,
                    site_sum_ds,
                    site_sum_gt,
                    site_sum_ds_gt,
                    site_sum_ds2,
                    site_sum_gt2,
                    site_conc_num,
                    site_conc_den,
                    site_nrc_matrix,
                )

        else:

            no_common_valid_observation += 1

        # ----------------------------------------------------
        # Advance
        # ----------------------------------------------------

        true_line = (
            true_proc.stdout.readline()
        )

        imputed_line = (
            imputed_proc.stdout.readline()
        )

        # ----------------------------------------------------
        # Progress
        # ----------------------------------------------------

        if (
            verbose
            and
            matched % 100000 == 0
        ):

            print(
                f"matched={matched:,} "
                f"raw_snps={sum(x.snp_count for x in raw):,} "
                f"filtered_snps="
                f"{sum(x.snp_count for x in filtered):,}",
                flush=True,
            )

    # ========================================================
    # Drain remaining pipes to prevent SIGPIPE
    # ========================================================

    if true_proc.stdout is not None:

        for _ in true_proc.stdout:
            pass

        true_proc.stdout.close()

    if imputed_proc.stdout is not None:

        for _ in imputed_proc.stdout:
            pass

        imputed_proc.stdout.close()

    true_rc = (
        true_proc.wait()
    )

    imputed_rc = (
        imputed_proc.wait()
    )

    if (
        true_rc != 0
        or
        imputed_rc != 0
    ):

        raise RuntimeError(
            "bcftools failed: "
            f"truth rc={true_rc}, "
            f"imputed rc={imputed_rc}"
        )

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

        "no_common_valid_observation":
            no_common_valid_observation,
    }

    return (
        raw,
        filtered,
        qc,
    )


# ============================================================
# Chunk output
# ============================================================

def write_chunk_stats(
    output,
    raw,
    filtered,
    region,
    info_threshold,
    qc,
):

    with open(
        output,
        "w",
    ) as out:

        out.write(
            "# Unified MAF-Stratified Accuracy Sufficient Statistics\n"
        )

        out.write(
            f"# Region: {region if region else 'whole file'}\n"
        )

        out.write(
            "# Site set: IMP-tagged imputed sites only\n"
        )

        out.write(
            "# Common observation set: truth GT + imputed GT + DS + GP\n"
        )

        out.write(
            f"# Filtered definition: INFO>={info_threshold}\n"
        )

        out.write(
            "MAF_Bin\t"
            "Mode\t"
            "SNPs\t"
            "Observations\t"
            "Sum_DS\t"
            "Sum_GT\t"
            "Sum_DS_GT\t"
            "Sum_DS2\t"
            "Sum_GT2\t"
            "Conc_Numerator\t"
            "Conc_Denominator\t"
            "N00\tN01\tN02\t"
            "N10\tN11\tN12\t"
            "N20\tN21\tN22\n"
        )

        for (
            bin_idx,
            label,
        ) in enumerate(
            BIN_LABELS
        ):

            for (
                mode,
                acc,
            ) in (
                ("RAW", raw[bin_idx]),
                ("FILTERED", filtered[bin_idx]),
            ):

                matrix = (
                    acc.nrc_matrix
                )

                out.write(
                    f"{label}\t"
                    f"{mode}\t"
                    f"{acc.snp_count}\t"
                    f"{acc.n}\t"
                    f"{acc.sum_ds:.15g}\t"
                    f"{acc.sum_gt:.15g}\t"
                    f"{acc.sum_ds_gt:.15g}\t"
                    f"{acc.sum_ds2:.15g}\t"
                    f"{acc.sum_gt2:.15g}\t"
                    f"{acc.conc_numerator:.15g}\t"
                    f"{acc.conc_denominator:.15g}\t"
                    +
                    "\t".join(
                        str(x)
                        for x in matrix
                    )
                    +
                    "\n"
                )

    print(
        "\nQC:",
        flush=True,
    )

    for (
        key,
        value,
    ) in qc.items():

        print(
            f"  {key}: {value:,}",
            flush=True,
        )

    print(
        f"\nChunk statistics written to: {output}",
        flush=True,
    )


# ============================================================
# CLI
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Unified MAF-stratified R2, "
            "GP-concordance and GT-NRC analysis."
        )
    )

    parser.add_argument(
        "-t",
        "--true-vcf",
        required=True,
    )

    parser.add_argument(
        "-i",
        "--imputed-vcf",
        required=True,
    )

    parser.add_argument(
        "-r",
        "--region",
        default=None,
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
    )

    parser.add_argument(
        "--info-threshold",
        type=float,
        default=0.8,
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
    )

    args = parser.parse_args()

    (
        raw,
        filtered,
        qc,
    ) = run_analysis(
        true_vcf=args.true_vcf,
        imputed_vcf=args.imputed_vcf,
        region=args.region,
        info_threshold=args.info_threshold,
        verbose=not args.quiet,
    )

    write_chunk_stats(
        output=args.output,
        raw=raw,
        filtered=filtered,
        region=args.region,
        info_threshold=args.info_threshold,
        qc=qc,
    )


if __name__ == "__main__":
    main()
