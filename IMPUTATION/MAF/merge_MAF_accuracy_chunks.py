#!/usr/bin/env python3

import csv
import glob
import math
from collections import defaultdict


INPUT_PATTERN = (
    "MAF_accuracy_chunks/"
    "accuracy_chunk_*.stats.tsv"
)

EXPECTED_CHUNKS = 8


FILES = sorted(
    glob.glob(
        INPUT_PATTERN
    )
)


if len(FILES) != EXPECTED_CHUNKS:

    raise RuntimeError(
        f"Expected {EXPECTED_CHUNKS} chunk files, "
        f"but found {len(FILES)}"
    )


# ============================================================
# Columns
# ============================================================

INTEGER_COLUMNS = [
    "SNPs",
    "Observations",

    "N00",
    "N01",
    "N02",

    "N10",
    "N11",
    "N12",

    "N20",
    "N21",
    "N22",
]


FLOAT_COLUMNS = [
    "Sum_DS",
    "Sum_GT",
    "Sum_DS_GT",
    "Sum_DS2",
    "Sum_GT2",

    "Conc_Numerator",
    "Conc_Denominator",
]


# ============================================================
# Combined statistics
#
# key = (MAF_Bin, Mode)
# ============================================================

combined = defaultdict(
    lambda: {
        **{
            col: 0
            for col in INTEGER_COLUMNS
        },
        **{
            col: 0.0
            for col in FLOAT_COLUMNS
        },
    }
)


# ============================================================
# Read all chunks
# ============================================================

for filename in FILES:

    print(
        "Reading:",
        filename,
    )

    with open(filename) as handle:

        rows = (
            line
            for line in handle
            if not line.startswith("#")
        )

        reader = csv.DictReader(
            rows,
            delimiter="\t",
        )

        for row in reader:

            key = (
                row["MAF_Bin"],
                row["Mode"],
            )

            acc = combined[
                key
            ]

            for col in INTEGER_COLUMNS:

                acc[col] += int(
                    row[col]
                )

            for col in FLOAT_COLUMNS:

                acc[col] += float(
                    row[col]
                )


# ============================================================
# Metric functions
# ============================================================

def calculate_r2(x):

    n = x["Observations"]

    if n < 3:
        return math.nan

    numerator = (
        n * x["Sum_DS_GT"]
        -
        x["Sum_DS"]
        *
        x["Sum_GT"]
    )

    denominator_ds = (
        n * x["Sum_DS2"]
        -
        x["Sum_DS"] ** 2
    )

    denominator_gt = (
        n * x["Sum_GT2"]
        -
        x["Sum_GT"] ** 2
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


def calculate_concordance(x):

    denominator = (
        x["Conc_Denominator"]
    )

    if denominator <= 0:

        return math.nan

    return (
        x["Conc_Numerator"]
        /
        denominator
    )


def calculate_nrc(x):

    total = (
        x["N00"]
        +
        x["N01"]
        +
        x["N02"]
        +
        x["N10"]
        +
        x["N11"]
        +
        x["N12"]
        +
        x["N20"]
        +
        x["N21"]
        +
        x["N22"]
    )

    evaluable = (
        total
        -
        x["N00"]
    )

    correct_nonref = (
        x["N11"]
        +
        x["N22"]
    )

    if evaluable <= 0:

        return math.nan

    return (
        correct_nonref
        /
        evaluable
    )


def format_value(value):

    if math.isnan(value):
        return "NA"

    return f"{value:.10f}"


def bin_start(label):

    return float(
        label.split(
            "-",
            1,
        )[0]
    )


# ============================================================
# Obtain MAF bins
# ============================================================

MAF_BINS = sorted(
    set(
        key[0]
        for key in combined
    ),
    key=bin_start,
)


# ============================================================
# Final summary output
# ============================================================

SUMMARY_OUTPUT = (
    "results_MAF-accuracy-unified.txt"
)


with open(
    SUMMARY_OUTPUT,
    "w",
) as out:

    out.write(
        "# Unified MAF-Stratified Imputation Accuracy\n"
    )

    out.write(
        "# Combined across 8 genomic chunks\n"
    )

    out.write(
        "# Site set: IMP-tagged imputed sites only\n"
    )

    out.write(
        "# Common observation set: truth GT + imputed GT + DS + GP\n"
    )

    out.write(
        "# Filtered definition: INFO>=0.8\n"
    )

    out.write(
        "# R2 = pooled Pearson correlation squared between truth GT and imputed DS\n"
    )

    out.write(
        "# Concordance = pooled GP probability assigned to the true genotype / total GP probability\n"
    )

    out.write(
        "# NRC = (N11+N22)/(Ntotal-N00)\n"
    )

    out.write(
        "MAF_Bin\t"
        "Raw_SNPs\t"
        "Raw_Observations\t"
        "Raw_R2\t"
        "Raw_Concordance\t"
        "Raw_NRC\t"
        "Filtered_SNPs\t"
        "Filtered_Observations\t"
        "Filtered_R2\t"
        "Filtered_Concordance\t"
        "Filtered_NRC\n"
    )


    for maf_bin in MAF_BINS:

        raw = combined[
            (
                maf_bin,
                "RAW",
            )
        ]

        filt = combined[
            (
                maf_bin,
                "FILTERED",
            )
        ]


        raw_r2 = calculate_r2(
            raw
        )

        raw_conc = calculate_concordance(
            raw
        )

        raw_nrc = calculate_nrc(
            raw
        )


        filt_r2 = calculate_r2(
            filt
        )

        filt_conc = calculate_concordance(
            filt
        )

        filt_nrc = calculate_nrc(
            filt
        )


        out.write(
            f"{maf_bin}\t"
            f"{raw['SNPs']}\t"
            f"{raw['Observations']}\t"
            f"{format_value(raw_r2)}\t"
            f"{format_value(raw_conc)}\t"
            f"{format_value(raw_nrc)}\t"
            f"{filt['SNPs']}\t"
            f"{filt['Observations']}\t"
            f"{format_value(filt_r2)}\t"
            f"{format_value(filt_conc)}\t"
            f"{format_value(filt_nrc)}\n"
        )


# ============================================================
# Detailed sufficient-statistics output
# ============================================================

DETAIL_OUTPUT = (
    "results_MAF-accuracy-unified-sufficient-stats.txt"
)


with open(
    DETAIL_OUTPUT,
    "w",
) as out:

    columns = (
        [
            "MAF_Bin",
            "Mode",
        ]
        +
        INTEGER_COLUMNS
        +
        FLOAT_COLUMNS
    )

    out.write(
        "\t".join(columns)
        +
        "\n"
    )


    for maf_bin in MAF_BINS:

        for mode in (
            "RAW",
            "FILTERED",
        ):

            x = combined[
                (
                    maf_bin,
                    mode,
                )
            ]

            values = [
                maf_bin,
                mode,
            ]

            for col in INTEGER_COLUMNS:

                values.append(
                    str(
                        x[col]
                    )
                )

            for col in FLOAT_COLUMNS:

                values.append(
                    f"{x[col]:.15g}"
                )

            out.write(
                "\t".join(values)
                +
                "\n"
            )


print()
print(
    "Merge completed successfully."
)

print(
    "Final summary:",
    SUMMARY_OUTPUT,
)

print(
    "Sufficient statistics:",
    DETAIL_OUTPUT,
)
