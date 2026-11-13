#!/usr/bin/env python3

import csv
import glob
import math

from collections import defaultdict


# ============================================================
# Input files
# ============================================================

SUMMARY_FILES = sorted(
    glob.glob(
        "NRC_chunks/NRC_chunk_*.summary.tsv"
    )
)

MATRIX_FILES = sorted(
    glob.glob(
        "NRC_chunks/NRC_chunk_*.matrix.tsv"
    )
)


EXPECTED_CHUNKS = 8


# ============================================================
# Check files
# ============================================================

if len(SUMMARY_FILES) != EXPECTED_CHUNKS:

    raise RuntimeError(
        f"Expected {EXPECTED_CHUNKS} summary files, "
        f"found {len(SUMMARY_FILES)}"
    )


if len(MATRIX_FILES) != EXPECTED_CHUNKS:

    raise RuntimeError(
        f"Expected {EXPECTED_CHUNKS} matrix files, "
        f"found {len(MATRIX_FILES)}"
    )


# ============================================================
# Merge summary counts
# ============================================================

SUM_COLS = [
    "Raw_SNPs",
    "Raw_Observations",
    "Raw_Correct_NonRef",
    "Raw_NRC_Evaluable",

    "Filtered_SNPs",
    "Filtered_Observations",
    "Filtered_Correct_NonRef",
    "Filtered_NRC_Evaluable",
]


summary = defaultdict(
    lambda: defaultdict(int)
)


for filename in SUMMARY_FILES:

    print(
        "Reading",
        filename,
    )

    with open(filename) as f:

        rows = (
            line
            for line in f
            if not line.startswith("#")
        )

        reader = csv.DictReader(
            rows,
            delimiter="\t",
        )

        for row in reader:

            maf_bin = row[
                "MAF_Bin"
            ]

            for col in SUM_COLS:

                summary[
                    maf_bin
                ][
                    col
                ] += int(
                    row[col]
                )


# ============================================================
# Merge 3 x 3 matrices
# ============================================================

MATRIX_COLS = [
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


matrices = defaultdict(
    lambda: defaultdict(
        lambda: [0] * 9
    )
)


for filename in MATRIX_FILES:

    print(
        "Reading",
        filename,
    )

    with open(filename) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t",
        )

        for row in reader:

            maf_bin = row[
                "MAF_Bin"
            ]

            mode = row[
                "Mode"
            ]

            for (
                j,
                col,
            ) in enumerate(
                MATRIX_COLS
            ):

                matrices[
                    maf_bin
                ][
                    mode
                ][
                    j
                ] += int(
                    row[col]
                )


# ============================================================
# Helpers
# ============================================================

def bin_start(label):

    return float(
        label.split(
            "-",
            1,
        )[0]
    )


def calc_nrc(
    correct,
    evaluable,
):

    if evaluable > 0:

        return (
            correct
            /
            evaluable
        )

    return math.nan


def fmt(x):

    if math.isnan(x):
        return "NA"

    return f"{x:.10f}"


# ============================================================
# Output filenames
# ============================================================

FINAL_SUMMARY = (
    "results_MAF-NRC-all.txt"
)

FINAL_MATRIX = (
    "results_MAF-NRC-all-matrices.txt"
)


# ============================================================
# Write final summary
# ============================================================

with open(
    FINAL_SUMMARY,
    "w",
) as out:

    out.write(
        "# MAF-Stratified NRC Results\n"
    )

    out.write(
        "# Combined across 8 genomic chunks\n"
    )

    out.write(
        "# Hard call source: imputed GT field\n"
    )

    out.write(
        "# Site set: IMP-tagged imputed sites only\n"
    )

    out.write(
        "# Filtered definition: INFO>=0.8\n"
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


    for maf_bin in sorted(
        summary,
        key=bin_start,
    ):

        x = summary[
            maf_bin
        ]


        raw_nrc = calc_nrc(
            x[
                "Raw_Correct_NonRef"
            ],
            x[
                "Raw_NRC_Evaluable"
            ],
        )


        filt_nrc = calc_nrc(
            x[
                "Filtered_Correct_NonRef"
            ],
            x[
                "Filtered_NRC_Evaluable"
            ],
        )


        # ----------------------------------------------------
        # Matrix-based QC
        # ----------------------------------------------------

        raw_matrix = matrices[
            maf_bin
        ][
            "RAW"
        ]

        filt_matrix = matrices[
            maf_bin
        ][
            "FILTERED"
        ]


        raw_total = sum(
            raw_matrix
        )

        filt_total = sum(
            filt_matrix
        )


        raw_correct = (
            raw_matrix[4]
            +
            raw_matrix[8]
        )

        filt_correct = (
            filt_matrix[4]
            +
            filt_matrix[8]
        )


        raw_eval = (
            raw_total
            -
            raw_matrix[0]
        )

        filt_eval = (
            filt_total
            -
            filt_matrix[0]
        )


        # ----------------------------------------------------
        # Consistency checks
        # ----------------------------------------------------

        if (
            raw_total
            !=
            x["Raw_Observations"]
        ):

            raise RuntimeError(
                f"Raw observation mismatch in {maf_bin}"
            )


        if (
            filt_total
            !=
            x["Filtered_Observations"]
        ):

            raise RuntimeError(
                f"Filtered observation mismatch in {maf_bin}"
            )


        if (
            raw_correct
            !=
            x["Raw_Correct_NonRef"]
            or
            raw_eval
            !=
            x["Raw_NRC_Evaluable"]
        ):

            raise RuntimeError(
                f"Raw matrix mismatch in {maf_bin}"
            )


        if (
            filt_correct
            !=
            x["Filtered_Correct_NonRef"]
            or
            filt_eval
            !=
            x["Filtered_NRC_Evaluable"]
        ):

            raise RuntimeError(
                f"Filtered matrix mismatch in {maf_bin}"
            )


        # ----------------------------------------------------
        # Write
        # ----------------------------------------------------

        out.write(

            f"{maf_bin}\t"

            f"{x['Raw_SNPs']}\t"

            f"{x['Raw_Observations']}\t"

            f"{x['Raw_Correct_NonRef']}\t"

            f"{x['Raw_NRC_Evaluable']}\t"

            f"{fmt(raw_nrc)}\t"

            f"{x['Filtered_SNPs']}\t"

            f"{x['Filtered_Observations']}\t"

            f"{x['Filtered_Correct_NonRef']}\t"

            f"{x['Filtered_NRC_Evaluable']}\t"

            f"{fmt(filt_nrc)}\n"

        )


# ============================================================
# Write merged matrices
# ============================================================

with open(
    FINAL_MATRIX,
    "w",
) as out:

    out.write(
        "MAF_Bin\t"
        "Mode\t"
        "N00\tN01\tN02\t"
        "N10\tN11\tN12\t"
        "N20\tN21\tN22\n"
    )


    for maf_bin in sorted(
        matrices,
        key=bin_start,
    ):

        for mode in (
            "RAW",
            "FILTERED",
        ):

            out.write(

                maf_bin
                +
                "\t"
                +
                mode
                +
                "\t"
                +
                "\t".join(
                    str(x)
                    for x in matrices[
                        maf_bin
                    ][
                        mode
                    ]
                )
                +
                "\n"

            )


print()
print(
    "Merged successfully."
)

print(
    "Summary:",
    FINAL_SUMMARY,
)

print(
    "Matrices:",
    FINAL_MATRIX,
)
