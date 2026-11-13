#!/usr/bin/env python3

import csv
import math
import statistics

files = {
    "Actual_array":
        "Chr11/maf_comparison/Actual_array.maf.tsv",
    "Random_rep1":
        "Chr11/maf_comparison/Random_rep1.maf.tsv",
    "Random_rep2":
        "Chr11/maf_comparison/Random_rep2.maf.tsv",
    "Random_rep3":
        "Chr11/maf_comparison/Random_rep3.maf.tsv",
}

header = [
    "Dataset",
    "N_markers",
    "N_valid_MAF",
    "MAF_eq_0",
    "MAF_eq_0_fraction",
    "MAF_0_0.01",
    "MAF_0_0.01_fraction",
    "MAF_0.01_0.05",
    "MAF_0.01_0.05_fraction",
    "MAF_0.05_0.10",
    "MAF_0.05_0.10_fraction",
    "MAF_0.10_0.20",
    "MAF_0.10_0.20_fraction",
    "MAF_ge_0.20",
    "MAF_ge_0.20_fraction",
    "Mean_MAF",
    "Median_MAF",
]

print("\t".join(header))

for name, path in files.items():

    mafs = []

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            if row["MAF"] == "NA":
                continue
            mafs.append(float(row["MAF"]))

    n = len(mafs)

    n0 = sum(x == 0 for x in mafs)
    n001 = sum(0 < x < 0.01 for x in mafs)
    n005 = sum(0.01 <= x < 0.05 for x in mafs)
    n010 = sum(0.05 <= x < 0.10 for x in mafs)
    n020 = sum(0.10 <= x < 0.20 for x in mafs)
    n20 = sum(x >= 0.20 for x in mafs)

    vals = [
        name,
        1493,
        n,

        n0,
        n0 / n if n else math.nan,

        n001,
        n001 / n if n else math.nan,

        n005,
        n005 / n if n else math.nan,

        n010,
        n010 / n if n else math.nan,

        n020,
        n020 / n if n else math.nan,

        n20,
        n20 / n if n else math.nan,

        statistics.mean(mafs) if mafs else math.nan,
        statistics.median(mafs) if mafs else math.nan,
    ]

    out = []
    for i, x in enumerate(vals):
        if isinstance(x, float):
            out.append(f"{x:.6f}")
        else:
            out.append(str(x))

    print("\t".join(out))
