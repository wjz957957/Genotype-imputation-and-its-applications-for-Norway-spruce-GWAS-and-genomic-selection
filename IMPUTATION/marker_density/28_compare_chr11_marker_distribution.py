#!/usr/bin/env python3

import numpy as np

CHR_LENGTH = 954_579_047
WINDOW = 10_000_000

datasets = {
    "Original_array_2889":
        "Chr11/array_real/Chr11.original_array.sites.tsv",

    "Actual_usable_array_1493":
        "Chr11/array_real/Chr11.ARRAY50.sites.tsv",

    "Random1493_rep1":
        "Chr11/random1493/rep1/Chr11.random1493.rep1.sites.tsv",

    "Random1493_rep2":
        "Chr11/random1493/rep2/Chr11.random1493.rep2.sites.tsv",

    "Random1493_rep3":
        "Chr11/random1493/rep3/Chr11.random1493.rep3.sites.tsv",
}


def read_positions(path):
    pos = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            x = line.split()
            pos.append(int(x[1]))
    return np.array(sorted(set(pos)), dtype=np.int64)


header = [
    "Dataset",
    "N_markers",
    "SNP_per_Mb",
    "Mean_gap_kb",
    "Median_gap_kb",
    "Q75_gap_kb",
    "Q90_gap_kb",
    "Q95_gap_kb",
    "Max_gap_Mb",
    "Window10Mb_mean",
    "Window10Mb_SD",
    "Window10Mb_CV",
    "Empty_10Mb_windows",
    "Empty_10Mb_fraction",
]

print("\t".join(header))

for name, path in datasets.items():

    pos = read_positions(path)
    gaps = np.diff(pos)

    # Include the final partial chromosome window
    bins = np.arange(0, CHR_LENGTH + WINDOW, WINDOW)

    if bins[-1] < CHR_LENGTH:
        bins = np.append(bins, CHR_LENGTH)

    counts, _ = np.histogram(pos, bins=bins)

    mean_count = counts.mean()
    sd_count = counts.std(ddof=1)
    cv = sd_count / mean_count if mean_count > 0 else np.nan

    empty = np.sum(counts == 0)

    result = [
        name,
        len(pos),
        len(pos) / (CHR_LENGTH / 1e6),
        np.mean(gaps) / 1e3,
        np.median(gaps) / 1e3,
        np.percentile(gaps, 75) / 1e3,
        np.percentile(gaps, 90) / 1e3,
        np.percentile(gaps, 95) / 1e3,
        np.max(gaps) / 1e6,
        mean_count,
        sd_count,
        cv,
        empty,
        empty / len(counts),
    ]

    formatted = [
        result[0],
        str(result[1]),
        f"{result[2]:.4f}",
        f"{result[3]:.3f}",
        f"{result[4]:.3f}",
        f"{result[5]:.3f}",
        f"{result[6]:.3f}",
        f"{result[7]:.3f}",
        f"{result[8]:.3f}",
        f"{result[9]:.3f}",
        f"{result[10]:.3f}",
        f"{result[11]:.4f}",
        str(result[12]),
        f"{result[13]:.4f}",
    ]

    print("\t".join(formatted))
