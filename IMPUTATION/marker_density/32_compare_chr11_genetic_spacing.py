#!/usr/bin/env python3

import numpy as np

MAPFILE = "/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map/LD.GenMap.chr11.txt"

datasets = {
    "Actual_array":
        "Chr11/array_real/Chr11.ARRAY50.sites.tsv",

    "Random_rep1":
        "Chr11/random1493/rep1/Chr11.random1493.rep1.sites.tsv",

    "Random_rep2":
        "Chr11/random1493/rep2/Chr11.random1493.rep2.sites.tsv",

    "Random_rep3":
        "Chr11/random1493/rep3/Chr11.random1493.rep3.sites.tsv",
}

GENETIC_WINDOW = 5.0   # cM


def read_map(path):
    pos = []
    cm = []

    with open(path) as f:
        header = next(f)

        for line in f:
            if not line.strip():
                continue

            x = line.split()
            pos.append(int(x[0]))
            cm.append(float(x[2]))

    pos = np.asarray(pos, dtype=np.int64)
    cm = np.asarray(cm, dtype=float)

    # Safety checks
    if np.any(np.diff(pos) < 0):
        raise RuntimeError("Genetic map positions are not sorted.")

    if np.any(np.diff(cm) < -1e-10):
        raise RuntimeError("Genetic map cM values are not monotonic.")

    return pos, cm


def read_markers(path):
    pos = []

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            x = line.split()
            pos.append(int(x[1]))

    return np.asarray(sorted(set(pos)), dtype=np.int64)


map_pos, map_cm = read_map(MAPFILE)

MAP_MIN_BP = map_pos[0]
MAP_MAX_BP = map_pos[-1]
MAP_MIN_CM = map_cm[0]
MAP_MAX_CM = map_cm[-1]

print(
    "\t".join([
        "Dataset",
        "N_markers",
        "Markers_inside_map",
        "Markers_outside_map",
        "Median_gap_cM",
        "Q75_gap_cM",
        "Q90_gap_cM",
        "Q95_gap_cM",
        "Max_gap_cM",
        "Mean_gap_cM",
        "Window5cM_mean",
        "Window5cM_SD",
        "Window5cM_CV",
        "Empty_5cM_windows",
        "Empty_5cM_fraction"
    ])
)

for name, path in datasets.items():

    pos = read_markers(path)

    inside = (pos >= MAP_MIN_BP) & (pos <= MAP_MAX_BP)

    n_inside = int(np.sum(inside))
    n_outside = int(np.sum(~inside))

    # For markers outside the map boundaries, np.interp assigns
    # the nearest boundary cM. We report their number explicitly.
    marker_cm = np.interp(
        pos,
        map_pos,
        map_cm
    )

    marker_cm = np.sort(marker_cm)

    gaps = np.diff(marker_cm)

    # 5-cM windows across the complete genetic map
    bins = np.arange(
        0,
        MAP_MAX_CM + GENETIC_WINDOW,
        GENETIC_WINDOW
    )

    if bins[-1] < MAP_MAX_CM:
        bins = np.append(bins, MAP_MAX_CM)

    counts, _ = np.histogram(marker_cm, bins=bins)

    mean_count = counts.mean()
    sd_count = counts.std(ddof=1)

    cv = (
        sd_count / mean_count
        if mean_count > 0
        else np.nan
    )

    empty = int(np.sum(counts == 0))

    values = [
        name,
        len(pos),
        n_inside,
        n_outside,
        np.median(gaps),
        np.percentile(gaps, 75),
        np.percentile(gaps, 90),
        np.percentile(gaps, 95),
        np.max(gaps),
        np.mean(gaps),
        mean_count,
        sd_count,
        cv,
        empty,
        empty / len(counts)
    ]

    formatted = [
        values[0],
        str(values[1]),
        str(values[2]),
        str(values[3]),
        f"{values[4]:.6f}",
        f"{values[5]:.6f}",
        f"{values[6]:.6f}",
        f"{values[7]:.6f}",
        f"{values[8]:.6f}",
        f"{values[9]:.6f}",
        f"{values[10]:.4f}",
        f"{values[11]:.4f}",
        f"{values[12]:.4f}",
        str(values[13]),
        f"{values[14]:.4f}"
    ]

    print("\t".join(formatted))
