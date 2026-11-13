#!/usr/bin/env python3
import gzip
import numpy as np
import sys

vcf = sys.argv[1]

max_sites = 10000

checked = 0
mismatch = 0
sites = 0

def gt_to_dosage(gt):
    if gt in ("0/0", "0|0"):
        return 0
    elif gt in ("0/1", "1/0", "0|1", "1|0"):
        return 1
    elif gt in ("1/1", "1|1"):
        return 2
    else:
        return None

with gzip.open(vcf, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        cols = line.rstrip().split("\t")

        fmt = cols[8].split(":")

        if "GT" not in fmt or "GP" not in fmt:
            continue

        gt_idx = fmt.index("GT")
        gp_idx = fmt.index("GP")

        for sample in cols[9:]:
            parts = sample.split(":")

            if gt_idx >= len(parts) or gp_idx >= len(parts):
                continue

            gt = gt_to_dosage(parts[gt_idx])

            if gt is None:
                continue

            try:
                gp = [float(x) for x in parts[gp_idx].split(",")]
            except:
                continue

            if len(gp) != 3:
                continue

            gp_call = int(np.argmax(gp))

            checked += 1

            if gt != gp_call:
                mismatch += 1

                if mismatch <= 10:
                    print(
                        "Mismatch:",
                        parts[gt_idx],
                        parts[gp_idx],
                        gt,
                        gp_call
                    )

        sites += 1

        if sites >= max_sites:
            break

print("Sites checked:", sites)
print("Genotypes checked:", checked)
print("Mismatches:", mismatch)

if checked > 0:
    print("Mismatch rate:", mismatch / checked)
