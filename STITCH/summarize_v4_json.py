#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path.cwd()
MANIFEST = ROOT / "config/manifest.tsv"
LENGTHS = ROOT / "config/chromosome_lengths.tsv"
RESULTS = ROOT / "per_file_v4"
OUTPUT = ROOT / "summary_v4"
OUTPUT.mkdir(exist_ok=True)

EXPECTED_SNPS = 2741098686

COUNT_KEYS = [
    "bcf_records_at_candidate_positions",
    "matched_records",
    "sample_count",
    "missing_gt",
    "called_diploid_gt",
    "hom_gt",
    "het_gt",
    "sites_with_called_alleles",
    "unmatched_bcf_records",
]

def write_tsv(path, header, rows):
    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)

def rate(a, b):
    return a / b if b else None

def fmt(x):
    return f"{x:.8f}" if x is not None else "NA"

# Read chromosome lengths.
chrom_lengths = {}
with LENGTHS.open() as f:
    for line in f:
        fields = line.strip().split()
        if len(fields) < 2:
            continue
        match = re.fullmatch(r"(?:chr)?0?([1-9]|1[0-2])",
                             fields[0], re.I)
        if match and fields[1].isdigit():
            chrom = f"chr{int(match.group(1)):02d}"
            chrom_lengths[chrom] = int(fields[1])

if len(chrom_lengths) != 12:
    raise SystemExit(
        f"ERROR: expected 12 chromosome lengths, "
        f"found {len(chrom_lengths)}: {chrom_lengths}"
    )

# Aggregate data.
total = Counter()
by_chrom = defaultdict(Counter)
af_total = Counter()
maf_total = Counter()
problems = []
chunk_rows = []
sample_counts_seen = set()

with MANIFEST.open(newline="") as f:
    manifest = list(csv.DictReader(f, delimiter="\t"))

if len(manifest) != 165:
    raise SystemExit(
        f"ERROR: expected 165 manifest rows, got {len(manifest)}"
    )

for row in manifest:
    chunk = row["chunk"].strip()
    chrom_match = re.match(r"PA_chr(\d{2})_", chunk)

    if not chrom_match:
        raise SystemExit(f"Invalid chunk name: {chunk}")

    chrom = "chr" + chrom_match.group(1)
    path = RESULTS / f"{chunk}.json"

    if not path.is_file():
        raise SystemExit(f"Missing JSON: {path}")

    with path.open() as f:
        data = json.load(f)

    if data.get("chunk") != chunk:
        problems.append(f"{chunk}: JSON chunk name mismatch")

    counts = data["counts"]

    for key in COUNT_KEYS:
        if key not in counts:
            raise SystemExit(f"{chunk}: missing field {key}")

    n = int(counts["matched_records"])
    samples = int(counts["sample_count"])
    sample_counts_seen.add(samples)

    called = int(counts["called_diploid_gt"])
    missing = int(counts["missing_gt"])
    hom = int(counts["hom_gt"])
    het = int(counts["het_gt"])
    expected_gt = n * samples

    # Fundamental genotype accounting.
    if called != hom + het:
        problems.append(
            f"{chunk}: called != hom + het"
        )

    if called + missing != expected_gt:
        problems.append(
            f"{chunk}: called + missing "
            f"!= matched_records * sample_count; "
            f"difference={expected_gt-called-missing}"
        )

    # Matching and AF/MAF checks.
    sites_called = int(counts["sites_with_called_alleles"])
    af = {int(k): int(v)
          for k, v in data["af_bins"].items()}
    maf = {int(k): int(v)
           for k, v in data["maf_bins"].items()}

    if set(af) != set(range(20)):
        problems.append(f"{chunk}: unexpected AF bins")
    if set(maf) != set(range(20)):
        problems.append(f"{chunk}: unexpected MAF bins")

    if sum(af.values()) != sites_called:
        problems.append(
            f"{chunk}: AF bin sum != sites_with_called_alleles"
        )
    if sum(maf.values()) != sites_called:
        problems.append(
            f"{chunk}: MAF bin sum != sites_with_called_alleles"
        )

    for key, value in counts.items():
        if key == "sample_count":
            continue
        if isinstance(value, int):
            total[key] += value
            by_chrom[chrom][key] += value

    total["expected_gt"] += expected_gt
    by_chrom[chrom]["expected_gt"] += expected_gt

    af_total.update(af)
    maf_total.update(maf)

    chunk_rows.append([
        chunk, chrom, n, samples, called, hom, het, missing,
        fmt(rate(het, called)),
        fmt(rate(missing, expected_gt)),
        sites_called,
        counts["unmatched_bcf_records"]
    ])

# Save individual chunk statistics.
write_tsv(
    OUTPUT / "per_chunk_summary.tsv",
    [
        "chunk", "chromosome", "matched_snps", "samples",
        "called_gt", "hom_gt", "het_gt", "missing_gt",
        "observed_heterozygosity",
        "missingness",
        "sites_with_called_alleles",
        "unmatched_bcf_records"
    ],
    chunk_rows
)

# Save chromosome statistics and SNP densities.
chrom_rows = []

for chrom in sorted(chrom_lengths):
    c = by_chrom[chrom]
    length_bp = chrom_lengths[chrom]
    n = c["matched_records"]

    chrom_rows.append([
        chrom,
        length_bp,
        n,
        fmt(n / (length_bp / 1_000_000)),
        c["called_diploid_gt"],
        c["het_gt"],
        c["missing_gt"],
        fmt(rate(c["het_gt"], c["called_diploid_gt"])),
        fmt(rate(c["missing_gt"], c["expected_gt"]))
    ])

write_tsv(
    OUTPUT / "per_chromosome_summary.tsv",
    [
        "chromosome", "length_bp", "candidate_snps",
        "snps_per_mb", "called_gt", "het_gt", "missing_gt",
        "observed_heterozygosity", "missingness"
    ],
    chrom_rows
)

# Save AF and MAF distributions.
write_tsv(
    OUTPUT / "allele_frequency_bins.tsv",
    ["bin", "lower", "upper", "af_count", "maf_count"],
    [
        [
            i, f"{i*0.05:.2f}", f"{(i+1)*0.05:.2f}",
            af_total[i], maf_total[i]
        ]
        for i in range(20)
    ]
)

# Global summary.
global_rows = [
    ["chunks", len(manifest)],
    ["sample_counts_observed",
     ",".join(map(str, sorted(sample_counts_seen)))],
    ["expected_candidate_snps", EXPECTED_SNPS],
    ["matched_candidate_snps", total["matched_records"]],
    ["bcf_records_at_candidate_positions",
     total["bcf_records_at_candidate_positions"]],
    ["unmatched_bcf_records", total["unmatched_bcf_records"]],
    ["sites_with_called_alleles",
     total["sites_with_called_alleles"]],
    ["expected_gt", total["expected_gt"]],
    ["called_diploid_gt", total["called_diploid_gt"]],
    ["hom_gt", total["hom_gt"]],
    ["het_gt", total["het_gt"]],
    ["missing_gt", total["missing_gt"]],
    ["unaccounted_gt",
     total["expected_gt"] -
     total["called_diploid_gt"] -
     total["missing_gt"]],
    ["observed_heterozygosity",
     fmt(rate(total["het_gt"], total["called_diploid_gt"]))],
    ["missingness",
     fmt(rate(total["missing_gt"], total["expected_gt"]))],
    ["af_bin_sum", sum(af_total.values())],
    ["maf_bin_sum", sum(maf_total.values())],
    ["total_chromosome_length_bp", sum(chrom_lengths.values())],
    ["overall_snps_per_mb",
     fmt(total["matched_records"] /
         (sum(chrom_lengths.values()) / 1_000_000))],
]

write_tsv(
    OUTPUT / "genome_summary.tsv",
    ["metric", "value"],
    global_rows
)

with (OUTPUT / "validation_issues.txt").open("w") as f:
    if problems:
        f.write("\n".join(problems) + "\n")
    else:
        f.write("No internal accounting issues detected.\n")

print("\n========== GENOME SUMMARY ==========")
for key, value in global_rows:
    print(f"{key}: {value}")

print("\n========== VALIDATION ==========")
print("Candidate SNP total matches expectation:",
      total["matched_records"] == EXPECTED_SNPS)
print("Sample counts:", sorted(sample_counts_seen))
print("Internal accounting issues:", len(problems))
print("Output directory:", OUTPUT)
