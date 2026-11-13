#!/usr/bin/env python3
import csv
from pathlib import Path

root = Path.cwd()

# Previously calculated Ti/Tv statistics from position files.
titv = {
    "chr01": (184727238, 110392063),
    "chr02": (144210263, 86254095),
    "chr03": (157152524, 94346051),
    "chr04": (145518613, 87318150),
    "chr05": (163480496, 98241771),
    "chr06": (116494034, 69755380),
    "chr07": (152521473, 91247144),
    "chr08": (149668748, 89757610),
    "chr09": (135796411, 81673215),
    "chr10": (115860044, 68683058),
    "chr11": (101450685, 60748741),
    "chr12": (147252718, 88548161),
}

with open(root / "summary_v4/per_chromosome_summary.tsv") as f:
    stats = {
        r["chromosome"]: r
        for r in csv.DictReader(f, delimiter="\t")
    }

assert set(stats) == set(titv), "Chromosome mismatch"

rows = []
for chrom in sorted(titv):
    d = stats[chrom]
    ti, tv = titv[chrom]
    n = int(d["candidate_snps"])

    assert ti + tv == n, f"SNP count mismatch: {chrom}"

    rows.append([
        chrom.replace("chr", "Chr"),
        n, ti, tv,
        f"{ti/tv:.4f}",
        f'{float(d["snps_per_mb"]):.1f}',
        f'{float(d["observed_heterozygosity"])*100:.2f}',
        f'{float(d["missingness"])*100:.2f}'
    ])

with open(root / "summary_v4/genome_summary.tsv") as f:
    genome = {
        r["metric"]: r["value"]
        for r in csv.DictReader(f, delimiter="\t")
    }

total_ti = sum(v[0] for v in titv.values())
total_tv = sum(v[1] for v in titv.values())
total_n = sum(int(r[1]) for r in rows)

assert total_n == int(genome["matched_candidate_snps"])
assert total_ti + total_tv == total_n

rows.append([
    "Genome-wide",
    total_n, total_ti, total_tv,
    f"{total_ti/total_tv:.4f}",
    f'{float(genome["overall_snps_per_mb"]):.1f}',
    f'{float(genome["observed_heterozygosity"])*100:.2f}',
    f'{float(genome["missingness"])*100:.2f}'
])

output = root / "summary_v4/supplementary_snp_statistics.tsv"

with output.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "CHR", "Candidate SNP positions",
        "Transitions", "Transversions", "Ts/Tv",
        "SNP density (positions/Mb)",
        "Observed heterozygosity (%)",
        "Genotype missingness (%)"
    ])
    writer.writerows(rows)

print("Saved:", output)
print("Chromosomes:", len(titv))
print("Total SNPs:", total_n)
print("Genome-wide Ts/Tv:", f"{total_ti/total_tv:.4f}")
print("All chromosome SNP counts matched.")
