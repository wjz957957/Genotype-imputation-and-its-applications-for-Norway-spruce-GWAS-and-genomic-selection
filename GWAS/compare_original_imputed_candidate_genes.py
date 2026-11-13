#!/usr/bin/env python3

import csv
import gzip
import argparse
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--root", required=True)
parser.add_argument("--work", required=True)
parser.add_argument("--gff", required=True)
parser.add_argument("--window", type=int, default=50000)
parser.add_argument("--nearest-window", type=int, default=300000)
args = parser.parse_args()

ROOT = Path(args.root)
WORK = Path(args.work)
WINDOW = args.window

FILES = {
    ("BB", "Original"): ROOT / "GWAS/GWAS-Bud/gwas_clumps_revised_origin_Bud.clumped",
    ("BB", "Imputed"): ROOT / "GWAS/GWAS-Bud/gwas_clumps_revised_imputed_Bud.clumped",
    ("DBH", "Original"): ROOT / "GWAS/GWAS-DBH/gwas_clumps_revised_origin_DBH.clumped",
    ("DBH", "Imputed"): ROOT / "GWAS/GWAS-DBH/gwas_clumps_revised_imputed_DBH.clumped",
    ("FD", "Original"): ROOT / "GWAS/GWAS-FD/gwas_clumps_revised_origin_FD.clumped",
    ("FD", "Imputed"): ROOT / "GWAS/GWAS-FD/gwas_clumps_revised_imputed_FD.clumped",
}

EXPECTED = {
    ("BB", "Original"): 6,
    ("BB", "Imputed"): 23,
    ("DBH", "Original"): 1,
    ("DBH", "Imputed"): 3,
    ("FD", "Original"): 7,
    ("FD", "Imputed"): 13,
}

def chromosome_name(value):
    value = value.replace("PA_chr", "").replace("chr", "")
    return f"PA_chr{int(value):02d}"

def parse_attributes(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            result[k] = v
    return result

def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=columns,
            delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

# --------------------------------------------------
# 1. Read six existing PLINK clumping results
# --------------------------------------------------

regions = defaultdict(list)

for (trait, dataset), path in FILES.items():

    if not path.exists():
        raise FileNotFoundError(path)

    with open(path) as f:
        for line in f:

            fields = line.split()

            if not fields or fields[0] == "CHR":
                continue

            if len(fields) < 5:
                continue

            chrom = chromosome_name(fields[0])
            snp = fields[2]
            bp = int(fields[3])
            p = fields[4]

            region_id = f"{trait}_{dataset}_{chrom}_{bp}"

            regions[(trait, dataset)].append({
                "Region_ID": region_id,
                "Trait": trait,
                "Dataset": dataset,
                "Chromosome": chrom,
                "Index_SNP": snp,
                "Position": bp,
                "P": p
            })

            if len(regions[(trait, dataset)]) > EXPECTED[(trait, dataset)]:
                raise ValueError(
                    f"Unexpected region count: {trait} {dataset}"
                )

for key, expected in EXPECTED.items():
    observed = len(regions[key])
    print(f"{key[0]} {key[1]}: {observed} regions")
    if observed != expected:
        raise ValueError(
            f"Expected {expected}, observed {observed}: {key}"
        )

region_columns = [
    "Region_ID", "Trait", "Dataset",
    "Chromosome", "Index_SNP", "Position", "P"
]

for dataset in ["Original", "Imputed"]:

    rows = [
        r for key in regions if key[1] == dataset
        for r in regions[key]
    ]

    folder = WORK / "01_region_definition" / dataset.lower()

    write_tsv(
        folder / f"{dataset.lower()}_clumped_regions.tsv",
        rows, region_columns
    )

# --------------------------------------------------
# 2. Read Norway spruce v2 GFF3 gene features
# --------------------------------------------------

print("\nLoading Norway spruce v2 GFF3...")

genes_by_chr = defaultdict(list)
gene_info = {}

with gzip.open(args.gff, "rt") as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9 or fields[2] != "gene":
            continue

        chrom, source = fields[0], fields[1]
        start, end = int(fields[3]), int(fields[4])
        attrs = parse_attributes(fields[8])

        gid = attrs.get("ID")

        if not gid:
            continue

        gene = {
            "Gene_ID": gid,
            "Chromosome": chrom,
            "Gene_start": start,
            "Gene_end": end,
            "Gene_source": source,
            "Gene_strand": fields[6],
            "Gene_name": attrs.get("Name", ""),
            "Gene_biotype": attrs.get(
                "biotype", attrs.get("gene_biotype", "")
            ),
            "Gene_description": attrs.get(
                "description", attrs.get("Note", "")
            )
        }

        genes_by_chr[chrom].append(gene)
        gene_info[gid] = gene

for chrom in genes_by_chr:
    genes_by_chr[chrom].sort(
        key=lambda x: x["Gene_start"]
    )

print("Total gene features:", len(gene_info))

# --------------------------------------------------
# 3. Map significant regions to genes
#    Main candidate definition: +/- 50 kb
# --------------------------------------------------

mapping = defaultdict(list)
nearest_rows = []

for key, region_list in regions.items():

    trait, dataset = key

    for region in region_list:

        chrom = region["Chromosome"]
        bp = region["Position"]

        genes = genes_by_chr.get(chrom, [])

        left = max(1, bp - WINDOW)
        right = bp + WINDOW

        matched = []

        for gene in genes:

            if gene["Gene_start"] > right:
                break

            if gene["Gene_end"] < left:
                continue

            if gene["Gene_start"] <= bp <= gene["Gene_end"]:
                distance = 0
                relationship = "LEAD_WITHIN_GENE"

            elif gene["Gene_end"] < bp:
                distance = bp - gene["Gene_end"]
                relationship = "GENE_LEFT_OF_LEAD"

            else:
                distance = gene["Gene_start"] - bp
                relationship = "GENE_RIGHT_OF_LEAD"

            row = {
                **region,
                **gene,
                "Distance_to_lead_bp": distance,
                "Lead_gene_relationship": relationship,
                "Search_window_bp": WINDOW
            }

            matched.append(row)

        if matched:
            mapping[key].extend(matched)

        else:
            mapping[key].append({
                **region,
                "Gene_ID": "",
                "Distance_to_lead_bp": "",
                "Lead_gene_relationship": "NO_GENE_WITHIN_WINDOW",
                "Search_window_bp": WINDOW
            })

        # Record nearest gene separately
        if genes:

            nearest = min(
                genes,
                key=lambda g: (
                    0 if g["Gene_start"] <= bp <= g["Gene_end"]
                    else min(abs(bp-g["Gene_start"]),
                             abs(bp-g["Gene_end"]))
                )
            )

            dist = (
                0 if nearest["Gene_start"] <= bp <= nearest["Gene_end"]
                else min(
                    abs(bp-nearest["Gene_start"]),
                    abs(bp-nearest["Gene_end"])
                )
            )

            nearest_rows.append({
                **region,
                "Nearest_gene_ID": nearest["Gene_ID"],
                "Nearest_gene_distance_bp": dist,
                "Within_nearest_window":
                    "YES" if dist <= args.nearest_window else "NO"
            })

        else:
            nearest_rows.append({
                **region,
                "Nearest_gene_ID": "",
                "Nearest_gene_distance_bp": "",
                "Within_nearest_window": "NO"
            })

# Save separate mapping results

map_columns = [
    "Region_ID", "Trait", "Dataset",
    "Index_SNP", "Chromosome", "Position", "P",
    "Gene_ID", "Gene_start", "Gene_end",
    "Gene_source", "Gene_strand",
    "Gene_name", "Gene_biotype",
    "Distance_to_lead_bp",
    "Lead_gene_relationship",
    "Search_window_bp"
]

for dataset in ["Original", "Imputed"]:

    rows = [
        r for key in mapping if key[1] == dataset
        for r in mapping[key]
    ]

    folder = (
        WORK / "02_gene_annotation" / dataset.lower()
    )

    write_tsv(
        folder / f"{dataset.lower()}_region_gene_map_50kb.tsv",
        rows, map_columns
    )

# --------------------------------------------------
# 4. Compare gene sets separately for each trait
# --------------------------------------------------

gene_regions = defaultdict(list)

for key, rows in mapping.items():

    trait, dataset = key

    for row in rows:

        gid = row.get("Gene_ID")

        if gid:
            gene_regions[(trait, dataset, gid)].append(
                row["Region_ID"]
            )

comparison = []

for trait in ["BB", "DBH", "FD"]:

    original_genes = {
        gid for t, d, gid in gene_regions
        if t == trait and d == "Original"
    }

    imputed_genes = {
        gid for t, d, gid in gene_regions
        if t == trait and d == "Imputed"
    }

    all_genes = original_genes | imputed_genes

    for gid in sorted(all_genes):

        in_original = gid in original_genes
        in_imputed = gid in imputed_genes

        if in_original and in_imputed:
            classification = "SHARED"

        elif in_imputed:
            classification = "IMPUTED_SPECIFIC"

        else:
            classification = "ORIGINAL_SPECIFIC"

        info = gene_info.get(gid, {})

        comparison.append({
            "Trait": trait,
            "Gene_ID": gid,
            "Classification": classification,
            "Original": "YES" if in_original else "NO",
            "Imputed": "YES" if in_imputed else "NO",
            "Original_association_regions": ";".join(
                gene_regions[(trait, "Original", gid)]
            ),
            "Imputed_association_regions": ";".join(
                gene_regions[(trait, "Imputed", gid)]
            ),
            "Chromosome": info.get("Chromosome", ""),
            "Gene_start": info.get("Gene_start", ""),
            "Gene_end": info.get("Gene_end", ""),
            "Gene_source": info.get("Gene_source", ""),
            "Gene_name": info.get("Gene_name", ""),
            "Gene_biotype": info.get("Gene_biotype", ""),
            "Gene_description": info.get("Gene_description", "")
        })

comparison_columns = list(comparison[0].keys())

outdir = WORK / "02_gene_annotation" / "gene_set_comparison_50kb"

write_tsv(
    outdir / "candidate_gene_set_comparison_by_trait.tsv",
    comparison, comparison_columns
)

imputed_specific = [
    r for r in comparison
    if r["Classification"] == "IMPUTED_SPECIFIC"
]

write_tsv(
    outdir / "imputed_specific_candidate_genes.tsv",
    imputed_specific, comparison_columns
)

# --------------------------------------------------
# 5. Summary
# --------------------------------------------------

summary = []

for trait in ["BB", "DBH", "FD"]:

    original = {
        gid for t, d, gid in gene_regions
        if t == trait and d == "Original"
    }

    imputed = {
        gid for t, d, gid in gene_regions
        if t == trait and d == "Imputed"
    }

    summary.append({
        "Trait": trait,
        "Original_regions": len(regions[(trait, "Original")]),
        "Imputed_regions": len(regions[(trait, "Imputed")]),
        "Original_candidate_genes": len(original),
        "Imputed_candidate_genes": len(imputed),
        "Shared_genes": len(original & imputed),
        "Original_specific_genes": len(original - imputed),
        "Imputed_specific_genes": len(imputed - original)
    })

write_tsv(
    outdir / "candidate_gene_set_summary.tsv",
    summary, list(summary[0].keys())
)

write_tsv(
    outdir / "nearest_gene_by_region.tsv",
    nearest_rows,
    list(nearest_rows[0].keys())
)

print("\n===== GENE SET COMPARISON =====")

for row in summary:
    print(row)

print("\nAnalysis completed.")
print("Output directory:", outdir)
