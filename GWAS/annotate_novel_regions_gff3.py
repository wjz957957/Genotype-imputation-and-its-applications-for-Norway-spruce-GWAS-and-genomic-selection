#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
from collections import defaultdict


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=columns,
            delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_attributes(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def load_gene_annotation(gff_path):
    genes = defaultdict(list)

    with gzip.open(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9 or fields[2] != "gene":
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])

            gene_id = attrs.get("ID")

            if not gene_id:
                continue

            genes[chrom].append({
                "Gene_ID": gene_id,
                "Chromosome": chrom,
                "Gene_start": start,
                "Gene_end": end,
                "Gene_strand": strand,
                "Gene_source": fields[1],
                "Gene_length_bp": end - start + 1
            })

    for chrom in genes:
        genes[chrom].sort(key=lambda x: x["Gene_start"])

    return genes


def load_longest_gene_ids(gff_path):
    ids = set()

    with gzip.open(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9 or fields[2] != "gene":
                continue

            gene_id = parse_attributes(fields[8]).get("ID")

            if gene_id:
                ids.add(gene_id)

    return ids


def gene_distance(gene, position):
    if gene["Gene_start"] <= position <= gene["Gene_end"]:
        return 0
    if position < gene["Gene_start"]:
        return gene["Gene_start"] - position
    return position - gene["Gene_end"]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--regions", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--longest-gff", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window", type=int, default=50000)

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading complete Norway spruce v2 GFF3...")
    genes = load_gene_annotation(args.gff)

    print("Loading longest_no_TE gene IDs...")
    longest_ids = load_longest_gene_ids(args.longest_gff)

    total_genes = sum(len(x) for x in genes.values())

    print("Total genes:", total_genes)
    print("Genes in longest_no_TE:", len(longest_ids))

    regions = read_tsv(args.regions)

    expected = {"BB": 9, "DBH": 2, "FD": 5}
    observed = defaultdict(int)

    for region in regions:
        observed[region["Trait"]] += 1

    if len(regions) != 16 or dict(observed) != expected:
        raise ValueError(
            f"Unexpected region counts: {dict(observed)}"
        )

    candidate_rows = []
    nearest_rows = []
    summary_rows = []

    for region in regions:

        trait = region["Trait"]
        lead = region["Imputed index SNP"]
        chrom_num = int(region["Chr"])
        position = int(region["Position (bp)"])

        chrom = f"PA_chr{chrom_num:02d}"

        if lead != f"{chrom_num}:{position}":
            raise ValueError(
                f"Lead SNP and coordinate mismatch: {lead}"
            )

        if chrom not in genes:
            raise ValueError(
                f"Chromosome not found in GFF3: {chrom}"
            )

        window_start = max(1, position - args.window)
        window_end = position + args.window

        chromosome_genes = genes[chrom]

        # All genes overlapping the +/-50 kb window
        candidates = [
            g for g in chromosome_genes
            if g["Gene_start"] <= window_end
            and g["Gene_end"] >= window_start
        ]

        # Distance to every gene on this chromosome
        distances = [
            (gene_distance(g, position), g)
            for g in chromosome_genes
        ]

        minimum_distance = min(d[0] for d in distances)

        nearest = [
            g for d, g in distances
            if d == minimum_distance
        ]

        for gene in candidates:

            distance = gene_distance(gene, position)

            if distance == 0:
                relation = "OVERLAPS_LEAD_SNP"
            elif gene["Gene_end"] < position:
                relation = "GENOMIC_LEFT_OF_LEAD"
            else:
                relation = "GENOMIC_RIGHT_OF_LEAD"

            row = {
                "Trait": trait,
                "Lead_SNP": lead,
                "Lead_P": region["Index P"],
                "Chromosome": chrom,
                "Lead_position": position,
                "Window_start": window_start,
                "Window_end": window_end,
                **gene,
                "Distance_to_lead_bp": distance,
                "Lead_gene_relationship": relation,
                "Present_in_longest_no_TE": (
                    "YES" if gene["Gene_ID"] in longest_ids
                    else "NO"
                )
            }

            candidate_rows.append(row)

        for gene in nearest:
            nearest_rows.append({
                "Trait": trait,
                "Lead_SNP": lead,
                "Lead_P": region["Index P"],
                "Chromosome": chrom,
                "Lead_position": position,
                **gene,
                "Distance_to_lead_bp": minimum_distance,
                "Present_in_longest_no_TE": (
                    "YES" if gene["Gene_ID"] in longest_ids
                    else "NO"
                )
            })

        summary_rows.append({
            "Trait": trait,
            "Lead_SNP": lead,
            "Lead_P": region["Index P"],
            "Chromosome": chrom,
            "Lead_position": position,
            "Search_window_bp": args.window,
            "N_genes_within_window": len(candidates),
            "Nearest_gene_ID": ";".join(
                g["Gene_ID"] for g in nearest
            ),
            "Nearest_gene_distance_bp": minimum_distance,
            "Nearest_gene_in_window": (
                "YES" if nearest[0] in candidates else "NO"
            )
        })

        print(
            f"{trait} {lead}: "
            f"{len(candidates)} genes; "
            f"nearest distance={minimum_distance} bp"
        )

    candidate_columns = [
        "Trait", "Lead_SNP", "Lead_P", "Chromosome",
        "Lead_position", "Window_start", "Window_end",
        "Gene_ID", "Gene_start", "Gene_end",
        "Gene_strand", "Gene_source", "Gene_length_bp",
        "Distance_to_lead_bp", "Lead_gene_relationship",
        "Present_in_longest_no_TE"
    ]

    nearest_columns = [
        "Trait", "Lead_SNP", "Lead_P", "Chromosome",
        "Lead_position", "Gene_ID", "Gene_start",
        "Gene_end", "Gene_strand", "Gene_source",
        "Gene_length_bp", "Distance_to_lead_bp",
        "Present_in_longest_no_TE"
    ]

    summary_columns = [
        "Trait", "Lead_SNP", "Lead_P", "Chromosome",
        "Lead_position", "Search_window_bp",
        "N_genes_within_window", "Nearest_gene_ID",
        "Nearest_gene_distance_bp", "Nearest_gene_in_window"
    ]

    write_tsv(
        os.path.join(args.outdir,
                     "new_region_gene_candidates_50kb.tsv"),
        candidate_rows, candidate_columns
    )

    write_tsv(
        os.path.join(args.outdir,
                     "new_region_nearest_gene.tsv"),
        nearest_rows, nearest_columns
    )

    write_tsv(
        os.path.join(args.outdir,
                     "new_region_gene_annotation_summary.tsv"),
        summary_rows, summary_columns
    )

    print("\nGene annotation completed.")
    print("Regions:", len(regions))
    print("Candidate gene records:", len(candidate_rows))
    print("Nearest gene records:", len(nearest_rows))
    print("Output:", args.outdir)


if __name__ == "__main__":
    main()
