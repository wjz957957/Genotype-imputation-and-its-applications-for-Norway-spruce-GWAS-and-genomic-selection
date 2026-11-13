#!/usr/bin/env python3

import csv
import argparse
from pathlib import Path
from collections import Counter, defaultdict

SOURCES = [
    ("BB", "Original", "GWAS-Bud",
     "gwas_clumps_revised_origin_Bud.clumped"),
    ("BB", "Imputed", "GWAS-Bud",
     "gwas_clumps_revised_imputed_Bud.clumped"),

    ("DBH", "Original", "GWAS-DBH",
     "gwas_clumps_revised_origin_DBH.clumped"),
    ("DBH", "Imputed", "GWAS-DBH",
     "gwas_clumps_revised_imputed_DBH.clumped"),

    ("FD", "Original", "GWAS-FD",
     "gwas_clumps_revised_origin_FD.clumped"),
    ("FD", "Imputed", "GWAS-FD",
     "gwas_clumps_revised_imputed_FD.clumped"),
]


def read_clumped(path, trait, dataset):

    records = []
    header = None

    with open(path) as f:

        for line_number, line in enumerate(f, 1):

            fields = line.strip().split()

            if not fields:
                continue

            if header is None:
                if all(x in fields for x in
                       ["CHR", "SNP", "BP", "P"]):
                    header = fields
                continue

            if len(fields) < len(header):
                raise ValueError(
                    f"Malformed row: {path}:{line_number}"
                )

            row = dict(zip(header, fields))

            chrom = row["CHR"]
            bp = int(row["BP"])
            snp = row["SNP"]

            if snp in ("NONE", ".", "NA"):
                raise ValueError(
                    f"Invalid index SNP: {path}:{line_number}"
                )

            records.append({
                "Region_ID":
                    f"{trait}_{dataset}_chr{chrom}_{bp}",
                "Trait": trait,
                "Dataset": dataset,
                "CHR": chrom,
                "BP": bp,
                "Lead_SNP": snp,
                "P": row["P"],
                "F": row.get("F", "NA"),
                "TOTAL": row.get("TOTAL", "NA"),
                "SP2": row.get("SP2", "NA"),
                "Source_file": path.name
            })

    if header is None:
        raise ValueError(f"No valid header: {path}")

    return records


def write_tsv(path, rows, columns):

    with open(path, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--root", required=True)
    parser.add_argument("--new-dir", required=True)

    args = parser.parse_args()

    root = Path(args.root)
    new = Path(args.new_dir)

    out = new / "01_regions"
    out.mkdir(parents=True, exist_ok=True)

    all_records = []

    for trait, dataset, subdir, filename in SOURCES:

        path = root / "GWAS" / subdir / filename

        if not path.is_file():
            raise FileNotFoundError(path)

        records = read_clumped(path, trait, dataset)

        all_records.extend(records)

        print(
            f"{trait:4s} {dataset:8s}: "
            f"{len(records)} association regions"
        )

    columns = [
        "Region_ID", "Trait", "Dataset",
        "CHR", "BP", "Lead_SNP", "P",
        "F", "TOTAL", "SP2", "Source_file"
    ]

    write_tsv(
        out / "all_index_snps.tsv",
        all_records,
        columns
    )

    # Region counts and duplicate positions
    summary = []
    duplicates = []

    grouped = defaultdict(list)

    for r in all_records:
        key = (r["Trait"], r["Dataset"])
        grouped[key].append(r)

    for key, records in grouped.items():

        trait, dataset = key

        positions = [
            (r["CHR"], r["BP"]) for r in records
        ]

        counts = Counter(positions)

        summary.append({
            "Trait": trait,
            "Dataset": dataset,
            "N_index_SNPs": len(records),
            "N_unique_positions": len(counts),
            "N_duplicate_positions":
                sum(v - 1 for v in counts.values() if v > 1)
        })

        for pos, count in counts.items():

            if count > 1:
                duplicates.append({
                    "Trait": trait,
                    "Dataset": dataset,
                    "CHR": pos[0],
                    "BP": pos[1],
                    "Count": count
                })

    write_tsv(
        out / "region_counts.tsv",
        summary,
        [
            "Trait", "Dataset",
            "N_index_SNPs",
            "N_unique_positions",
            "N_duplicate_positions"
        ]
    )

    write_tsv(
        out / "duplicate_index_positions.tsv",
        duplicates,
        ["Trait", "Dataset", "CHR", "BP", "Count"]
    )

    # Exact positional overlap between Original and Imputed
    overlap = []

    for trait in ["BB", "DBH", "FD"]:

        original = {
            (r["CHR"], r["BP"])
            for r in grouped[(trait, "Original")]
        }

        imputed = {
            (r["CHR"], r["BP"])
            for r in grouped[(trait, "Imputed")]
        }

        overlap.append({
            "Trait": trait,
            "Original_regions": len(original),
            "Imputed_regions": len(imputed),
            "Shared_positions": len(original & imputed),
            "Original_only": len(original - imputed),
            "Imputed_only": len(imputed - original)
        })

    write_tsv(
        out / "index_snp_overlap_summary.tsv",
        overlap,
        [
            "Trait", "Original_regions",
            "Imputed_regions", "Shared_positions",
            "Original_only", "Imputed_only"
        ]
    )

    print("\nExtraction completed.")
    print("Output:", out)
    print("Total index SNP records:", len(all_records))


if __name__ == "__main__":
    main()
