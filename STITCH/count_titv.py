#!/usr/bin/env python3
import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("ROOT", ".")).resolve()
OUT = ROOT / "titv_stats"
OUT.mkdir(exist_ok=True)

with (ROOT / "config/manifest.tsv").open(newline="") as f:
    manifest = list(csv.DictReader(f, delimiter="\t"))

TI = {b"AG", b"GA", b"CT", b"TC"}
BASES = {b"A", b"C", b"G", b"T"}

def count_chunk(task_id):
    if not 1 <= task_id <= len(manifest):
        raise SystemExit(f"Invalid task ID: {task_id}")

    chunk = manifest[task_id - 1]["chunk"].strip()
    source = ROOT / "work" / f"{chunk}.keys.tsv"
    result = OUT / f"{chunk}.tsv"

    if not source.is_file():
        raise SystemExit(f"Missing input: {source}")

    n = ti = tv = invalid = 0

    with source.open("rb", buffering=1024 * 1024) as f:
        for line in f:
            fields = line.rstrip(b"\r\n").split(b"\t")

            if len(fields) != 4:
                raise ValueError(
                    f"{chunk}: invalid column count at line {n+1}"
                )

            if fields[0] != chunk.encode():
                raise ValueError(
                    f"{chunk}: chromosome/contig mismatch at line {n+1}"
                )

            ref, alt = fields[2].upper(), fields[3].upper()
            n += 1

            if ref not in BASES or alt not in BASES or ref == alt:
                invalid += 1
            elif ref + alt in TI:
                ti += 1
            else:
                tv += 1

    with (ROOT / "per_file_v4" / f"{chunk}.json").open() as f:
        expected = json.load(f)["counts"]["matched_records"]

    if n != expected:
        raise SystemExit(
            f"ERROR {chunk}: keys={n}, matched_records={expected}"
        )

    if n != ti + tv + invalid:
        raise SystemExit(f"ERROR {chunk}: Ti/Tv count mismatch")

    tmp = result.with_name(result.name + f".{os.getpid()}.tmp")
    try:
        with tmp.open("w") as f:
            f.write("chunk\ttotal_snps\tTi\tTv\tinvalid\n")
            f.write(f"{chunk}\t{n}\t{ti}\t{tv}\t{invalid}\n")
        os.replace(tmp, result)
    finally:
        if tmp.exists():
            tmp.unlink()

    print(f"{chunk}: SNPs={n} Ti={ti} Tv={tv} invalid={invalid}")

def merge():
    totals = {"total_snps": 0, "Ti": 0, "Tv": 0, "invalid": 0}
    rows = []

    for entry in manifest:
        chunk = entry["chunk"].strip()
        path = OUT / f"{chunk}.tsv"

        if not path.is_file():
            raise SystemExit(f"ERROR: missing {path}")

        with path.open(newline="") as f:
            data = list(csv.DictReader(f, delimiter="\t"))

        if len(data) != 1 or data[0]["chunk"] != chunk:
            raise SystemExit(f"ERROR: invalid result for {chunk}")

        row = data[0]
        for k in totals:
            value = int(row[k])
            if value < 0:
                raise SystemExit(f"ERROR: negative {k} for {chunk}")
            totals[k] += value

        rows.append(row)

    expected = 2741098686
    if totals["total_snps"] != expected:
        raise SystemExit(
            f"ERROR: observed {totals['total_snps']}, expected {expected}"
        )

    if totals["total_snps"] != (
        totals["Ti"] + totals["Tv"] + totals["invalid"]
    ):
        raise SystemExit("ERROR: global Ti/Tv accounting mismatch")

    with (OUT / "per_chunk_titv.tsv").open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["chunk", "total_snps", "Ti", "Tv", "invalid"],
            delimiter="\t"
        )
        w.writeheader()
        w.writerows(rows)

    ratio = totals["Ti"] / totals["Tv"] if totals["Tv"] else float("nan")

    with (OUT / "genome_titv.tsv").open("w") as f:
        f.write("metric\tvalue\n")
        for k, v in totals.items():
            f.write(f"{k}\t{v}\n")
        f.write(f"Ti_Tv_ratio\t{ratio:.8f}\n")

    print("========== GENOME Ti/Tv ==========")
    for k, v in totals.items():
        print(f"{k}: {v}")
    print(f"Ti/Tv: {ratio:.8f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: count_titv.py TASK_ID | merge")
    if sys.argv[1] == "merge":
        merge()
    else:
        count_chunk(int(sys.argv[1]))
