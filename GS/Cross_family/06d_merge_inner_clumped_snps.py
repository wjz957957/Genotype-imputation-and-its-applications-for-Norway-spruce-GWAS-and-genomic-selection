#!/usr/bin/env python3
"""Merge the 12 chromosome-wise clumped index SNP lists into fold-wide P-ranked pools."""
import argparse
import math
from pathlib import Path
import os


def need(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty input: {path}")
    return path


def read_ranked(path: Path, expected_chr: int):
    need(path)
    rows = []
    with path.open() as fh:
        header = fh.readline().strip().split()
        if header != ["SNP", "P"]:
            raise RuntimeError(f"Unexpected header in {path}: {header}")
        seen = set()
        for lineno, line in enumerate(fh, 2):
            f = line.split()
            if len(f) != 2:
                raise RuntimeError(f"Malformed row {path}:{lineno}")
            sid = f[0]
            try:
                p = float(f[1])
            except ValueError as e:
                raise RuntimeError(f"Invalid P-value {path}:{lineno}") from e
            if not math.isfinite(p) or not 0 <= p <= 1:
                raise RuntimeError(f"P-value outside [0,1] {path}:{lineno}")
            if sid in seen:
                raise RuntimeError(f"Duplicate SNP within {path}: {sid}")
            seen.add(sid)
            rows.append((sid, p, expected_chr))
    if not rows:
        raise RuntimeError(f"No clumped index SNPs in {path}")
    return rows


def atomic_write(path: Path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w") as fh:
        for line in lines:
            fh.write(line + "\n")
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--outer", required=True)
    ap.add_argument("--inner", required=True)
    a = ap.parse_args()
    work = Path(a.work)
    outer = f"{int(a.outer):02d}"
    inner = f"{int(a.inner):02d}"
    fold = f"outer_{outer}/inner_{inner}"
    by_chr_root = work / "Cross_family/results/06_inner_imputed_clumping_by_chr" / fold
    out_root = work / "Cross_family/results/06_inner_imputed_gwas_clumping" / fold
    all_summary = []
    for trait in ("BB", "H6"):
        rows = []
        for chr_num in range(1, 13):
            pad = f"{chr_num:02d}"
            chr_dir = by_chr_root / f"chr_{pad}"
            flag = need(chr_dir / "complete.flag")
            flag_text = flag.read_text().splitlines()
            if "status=PASS" not in flag_text:
                raise RuntimeError(f"Chromosome clumping flag is not PASS: {flag}")
            p = chr_dir / f"{trait}_chr{pad}.index_ranked.tsv"
            rows.extend(read_ranked(p, chr_num))
        ids = [r[0] for r in rows]
        if len(ids) != len(set(ids)):
            raise RuntimeError(f"Duplicate SNP IDs across chromosome lists for {fold}/{trait}")
        rows.sort(key=lambda r: (r[1], r[0]))
        out = out_root / trait
        atomic_write(out / "ranked_extra_snps.tsv", ["SNP\tP"] + [f"{sid}\t{p:.12g}" for sid, p, _ in rows])
        atomic_write(out / "rank_summary.tsv", [
            "Item\tValue",
            f"outer\t{outer}", f"inner\t{inner}", f"trait\t{trait}",
            "gwas_model\tGCTA_MLMA_with_mlma_subtract_grm",
            "clump_p1\t1", "clump_p2\t1", "clump_r2\t0.2", "clump_kb\t1000",
            f"chromosomes_merged\t12", f"ranked_clump_index_snps\t{len(rows)}",
            f"ranked_file\t{out / 'ranked_extra_snps.tsv'}"
        ])
        atomic_write(out / "complete.flag", [
            "version=2", f"outer={outer}", f"inner={inner}", f"trait={trait}",
            "gwas_model=GCTA_MLMA_with_mlma_subtract_grm", "clumping=PASS",
            "chromosomes=12/12", f"ranked_clump_index_snps={len(rows)}",
            "training_only=PASS", "status=PASS"
        ])
        all_summary.append((trait, len(rows)))
        print(f"PASS {fold}/{trait}: merged {len(rows)} clump index SNPs")
    atomic_write(out_root / "complete.flag", [
        "version=2", f"outer={outer}", f"inner={inner}", "traits=BB,H6",
        "chromosome_clumping=12/12", "ranking=GWAS_P_ascending",
        "training_only=PASS", "status=PASS"
    ])
    print(f"OVERALL PASS {fold}: {all_summary}")

if __name__ == "__main__":
    main()
