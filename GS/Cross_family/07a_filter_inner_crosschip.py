#!/usr/bin/env python3
"""Cross-chip LD filtering for inner-CV imputation-derived QTN candidates.

Workflow:
  rank:    remove exact overlap with the selected chip, preserving GWAS P order
  prepare: stage the top candidate pool and inner-training sample keep file
  report:  remove candidate SNPs with r2 >= threshold to any selected chip SNP
  publish: publish all QTN-count marker sets supported by the retained pool

Only inner-training genotypes are used to estimate LD. Inner-validation and
outer-test samples are not used for LD filtering. Sparse PLINK LD files report
only pairs at or above the requested threshold; absence from those files means
no reported pair reached the threshold, not that r2 is exactly zero.
"""
import argparse
import gzip
import math
import os
import shutil
from pathlib import Path

COUNTS = (1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500, 600, 700)


def required(path):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {p}")
    return p


def write(path, lines):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="") as fh:
        for line in lines:
            fh.write(str(line) + "\n")
    os.replace(tmp, path)


def read_ids(path):
    p = required(path)
    out = []
    seen = set()
    with p.open() as fh:
        for line_no, line in enumerate(fh, 1):
            fields = line.split()
            if len(fields) != 1:
                raise RuntimeError(f"Expected one ID per line at {p}:{line_no}")
            sid = fields[0]
            if sid in seen:
                raise RuntimeError(f"Duplicate ID {sid} in {p}")
            seen.add(sid)
            out.append(sid)
    return out


def read_train_pheno(path):
    p = required(path)
    ids = []
    seen = set()
    with p.open() as fh:
        for line_no, line in enumerate(fh, 1):
            f = line.split()
            if len(f) != 3:
                raise RuntimeError(f"Expected headerless FID IID Y at {p}:{line_no}")
            if f[0] != f[1]:
                raise RuntimeError(f"FID/IID mismatch at {p}:{line_no}: {f[0]} {f[1]}")
            if f[1] in seen:
                raise RuntimeError(f"Duplicate sample ID in {p}: {f[1]}")
            try:
                y = float(f[2])
            except ValueError:
                raise RuntimeError(f"Non-numeric phenotype at {p}:{line_no}")
            if not math.isfinite(y):
                raise RuntimeError(f"Non-finite phenotype at {p}:{line_no}")
            seen.add(f[1])
            ids.append(f[1])
    if not ids:
        raise RuntimeError(f"No samples in phenotype file: {p}")
    return ids


def read_ranked(path, limit=None):
    p = required(path)
    rows = []
    seen = set()
    prev = -1.0
    with p.open() as fh:
        header = fh.readline().split()
        if [x.upper() for x in header] != ["SNP", "P"]:
            raise RuntimeError(f"Expected header SNP P in {p}; found {header}")
        for line_no, line in enumerate(fh, 2):
            f = line.split()
            if len(f) != 2:
                raise RuntimeError(f"Malformed ranking row at {p}:{line_no}")
            sid = f[0]
            try:
                pv = float(f[1])
            except ValueError:
                raise RuntimeError(f"Invalid P value at {p}:{line_no}")
            if not math.isfinite(pv) or pv < 0 or pv > 1:
                raise RuntimeError(f"Out-of-range P value at {p}:{line_no}")
            if pv + 1e-15 < prev:
                raise RuntimeError(f"Ranking is not sorted by ascending P at {p}:{line_no}")
            if sid in seen:
                raise RuntimeError(f"Duplicate ranked SNP at {p}:{line_no}: {sid}")
            seen.add(sid)
            rows.append((sid, pv, f[1]))
            prev = pv
            if limit is not None and len(rows) >= limit:
                break
    return rows


def rank_nonchip(a):
    chip = set(read_ids(a.chip))
    ranked_path = required(a.ranked)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    n_in = n_overlap = n_out = 0
    seen = set()
    prev = -1.0
    with ranked_path.open() as src, tmp.open("w") as dst:
        header = src.readline().split()
        if [x.upper() for x in header] != ["SNP", "P"]:
            raise RuntimeError(f"Expected header SNP P in {ranked_path}; found {header}")
        dst.write("SNP\tP\n")
        for line_no, line in enumerate(src, 2):
            f = line.split()
            if len(f) != 2:
                raise RuntimeError(f"Malformed ranking row at {ranked_path}:{line_no}")
            sid = f[0]
            try:
                pv = float(f[1])
            except ValueError:
                raise RuntimeError(f"Invalid P value at {ranked_path}:{line_no}")
            if not math.isfinite(pv) or pv < 0 or pv > 1:
                raise RuntimeError(f"Out-of-range P value at {ranked_path}:{line_no}")
            if pv + 1e-15 < prev:
                raise RuntimeError(f"Ranking is not sorted by P at {ranked_path}:{line_no}")
            if sid in seen:
                raise RuntimeError(f"Duplicate SNP in ranking: {sid}")
            seen.add(sid)
            prev = pv
            n_in += 1
            if sid in chip:
                n_overlap += 1
                continue
            dst.write(f"{sid}\t{f[1]}\n")
            n_out += 1
    os.replace(tmp, out)
    summary = Path(a.summary)
    write(summary, ["Item\tValue", f"Input_clump_indices\t{n_in}",
                   f"Exact_chip_overlaps_removed\t{n_overlap}",
                   f"Nonchip_candidates_retained\t{n_out}"])
    if n_out == 0:
        raise RuntimeError("No non-chip GWAS clump candidates remain")
    print(f"PASS: {n_in} clump indices; removed {n_overlap} exact chip overlaps; {n_out} candidates remain")


def prepare(a):
    chip = read_ids(a.chip)
    train_ids = read_train_pheno(a.train_pheno)
    rows = read_ranked(a.ranked, limit=a.pool)
    if len(rows) < a.pool:
        # The requested pool is expected to be capped to available candidates by the submit script.
        raise RuntimeError(f"Requested pool {a.pool}, but only {len(rows)} ranked candidates are available")
    candidate_ids = [r[0] for r in rows]
    if set(candidate_ids) & set(chip):
        raise RuntimeError("Exact chip overlaps remain in ranked_nonchip input")
    if len(set(chip + candidate_ids)) != len(chip) + len(candidate_ids):
        raise RuntimeError("Duplicate SNP in chip plus QTN union")

    d = Path(a.stage)
    d.mkdir(parents=True, exist_ok=True)
    write(d / "selected_chip.snps", chip)
    write(d / "requested_qtn.snps", candidate_ids)
    write(d / "union.snps", chip + candidate_ids)
    write(d / "candidate_rank.tsv", ["Original_Rank\tSNP\tP"] +
          [f"{i}\t{sid}\t{pstr}" for i, (sid, pv, pstr) in enumerate(rows, 1)])
    write(d / "train.keep", [f"{sid}\t{sid}" for sid in train_ids])
    print(f"[PREP] chip={len(chip)} candidates={len(candidate_ids)} train_n={len(train_ids)} stage={d}")


def read_bim(path):
    out = {}
    with required(path).open() as fh:
        for line_no, line in enumerate(fh, 1):
            f = line.split()
            if len(f) != 6:
                raise RuntimeError(f"Malformed BIM line {path}:{line_no}")
            if f[1] in out:
                raise RuntimeError(f"Duplicate SNP in BIM: {f[1]}")
            out[f[1]] = f[0]
    return out


def report(a):
    d = Path(a.stage)
    chip = read_ids(d / "selected_chip.snps")
    train_ids = read_train_pheno(a.train_pheno)
    rows = read_ranked(a.ranked, limit=a.pool)
    staged = read_ids(d / "requested_qtn.snps")
    if staged != [sid for sid, pv, pstr in rows]:
        raise RuntimeError("Candidate ranking changed after LD genotype extraction")

    bim_chr = read_bim(a.bim)
    expected = set(chip) | set(staged)
    missing = expected - set(bim_chr)
    extra = set(bim_chr) - expected
    if missing or extra:
        raise RuntimeError(f"LD BIM does not match requested union: missing={len(missing)} {sorted(missing)[:5]}, extra={len(extra)} {sorted(extra)[:5]}")

    fam_ids = []
    with required(a.fam).open() as fh:
        for line_no, line in enumerate(fh, 1):
            f = line.split()
            if len(f) < 2 or f[0] != f[1]:
                raise RuntimeError(f"Invalid FID/IID in LD FAM at line {line_no}")
            fam_ids.append(f[1])
    if len(fam_ids) != len(train_ids) or set(fam_ids) != set(train_ids) or len(set(fam_ids)) != len(fam_ids):
        raise RuntimeError("LD genotype FAM samples do not match inner-training phenotype IDs exactly")

    chip_set = set(chip)
    cand_set = set(staged)
    chip_by_chr = {}
    for sid in chip:
        chip_by_chr[bim_chr[sid]] = chip_by_chr.get(bim_chr[sid], 0) + 1

    max_r2 = {}
    best_chip = {}
    n_pairs = 0
    for ld_name in a.ld:
        ld = required(ld_name)
        opener = gzip.open if ld.suffix == ".gz" else open
        with opener(ld, "rt") as fh:
            header = fh.readline().split()
            if not {"SNP_A", "SNP_B", "R2"}.issubset(header):
                raise RuntimeError(f"Unexpected PLINK LD header in {ld}: {header}")
            ia, ib, ir = header.index("SNP_A"), header.index("SNP_B"), header.index("R2")
            for line_no, line in enumerate(fh, 2):
                f = line.split()
                if len(f) <= max(ia, ib, ir):
                    raise RuntimeError(f"Malformed LD row {ld}:{line_no}")
                sa, sb = f[ia], f[ib]
                if sa in chip_set and sb in cand_set:
                    chip_id, qtn_id = sa, sb
                elif sb in chip_set and sa in cand_set:
                    chip_id, qtn_id = sb, sa
                else:
                    continue
                if bim_chr[chip_id] != bim_chr[qtn_id]:
                    raise RuntimeError(f"Unexpected cross-chromosome LD pair: {chip_id} {qtn_id}")
                try:
                    r2 = float(f[ir])
                except ValueError:
                    raise RuntimeError(f"Invalid R2 in {ld}:{line_no}")
                if not math.isfinite(r2) or r2 < 0 or r2 > 1.00001:
                    raise RuntimeError(f"Out-of-range R2 in {ld}:{line_no}: {r2}")
                r2 = min(1.0, r2)
                if r2 + 1e-8 < a.threshold:
                    raise RuntimeError(f"Sparse LD output contains R2={r2} below threshold={a.threshold}")
                n_pairs += 1
                if qtn_id not in max_r2 or r2 > max_r2[qtn_id] or (r2 == max_r2[qtn_id] and chip_id < best_chip[qtn_id]):
                    max_r2[qtn_id] = r2
                    best_chip[qtn_id] = chip_id

    annotated = ["Original_Rank\tSNP\tP\tChr\tLD_status\tMax_reported_high_r2\tTop_linked_CHIP\tKeep\tFiltered_Rank"]
    retained = []
    no_same_chr = 0
    for i, (sid, pv, pstr) in enumerate(rows, 1):
        if chip_by_chr.get(bim_chr[sid], 0) == 0:
            status, r2s, linked, keep = "NO_SELECTED_CHIP_ON_CHR", "NA", "NA", True
            no_same_chr += 1
        elif sid in max_r2:
            status, r2s, linked, keep = "FILTERED_HIGH_LD", f"{max_r2[sid]:.7f}", best_chip[sid], False
        else:
            status, r2s, linked, keep = "BELOW_THRESHOLD", f"<{a.threshold:g}", "NA", True
        if keep:
            retained.append((sid, pstr))
        annotated.append(f"{i}\t{sid}\t{pstr}\t{bim_chr[sid]}\t{status}\t{r2s}\t{linked}\t{int(keep)}\t{len(retained) if keep else 'NA'}")

    write(d / "candidate_ld_annotated.tsv", annotated)
    write(d / "filtered_ranked_extra_snps.tsv", ["SNP\tP"] + [f"{sid}\t{p}" for sid, p in retained])
    write(d / "excluded_high_ld.snps", [sid for sid, pv, pstr in rows if sid in max_r2])
    write(d / "filter_summary.tsv", [
        "Item\tValue", f"Candidate_pool_size\t{len(rows)}", f"CHIP_count\t{len(chip)}",
        f"Inner_training_size\t{len(train_ids)}", f"LD_threshold\t{a.threshold}",
        f"Filtered_high_LD\t{len(rows) - len(retained)}", f"Retained\t{len(retained)}",
        f"No_selected_chip_on_chr\t{no_same_chr}", f"High_LD_pairs_reported\t{n_pairs}",
        f"Can_supply_700_QTN\t{int(len(retained) >= 700)}"])
    write(d / "stage.complete", [f"pool={len(rows)}\tretained={len(retained)}\tthreshold={a.threshold}"])
    print(f"[RESULT] pool={len(rows)} retained={len(retained)} filtered={len(rows)-len(retained)} high_LD_pairs={n_pairs}")


def publish(a):
    stage = Path(a.stage)
    required(stage / "stage.complete")
    rows = read_ranked(stage / "filtered_ranked_extra_snps.tsv")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    marker_dir = out / "marker_sets"
    marker_dir.mkdir(parents=True, exist_ok=True)
    for old in marker_dir.glob("extra_*.snps"):
        old.unlink()

    write(out / "filtered_ranked_extra_snps.tsv", ["SNP\tP"] + [f"{sid}\t{pstr}" for sid, pv, pstr in rows])
    for name in ("candidate_ld_annotated.tsv", "filter_summary.tsv", "excluded_high_ld.snps"):
        src = required(stage / name)
        shutil.copyfile(src, out / name)

    write(marker_dir / "extra_0.snps", [])
    availability = ["N_QTN\tAvailable\tReason"]
    for n in (0,) + COUNTS:
        ok = n <= len(rows)
        availability.append(f"{n}\t{int(ok)}\t{'AVAILABLE' if ok else 'INSUFFICIENT_RETAINED_CANDIDATES'}")
        if n > 0 and ok:
            write(marker_dir / f"extra_{n}.snps", [sid for sid, pv, pstr in rows[:n]])
    write(out / "candidate_count_availability.tsv", availability)

    meta = {
        "outer": a.outer, "inner": a.inner, "trait": a.trait,
        "chip_strategy": a.chip_set, "ld_threshold": str(a.threshold),
        "retained_candidates": str(len(rows)), "max_qtn_available": str(max([n for n in COUNTS if n <= len(rows)] + [0])),
        "all_700_available": "YES" if len(rows) >= 700 else "NO",
        "training_only_ld": "PASS", "status": "PASS"
    }
    write(out / "complete.flag", [f"{k}={v}" for k, v in meta.items()])
    print(f"[PUBLISHED] {out}; retained={len(rows)}; max_qtn_available={meta['max_qtn_available']}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="mode", required=True)
    x = sub.add_parser("rank")
    x.add_argument("--ranked", required=True)
    x.add_argument("--chip", required=True)
    x.add_argument("--out", required=True)
    x.add_argument("--summary", required=True)

    y = sub.add_parser("prepare")
    y.add_argument("--ranked", required=True)
    y.add_argument("--chip", required=True)
    y.add_argument("--train-pheno", required=True)
    y.add_argument("--pool", type=int, required=True)
    y.add_argument("--stage", required=True)

    z = sub.add_parser("report")
    z.add_argument("--ranked", required=True)
    z.add_argument("--stage", required=True)
    z.add_argument("--bim", required=True)
    z.add_argument("--fam", required=True)
    z.add_argument("--train-pheno", required=True)
    z.add_argument("--ld", required=True, action="append")
    z.add_argument("--pool", type=int, required=True)
    z.add_argument("--threshold", type=float, default=0.8)

    q = sub.add_parser("publish")
    q.add_argument("--stage", required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--outer", required=True)
    q.add_argument("--inner", required=True)
    q.add_argument("--trait", required=True)
    q.add_argument("--chip-set", required=True)
    q.add_argument("--threshold", type=float, default=0.8)

    a = p.parse_args()
    if a.mode in ("prepare", "report") and a.pool < 1:
        p.error("pool must be >= 1")
    if a.mode == "report" and not 0 < a.threshold <= 1:
        p.error("threshold must be in (0,1]")
    if a.mode == "publish" and not 0 < a.threshold <= 1:
        p.error("threshold must be in (0,1]")
    {"rank": rank_nonchip, "prepare": prepare, "report": report, "publish": publish}[a.mode](a)


if __name__ == "__main__":
    main()
