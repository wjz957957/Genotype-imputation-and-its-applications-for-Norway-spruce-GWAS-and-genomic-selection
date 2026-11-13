#!/usr/bin/env python3
"""BB outer_01 top1000 QTN vs selected CHIP LD filtering. No GWAS/model fitting."""
import argparse
import math
from pathlib import Path

COUNTS = (1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75,
          100, 150, 200, 300, 400, 500, 750, 1000)


def require_path(path):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty input: {p}")
    return p


def read_chip(path):
    path = require_path(path)
    chip = []
    with path.open(encoding="utf-8") as fh:
        for i, line in enumerate(fh, 1):
            val = line.strip().split()
            if len(val) != 1:
                raise RuntimeError(f"CHIP SNP list must have exactly one ID per line: {path}:{i}")
            chip.append(val[0])
    if len(chip) != 100 or len(chip) != len(set(chip)):
        raise RuntimeError(f"Expected exactly 100 unique CHIP markers; got {len(chip)} entries")
    return chip


def read_ranked(path, count=1000):
    path = require_path(path)
    rows = []
    seen = set()
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().split()
        if header != ["SNP", "P"]:
            raise RuntimeError(f"Expected 'SNP P' ranked header, got {header}: {path}")
        previous_p = -1.0
        for i, line in enumerate(fh, 2):
            if len(rows) >= count:
                break
            fields = line.split()
            if len(fields) != 2:
                raise RuntimeError(f"Invalid ranked row at line {i}: {line!r}")
            sid, raw_p = fields
            try:
                p = float(raw_p)
            except ValueError:
                raise RuntimeError(f"Non-numeric P at line {i}: {raw_p}")
            if not math.isfinite(p) or not (0 <= p <= 1):
                raise RuntimeError(f"Invalid P at line {i}: {raw_p}")
            if p < previous_p:
                raise RuntimeError(f"Ranked file is not P-ascending at line {i}")
            if sid in seen:
                raise RuntimeError(f"Duplicate QTN in top {count}: {sid}")
            seen.add(sid)
            previous_p = p
            rows.append((sid, raw_p))
    if len(rows) != count:
        raise RuntimeError(f"Expected {count} candidate QTN; got {len(rows)}")
    return rows


def output_lines(path, lines):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_name(p.name + ".tmp")
    with t.open("w", encoding="utf-8", newline="\n") as fh:
        for line in lines:
            fh.write(line + "\n")
    t.replace(p)


def prep(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    chip = read_chip(args.chip)
    top = read_ranked(args.ranked)
    candidates = [s for s, _ in top]
    intersect = sorted(set(chip) & set(candidates))
    if intersect:
        raise RuntimeError(f"CHIP and imputed top1000 overlap: {intersect[:8]}")

    output_lines(out / "selected_chip_100.snps", chip)
    output_lines(out / "original_top1000.snps", candidates)
    output_lines(out / "original_top1000_ranked.tsv",
                 ["Original_Rank\tSNP\tP"] +
                 [f"{i}\t{s}\t{p}" for i, (s, p) in enumerate(top, 1)])
    output_lines(out / "all_1100_for_LD.snps", chip + candidates)
    print(f"[PREP] 100 CHIP, 1000 imputed candidates, unique=1100; GWAS source: {args.ranked}")


def read_bim(path):
    by_snp = {}
    with require_path(path).open(encoding="utf-8") as fh:
        for i, line in enumerate(fh, 1):
            fields = line.split()
            if len(fields) != 6:
                raise RuntimeError(f"Invalid BIM line {i}: {path}")
            chrom, sid = fields[0], fields[1]
            if sid in by_snp:
                raise RuntimeError(f"Duplicate BIM SNP: {sid}")
            by_snp[sid] = chrom
    return by_snp


def report(args):
    out = Path(args.out)
    chip = read_chip(out / "selected_chip_100.snps")
    top = read_ranked(args.ranked)
    qtn_ids = [s for s, _ in top]
    chip_set = set(chip)
    qtn_set = set(qtn_ids)

    bim = read_bim(args.bim)
    missing = (chip_set | qtn_set) - bim.keys()
    unexpected = set(bim) - chip_set - qtn_set
    if missing or unexpected:
        raise RuntimeError(f"LD data must contain exactly 1100 expected markers. "
                           f"missing={len(missing)} {sorted(missing)[:6]} "
                           f"unexpected={len(unexpected)} {sorted(unexpected)[:6]}")
    n_fam = sum(1 for line in require_path(args.fam).open() if line.strip())
    if n_fam != 813:
        raise RuntimeError(f"Expected 813 outer-training individuals, got {n_fam}")

    chip_chr_counts = {}
    for sid in chip:
        chrom = bim[sid]
        chip_chr_counts[chrom] = chip_chr_counts.get(chrom, 0) + 1

    maxr2 = {}
    linked = {}
    pair_counts = {}
    with require_path(args.ld).open(encoding="utf-8") as fh:
        header = fh.readline().split()
        try:
            ai, bi, ri = header.index("SNP_A"), header.index("SNP_B"), header.index("R2")
        except ValueError:
            raise RuntimeError(f"Unexpected PLINK --r2 header: {header}")
        for ln, line in enumerate(fh, 2):
            row = line.split()
            if len(row) <= max(ai, bi, ri):
                raise RuntimeError(f"Bad PLINK LD row at {ln}")
            a, b = row[ai], row[bi]
            if a in chip_set and b in qtn_set:
                q, c = b, a
            elif b in chip_set and a in qtn_set:
                q, c = a, b
            else:
                continue
            if bim[q] != bim[c]:
                raise RuntimeError(f"PLINK returned interchromosome pair: {q}, {c}")
            try:
                r = float(row[ri])
            except ValueError:
                raise RuntimeError(f"Unusable CHIP/QTN r2 at line {ln}: {row[ri]}")
            if not math.isfinite(r) or r < 0 or r > 1.00001:
                raise RuntimeError(f"Invalid CHIP/QTN r2 at line {ln}: {row[ri]}")
            r = min(1.0, r)
            pair_counts[q] = pair_counts.get(q, 0) + 1
            if q not in maxr2 or r > maxr2[q] or (r == maxr2[q] and c < linked[q]):
                maxr2[q] = r
                linked[q] = c

    # --r2 with --ld-window-r2 0 and a chromosome-spanning window must output
    # every same-chromosome CHIP/QTN pair. Never silently call absent data r2=0.
    for q in qtn_ids:
        expected = chip_chr_counts.get(bim[q], 0)
        observed = pair_counts.get(q, 0)
        if expected != observed:
            raise RuntimeError(
                f"Missing CHIP/QTN LD pairs: {q} chr={bim[q]} "
                f"expected={expected}, observed={observed}. "
                "Check LD-window arguments or genotype QC."
            )

    threshold = args.threshold
    all_lines = ["Original_Rank\tSNP\tP\tChromosome\tMax_r2_with_selected_CHIP\tTop_linked_CHIP\tLD_status\tKeep_r2_lt_threshold\tFiltered_Rank"]
    kept = []
    n_removed = 0
    n_no_chip = 0
    no_chip_chrs = set()
    for i, (sid, p) in enumerate(top, 1):
        if sid not in maxr2:
            status = "NO_SELECTED_CHIP_ON_CHR"
            maxstr = "NA"
            chipstr = "NA"
            keep = True
            n_no_chip += 1
            no_chip_chrs.add(bim[sid])
        else:
            status = "EVALUATED"
            maxstr = f"{maxr2[sid]:.7f}"
            chipstr = linked[sid]
            keep = maxr2[sid] < threshold
        if keep:
            kept.append((sid, p, i, maxstr, chipstr))
            new_rank = str(len(kept))
        else:
            n_removed += 1
            new_rank = "NA"
        all_lines.append(
            f"{i}\t{sid}\t{p}\t{bim[sid]}\t{maxstr}\t{chipstr}\t{status}\t{1 if keep else 0}\t{new_rank}"
        )

    output_lines(out / "top1000_QTN_vs_selected_CHIP_LD.tsv", all_lines)
    output_lines(out / "filtered_ranked_extra_snps.tsv",
                 ["SNP\tP"] + [f"{s}\t{p}" for s, p, *_ in kept])
    output_lines(out / "filtered_ranked_annotated.tsv",
                 ["Filtered_Rank\tOriginal_Rank\tSNP\tP\tMax_r2_with_selected_CHIP\tTop_linked_CHIP"] +
                 [f"{i}\t{old_i}\t{s}\t{p}\t{r}\t{c}" for i, (s, p, old_i, r, c) in enumerate(kept, 1)])
    output_lines(out / "original_top20_LD_diagnostic.tsv", all_lines[:21])

    marker_root = out / "marker_sets"
    output_lines(marker_root / "extra_0.snps", [])
    candidate_status = ["N_QTN\tAvailability\tReason"]
    for n in COUNTS:
        if len(kept) >= n:
            output_lines(marker_root / f"extra_{n}.snps", [s for s, *_ in kept[:n]])
            candidate_status.append(f"{n}\tAVAILABLE\tTop{n}_from_filtered_top1000")
        else:
            # Remove an old diagnostic file if a previously successful partial
            # run is resumed with new inputs. Only our own output directory.
            old = marker_root / f"extra_{n}.snps"
            if old.exists():
                old.unlink()
            candidate_status.append(f"{n}\tUNAVAILABLE\tOnly_{len(kept)}_survivors_in_tested_top1000")
    output_lines(out / "candidate_count_availability.tsv", candidate_status)
    total = len(top)
    summary = [
        ("Trait", "BB"), ("Outer", "01"),
        ("N_train", n_fam), ("N_selected_CHIP", len(chip)),
        ("N_original_QTN_candidates", total),
        ("LD_threshold_remove_r2_ge", threshold),
        ("N_QTN_filtered_out", n_removed),
        ("N_QTN_retained", len(kept)),
        ("N_QTN_without_selected_CHIP_on_same_chr", n_no_chip),
        ("Chromosomes_without_selected_CHIP", ",".join(sorted(no_chip_chrs)) or "none"),
        ("N_original_top20_filtered_out", sum(1 for s, *_ in top[:20] if s not in set(s for s, *_ in kept))),
        ("N_available_positive_candidate_counts", sum(n <= len(kept) for n in COUNTS)),
        ("GWAS_source", args.ranked),
    ]
    output_lines(out / "filter_summary.tsv", ["Item\tValue"] + [f"{key}\t{v}" for key, v in summary])
    print(f"[PASS] {total} QTN screened; r2 >= {threshold}: removed {n_removed}; "
          f"retained {len(kept)}; no same-chromosome selected CHIP: {n_no_chip}")
    print(f"[PASS] Largest available marker set: {max((x for x in COUNTS if x <= len(kept)), default=0)}")
    print(f"[RESULT] {out / 'filter_summary.tsv'}")
    print(f"[RESULT] {out / 'filtered_ranked_extra_snps.tsv'}")
    print(f"[RESULT] {out / 'top1000_QTN_vs_selected_CHIP_LD.tsv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="step", required=True)
    p = subparsers.add_parser("prepare")
    p.add_argument("--ranked", required=True)
    p.add_argument("--chip", required=True)
    p.add_argument("--out", required=True)
    r = subparsers.add_parser("report")
    r.add_argument("--ranked", required=True)
    r.add_argument("--bim", required=True)
    r.add_argument("--fam", required=True)
    r.add_argument("--ld", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--threshold", type=float, default=0.8)
    args = parser.parse_args()
    if args.step == "prepare":
        prep(args)
    else:
        if not 0 < args.threshold <= 1:
            parser.error("LD threshold must be > 0 and <= 1")
        report(args)


if __name__ == "__main__":
    main()
