#!/usr/bin/env python3
"""Nine lead-SNP to candidate-gene-body LD comparisons in Norway spruce.

Only BB and DBH, and only seven trait-gene pairs (nine lead-gene links),
pre-selected using the same-trait original/imputed +/-50-kb GWAS comparison.
Input is a validated functional annotation TSV (S16 seed) and existing GWAS
outputs. No GWAS recalculation and no invented LD results.

PLINK 1.9 full input -> temporary small SNP subset per trait -> pairwise r2.
"""
from __future__ import annotations
import argparse
import csv
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

EXPECTED = {
    ("BB", "PA_chr02_G002738", "2:539992175"),
    ("BB", "PA_chr02_G002738", "2:540097312"),
    ("BB", "PA_chr02_G002738", "2:540184241"),
    ("BB", "PA_chr12_G001028", "12:281962397"),
    ("BB", "PA_chr12_G001029", "12:281962397"),
    ("DBH", "PA_chr12_G001011", "12:279011806"),
    ("DBH", "PA_chr12_G001012", "12:279011806"),
    ("DBH", "PA_chr12_G001028", "12:281941149"),
    ("DBH", "PA_chr12_G001029", "12:281941149"),
}
SUMMARY_COLS = [
    "Trait", "Candidate_gene_ID", "Lead_index_SNP",
    "Gene_body_SNPs_excluding_index", "N_valid_LD_pairs",
    "Maximum_gene_body_r2", "Mean_gene_body_r2", "Top_gene_body_SNP",
    "Top_gene_body_SNP_in_original_BIM", "N_gene_body_SNPs_r2_ge_0.2",
    "N_gene_body_SNPs_r2_ge_0.8", "LD_status",
]
APPENDED = SUMMARY_COLS[3:]


def load_tsv(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def save_tsv(path: Path, cols, rows):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def check_seed(seed: Path, candidates: Path, regions: Path):
    cols, rows = load_tsv(seed)
    required = {"Trait", "Candidate_gene_ID", "Gene_coordinates_v2.0_bp", "Lead_index_SNP", "Lead_SNP_P_value"}
    if required - set(cols):
        raise ValueError(f"Seed missing columns: {sorted(required - set(cols))}")
    keys = [(r["Trait"], r["Candidate_gene_ID"], r["Lead_index_SNP"]) for r in rows]
    if len(rows) != 9 or len(set(keys)) != 9 or set(keys) != EXPECTED:
        raise ValueError(f"Seed must contain the exact 9 confirmed links. Missing={EXPECTED-set(keys)}; unexpected={set(keys)-EXPECTED}")
    for r in rows:
        m = re.fullmatch(r"PA_chr(\d{2}):(\d+)-(\d+)", r["Gene_coordinates_v2.0_bp"])
        if not m:
            raise ValueError(f"Invalid gene coordinates for {r['Candidate_gene_ID']}")
        chrom, start, end = map(int, m.groups())
        if start < 1 or end < start or int(r["Lead_index_SNP"].split(":")[0]) != chrom:
            raise ValueError(f"Invalid chromosome/genomic interval: {r}")
        r["_chr"], r["_start"], r["_end"] = chrom, start, end
        try:
            p = float(r["Lead_SNP_P_value"])
        except ValueError as e:
            raise ValueError(f"Invalid P value for {r}") from e
        if not 0 < p <= 1:
            raise ValueError(f"Out-of-range P value for {r}")
    # Cross-check the analysis source, not just a manually transcribed seed.
    ccols, crows = load_tsv(candidates)
    required_c = {"Trait", "Gene_ID", "Classification", "Gene_start", "Gene_end", "Imputed_lead_SNPs"}
    if required_c - set(ccols):
        raise ValueError(f"Gene comparison source missing columns: {required_c - set(ccols)}")
    reference = {}
    for r in crows:
        if r["Classification"] != "Imputed-specific":
            continue
        key = (r["Trait"], r["Gene_ID"])
        reference[key] = r
    if set(reference) != {(t, g) for t, g, snp in EXPECTED}:
        raise ValueError("The trait-specific candidate gene set does not match the S16 seed")
    for r in rows:
        source = reference[(r["Trait"], r["Candidate_gene_ID"])]
        if r["Lead_index_SNP"] not in source["Imputed_lead_SNPs"].split(";"):
            raise ValueError(f"Lead SNP not matched in GWAS candidate-gene table: {r}")
        if (r["_start"], r["_end"]) != (int(source["Gene_start"]), int(source["Gene_end"])):
            raise ValueError(f"Gene coordinates differ between sources for {r['Candidate_gene_ID']}")
    pcols, prows = load_tsv(regions)
    required_p = {"Trait", "Dataset", "Lead_SNP", "P"}
    if required_p - set(pcols):
        raise ValueError(f"Region manifest missing columns: {required_p - set(pcols)}")
    region_map = {(r["Trait"], r["Lead_SNP"]): float(r["P"]) for r in prows if r["Dataset"] == "Imputed"}
    for r in rows:
        refp = region_map.get((r["Trait"], r["Lead_index_SNP"]))
        p = float(r["Lead_SNP_P_value"])
        if refp is None or not math.isclose(p, refp, rel_tol=0.005):
            raise ValueError(f"GWAS P-value mismatch for {r['Trait']} {r['Lead_index_SNP']}: seed={p}, source={refp}")
    print("Verified: 9 links / 7 gene-trait combinations / 5 unique genes / 6 distinct lead SNPs")
    return cols, rows


def prefix_exists(prefix: Path, label: str, samples: int):
    for ext in (".bed", ".bim", ".fam"):
        p = Path(str(prefix) + ext)
        if not p.is_file():
            raise FileNotFoundError(f"Missing {label} PLINK input: {p}")
    with Path(str(prefix) + ".fam").open() as fh:
        count = sum(1 for _ in fh)
    if count != samples:
        raise ValueError(f"{label} PLINK sample count {count} differs from expected {samples}; check input prefix")
    print(f"Verified: {label} n={samples}, PLINK prefix={prefix}")


def parse_chr(s):
    m = re.fullmatch(r"(?:PA_chr|chr)?0*(\d+)", s, re.IGNORECASE)
    return int(m.group(1)) if m else None


def variants_in_gene_intervals(bim: Path, trait: str, rows: list[dict]):
    intervals = {}
    for r in rows:
        if r["Trait"] == trait:
            intervals[r["Candidate_gene_ID"]] = (r["_chr"], r["_start"], r["_end"])
    by_chr = defaultdict(list)
    for gid, (chrom, s, e) in intervals.items():
        by_chr[chrom].append((gid, s, e))
    wanted_leads = {r["Lead_index_SNP"] for r in rows if r["Trait"] == trait}
    found_leads = set()
    gene_ids = defaultdict(set)
    selected_ids = set()
    with bim.open(encoding="utf-8") as fh:
        for ln, line in enumerate(fh, 1):
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Malformed BIM at {bim}:{ln}")
            ch = parse_chr(parts[0])
            if ch is None or ch not in by_chr:
                continue
            sid, pos_str = parts[1], parts[3]
            try:
                bp = int(pos_str)
            except ValueError:
                raise ValueError(f"Invalid BIM position at {bim}:{ln}")
            if sid in wanted_leads:
                found_leads.add(sid)
                selected_ids.add(sid)
            for gid, start, end in by_chr[ch]:
                if start <= bp <= end:
                    gene_ids[gid].add(sid)
                    selected_ids.add(sid)
    if found_leads != wanted_leads:
        raise ValueError(f"{trait}: required lead SNP(s) missing from imputed BIM: {wanted_leads-found_leads}")
    print(f"{trait}: retained {len(selected_ids)} unique lead/gene-body marker IDs before PLINK extraction")
    for gid in sorted(intervals):
        print(f"    {gid}: {len(gene_ids[gid])} gene-body markers")
    return dict(gene_ids), selected_ids


def plink_run(args, prefix, logs_dir: Path):
    print("$", " ".join(map(str, args)), flush=True)
    result = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"PLINK failed on {prefix}; log:\n{result.stdout[-3500:]}")
    # Standard PLINK log is saved independently as --out prefix.log.
    return result.stdout


def file_pairwise_r2(path: Path, lead: str, valid_gene_ids: set[str]):
    if not path.is_file():
        raise FileNotFoundError(f"PLINK LD file missing: {path}")
    found = {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().split()
        positions = {h.upper(): i for i, h in enumerate(header)}
        if not {"SNP_A", "SNP_B", "R2"}.issubset(positions):
            raise ValueError(f"Unexpected PLINK 1.9 --r2 header in {path}: {header}")
        for ln, line in enumerate(fh, 2):
            tok = line.split()
            if len(tok) < len(header):
                raise ValueError(f"Bad LD line at {path}:{ln}")
            a, b = tok[positions["SNP_A"]], tok[positions["SNP_B"]]
            target = b if a == lead else (a if b == lead else None)
            if not target or target == lead or target not in valid_gene_ids:
                continue
            value = float(tok[positions["R2"]])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"Invalid r2 in {path}:{ln}: {value}")
            if target in found:
                raise ValueError(f"Duplicate lead-gene SNP pair in {path}: {target}")
            found[target] = value
    return found


def original_lookup(bim_path: str, tracked_ids: set[str]):
    if not bim_path:
        return None
    bim = Path(bim_path)
    if not bim.is_file():
        raise FileNotFoundError(f"Original BIM missing: {bim}")
    detected = set()
    with bim.open(encoding="utf-8") as fh:
        for line in fh:
            tok = line.split()
            if len(tok) >= 2 and tok[1] in tracked_ids:
                detected.add(tok[1])
    return detected


def calc_one(args, reduced_prefix: Path, row: dict, gene_ids: dict, out: Path, original_ids):
    trait, gid, lead = row["Trait"], row["Candidate_gene_ID"], row["Lead_index_SNP"]
    eligible = gene_ids.get(gid, set()) - {lead}
    tag = f"{trait}_{gid}_{lead.replace(':','_')}"
    result = {"Trait": trait, "Candidate_gene_ID": gid, "Lead_index_SNP": lead,
              "Gene_body_SNPs_excluding_index": len(eligible)}
    blank = {"N_valid_LD_pairs": 0, "Maximum_gene_body_r2": "NA", "Mean_gene_body_r2": "NA",
             "Top_gene_body_SNP": "NA", "Top_gene_body_SNP_in_original_BIM": "Not_assessed" if original_ids is None else "NA",
             "N_gene_body_SNPs_r2_ge_0.2": 0, "N_gene_body_SNPs_r2_ge_0.8": 0}
    if not eligible:
        result.update({**blank, "LD_status": "NO_GENE_BODY_SNP"})
        return result, []
    extract = out / f"{tag}.extract.txt"
    extract.write_text("\n".join(sorted(eligible | {lead})) + "\n", encoding="utf-8")
    output_prefix = out / tag
    cmd = [args.plink, "--bfile", str(reduced_prefix), "--extract", str(extract),
           "--ld-snp", lead, "--r2", "--ld-window", "999999", "--ld-window-kb", "1000",
           "--ld-window-r2", "0", "--threads", str(args.threads),
           "--memory", str(args.memory), "--out", str(output_prefix)]
    plink_run(cmd, output_prefix, out)
    pairs = file_pairwise_r2(Path(str(output_prefix) + ".ld"), lead, eligible)
    if not pairs:
        result.update({**blank, "LD_status": "NO_VALID_LD_PAIR"})
        return result, []
    sorted_pairs = sorted(pairs.items(), key=lambda x: (-x[1], x[0]))
    top_snp, top_val = sorted_pairs[0]
    result.update({
        "N_valid_LD_pairs": len(pairs),
        "Maximum_gene_body_r2": f"{top_val:.6f}",
        "Mean_gene_body_r2": f"{sum(pairs.values())/len(pairs):.6f}",
        "Top_gene_body_SNP": top_snp,
        "Top_gene_body_SNP_in_original_BIM": "Not_assessed" if original_ids is None else ("Yes" if top_snp in original_ids else "No"),
        "N_gene_body_SNPs_r2_ge_0.2": sum(v >= 0.2 for v in pairs.values()),
        "N_gene_body_SNPs_r2_ge_0.8": sum(v >= 0.8 for v in pairs.values()),
        "LD_status": "OK",
    })
    details = [
        {"Trait": trait, "Candidate_gene_ID": gid, "Lead_index_SNP": lead,
         "Gene_body_SNP": sid, "r2": f"{r2:.6f}",
         "Gene_body_SNP_in_original_BIM": "Not_assessed" if original_ids is None else ("Yes" if sid in original_ids else "No")}
        for sid, r2 in sorted_pairs
    ]
    print(f"  LD result: {trait}, {gid}, {lead}: max r2={top_val:.6f} (N={len(pairs)})", flush=True)
    return result, details


def report_insertions(path: Path, records):
    d = {(r["Trait"], r["Candidate_gene_ID"], r["Lead_index_SNP"]): r for r in records}
    def fmt(trait, gene, lead):
        item = d[(trait, gene, lead)]
        return f"{float(item['Maximum_gene_body_r2']):.3f}" if item["LD_status"] == "OK" else "not estimable"
    bb = fmt("BB", "PA_chr02_G002738", "2:540097312")
    dbh_a = fmt("DBH", "PA_chr12_G001011", "12:279011806")
    dbh_b = fmt("DBH", "PA_chr12_G001012", "12:279011806")
    methods = (
        "For additional candidate genes identified from the imputed GWAS, pairwise linkage "
        "disequilibrium (LD, r²) was calculated using PLINK v1.9 between each significant "
        "index SNP and variants within the corresponding gene body in the trait-specific "
        "imputed hard-call genotype data. Index-SNP self-comparisons were excluded."
    )
    results = (
        "For the representative BB association signal on chromosome 2, the index SNP "
        f"2:540097312 showed a maximum pairwise r² of {bb} with variants within PA_chr02_G002738. "
        "For DBH, the index SNP 12:279011806 showed maximum r² values of "
        f"{dbh_a} and {dbh_b} with SNPs within PA_chr12_G001011 and PA_chr12_G001012, "
        "respectively (Supplementary Table S16)."
    )
    discussion = (
        "The local LD analysis provided additional genomic context for selected newly detected "
        "association signals and nearby annotated candidate genes, without establishing their "
        "causal roles. In the present study, the main contribution of genotype imputation to GWAS "
        "appears to lie in increased marker coverage and the detection of additional association "
        "signals and trait-specific candidate genes."
    )
    caveat = ("NOTE: Verify all observed LD values in candidate_gene_LD_summary.tsv before inserting Results. "
              "If any value is 'not estimable', revise or omit the affected clause. "
              "The mean and maximum r² depend on tested SNP sets and should not be interpreted as gene-level association evidence.")
    path.write_text("METHODS\n" + methods + "\n\nRESULTS (two illustrative regions only)\n" + results +
                    "\n\nDISCUSSION (integrate; do not append repeated conclusion)\n" + discussion +
                    "\n\n" + caveat + "\n", encoding="utf-8")


def main():
    a = argparse.ArgumentParser()
    for key in ("seed", "candidates", "regions", "bb_bfile", "dbh_bfile", "outdir"):
        a.add_argument("--" + key.replace("_", "-"), required=True)
    a.add_argument("--plink", default="plink")
    a.add_argument("--threads", type=int, default=4)
    a.add_argument("--memory", type=int, default=16000, help="PLINK memory allocation in MB")
    a.add_argument("--bb-original-bim", default="", help="Optional original (non-imputed) BB BIM, for marker-presence flag")
    a.add_argument("--dbh-original-bim", default="", help="Optional original (non-imputed) DBH BIM, for marker-presence flag")
    a.add_argument("--validate-only", action="store_true", help="Validate inputs and read BIM, no LD computation")
    args = a.parse_args()
    if args.threads < 1 or args.memory < 1000:
        raise ValueError("Require --threads >=1 and --memory >=1000 MB")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    basecols, rows = check_seed(Path(args.seed), Path(args.candidates), Path(args.regions))
    prefixes = {"BB": Path(args.bb_bfile), "DBH": Path(args.dbh_bfile)}
    for trait, n in (("BB", 4138), ("DBH", 2674)):
        prefix_exists(prefixes[trait], trait, n)
    if not args.validate_only and not (shutil.which(args.plink) or Path(args.plink).is_file()):
        raise FileNotFoundError(f"PLINK not available: {args.plink}")
    variants = {}
    tracked = {}
    for trait in ("BB", "DBH"):
        variants[trait], tracked[trait] = variants_in_gene_intervals(
            Path(str(prefixes[trait]) + ".bim"), trait, rows)
    original = {
        "BB": original_lookup(args.bb_original_bim, tracked["BB"]),
        "DBH": original_lookup(args.dbh_original_bim, tracked["DBH"]),
    }
    if args.validate_only:
        print("Validation finished; no GWAS or LD jobs were run")
        return
    # Convert each full imputed genotype to a small, local SNP subset ONCE.
    # This avoids repeatedly reading multi-million-variant BEDs for every gene.
    reduced = {}
    for trait in ("BB", "DBH"):
        extract = out / f"{trait}_all_target_variants.extract.txt"
        extract.write_text("\n".join(sorted(tracked[trait])) + "\n", encoding="utf-8")
        small = out / f"{trait}_LD_target_markers"
        cmd = [args.plink, "--bfile", str(prefixes[trait]), "--extract", str(extract),
               "--make-bed", "--keep-allele-order", "--threads", str(args.threads),
               "--memory", str(args.memory), "--out", str(small)]
        plink_run(cmd, small, out)
        for ext in (".bed", ".bim", ".fam"):
            if not Path(str(small) + ext).is_file():
                raise RuntimeError(f"PLINK extraction did not produce {small}{ext}")
        with Path(str(small)+".bim").open() as h:
            observed = {x.split()[1] for x in h if len(x.split()) >= 2}
        if not tracked[trait].issubset(observed):
            raise ValueError(f"{trait}: variants missing after PLINK extraction: {sorted(tracked[trait]-observed)[:10]}")
        reduced[trait] = small
    results, pairs = [], []
    for r in rows:
        trait = r["Trait"]
        summary, details = calc_one(args, reduced[trait], r, variants[trait], out, original[trait])
        results.append(summary)
        pairs.extend(details)
    if len(results) != 9 or set((r["Trait"], r["Candidate_gene_ID"], r["Lead_index_SNP"]) for r in results) != EXPECTED:
        raise RuntimeError("Output row validation failed; refusing to make final supplementary table")
    save_tsv(out / "candidate_gene_LD_summary.tsv", SUMMARY_COLS, results)
    save_tsv(out / "candidate_gene_LD_pairs.tsv",
             ["Trait", "Candidate_gene_ID", "Lead_index_SNP", "Gene_body_SNP", "r2", "Gene_body_SNP_in_original_BIM"], pairs)
    lmap = {(r["Trait"], r["Candidate_gene_ID"], r["Lead_index_SNP"]): r for r in results}
    combined = []
    for r in rows:
        item = {k: r[k] for k in basecols}
        value = lmap[(r["Trait"], r["Candidate_gene_ID"], r["Lead_index_SNP"])]
        item.update({k: value[k] for k in APPENDED})
        combined.append(item)
    save_tsv(out / "Supplementary_Table_S16_with_LD.tsv", basecols + APPENDED, combined)
    report_insertions(out / "Manuscript_LD_insertions.txt", results)
    print("\nFinished:")
    for r in results:
        print(f"  {r['Trait']:3} {r['Candidate_gene_ID']:17} {r['Lead_index_SNP']:15} r2={r['Maximum_gene_body_r2']} status={r['LD_status']}")
    for name in ("candidate_gene_LD_summary.tsv", "candidate_gene_LD_pairs.tsv", "Supplementary_Table_S16_with_LD.tsv", "Manuscript_LD_insertions.txt"):
        print(" ", out / name)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
