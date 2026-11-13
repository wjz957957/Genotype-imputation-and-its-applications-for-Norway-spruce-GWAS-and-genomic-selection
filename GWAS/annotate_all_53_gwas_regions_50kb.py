#!/usr/bin/env python3
"""Fresh 50-kb gene annotation for all Original and Imputed GWAS index SNPs."""
import argparse, csv, gzip, re, sys
from pathlib import Path
from collections import defaultdict
from bisect import bisect_left, bisect_right

EXPECTED = {
    ("BB","Original"):6, ("BB","Imputed"):23,
    ("DBH","Original"):1, ("DBH","Imputed"):3,
    ("FD","Original"):7, ("FD","Imputed"):13,
}
REGION_REQUIRED = ["Region_ID","Trait","Dataset","CHR","BP","Lead_SNP","P"]

def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def gene_id_from_attr(text):
    m = re.search(r"(?:^|;)ID=([^;]+)", text.strip().strip(";"))
    return m.group(1) if m else None

def read_gff(path):
    by_chr, seen = defaultdict(list), {}
    total = skipped = no_id = 0
    with opener(path) as f:
        for ln, line in enumerate(f, 1):
            if not line or line.startswith("#"): continue
            x = line.rstrip("\n").split("\t")
            if len(x) != 9 or x[2] != "gene": continue
            total += 1
            m = re.fullmatch(r"PA_chr(\d{2})", x[0])
            if not m:
                skipped += 1
                continue
            ch = int(m.group(1))
            if not 1 <= ch <= 12:
                skipped += 1
                continue
            gid = gene_id_from_attr(x[8])
            if not gid:
                no_id += 1
                continue
            try:
                start, end = int(x[3]), int(x[4])
            except ValueError:
                raise ValueError(f"Invalid coordinates: {path}:{ln}")
            if start < 1 or end < start:
                raise ValueError(f"Invalid interval: {path}:{ln}")
            sig = (start, end, x[6], x[1])
            key = (ch, gid)
            if key in seen:
                if seen[key] != sig:
                    raise ValueError(f"Conflicting gene records for {gid}")
                continue
            seen[key] = sig
            by_chr[ch].append({
                "Gene_ID":gid, "Chromosome":x[0],
                "Gene_start":start, "Gene_end":end,
                "Gene_strand":x[6], "Gene_length_bp":end-start+1,
                "Gene_source":x[1], "Gene_attributes":x[8]
            })
    for ch in by_chr:
        by_chr[ch].sort(key=lambda g:(g["Gene_start"],g["Gene_end"],g["Gene_ID"]))
    print(f"GFF3 gene features={total}; canonical chr genes={sum(map(len,by_chr.values()))}; noncanonical skipped={skipped}; missing IDs={no_id}")
    return by_chr

def read_longest_ids(path):
    ids=set()
    with opener(path) as f:
        for line in f:
            if not line or line.startswith("#"): continue
            x=line.rstrip("\n").split("\t")
            if len(x)==9 and x[2]=="gene":
                gid=gene_id_from_attr(x[8])
                if gid: ids.add(gid)
    print(f"longest_no_TE gene IDs={len(ids)}")
    return ids

def read_regions(path, skip_count_check=False):
    with open(path,newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        missing=[c for c in REGION_REQUIRED if c not in (reader.fieldnames or [])]
        if missing: raise ValueError("Missing region columns: "+",".join(missing))
        rows=[]; seen=set(); counts=defaultdict(int)
        for ln,r in enumerate(reader,2):
            rid=r["Region_ID"].strip()
            if not rid or rid in seen: raise ValueError(f"Empty/duplicate Region_ID at line {ln}: {rid}")
            seen.add(rid)
            trait=r["Trait"].strip(); ds=r["Dataset"].strip()
            if trait not in {"BB","DBH","FD"} or ds not in {"Original","Imputed"}:
                raise ValueError(f"Unexpected trait/dataset at line {ln}")
            try: ch,bp=int(r["CHR"]),int(r["BP"])
            except ValueError: raise ValueError(f"Invalid CHR/BP at line {ln}")
            if ch not in range(1,13) or bp<1: raise ValueError(f"Out-of-range CHR/BP at line {ln}")
            if r["Lead_SNP"].strip()!=f"{ch}:{bp}":
                raise ValueError(f"Lead SNP coordinate mismatch at line {ln}: {r['Lead_SNP']}")
            r["_chr"]=ch; r["_bp"]=bp; r["_chrom"]=f"PA_chr{ch:02d}"
            rows.append(r); counts[(trait,ds)]+=1
    if not skip_count_check:
        if len(rows)!=53: raise ValueError(f"Expected 53 index SNP rows, found {len(rows)}")
        for k,v in EXPECTED.items():
            if counts[k]!=v: raise ValueError(f"Expected {v} rows for {k}, found {counts[k]}")
    print("Region counts:")
    for t in ("BB","DBH","FD"):
        print(f"  {t}: Original={counts[(t,'Original')]} Imputed={counts[(t,'Imputed')]}")
    print(f"  Total={len(rows)}")
    return rows

def build_index(by_chr):
    idx={}
    for ch,genes in by_chr.items():
        starts=[g["Gene_start"] for g in genes]
        prefix=[]; mx=0
        for g in genes:
            mx=max(mx,g["Gene_end"]); prefix.append(mx)
        idx[ch]=(genes,starts,prefix)
    return idx

def overlap_hits(ch,lo,hi,idx):
    if ch not in idx: return []
    genes,starts,prefix=idx[ch]
    left=bisect_left(prefix,lo)
    right=bisect_right(starts,hi)
    return [g for g in genes[left:right] if g["Gene_end"]>=lo]

def dist(bp,g):
    if bp<g["Gene_start"]: return g["Gene_start"]-bp
    if bp>g["Gene_end"]: return bp-g["Gene_end"]
    return 0

def rel(bp,g):
    if g["Gene_start"]<=bp<=g["Gene_end"]: return "OVERLAPS_LEAD_SNP"
    return "GENOMIC_LEFT_OF_LEAD" if g["Gene_end"]<bp else "GENOMIC_RIGHT_OF_LEAD"

def write_tsv(path,rows,cols):
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore",lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def aggregate(candidate, dataset):
    grouped=defaultdict(list)
    for r in candidate:
        if r["Dataset"]==dataset: grouped[r["Gene_ID"]].append(r)
    out=[]
    for gid,rs in sorted(grouped.items()):
        g=rs[0]
        out.append({**{k:g[k] for k in GENE_FIELDS},
            "Traits":";".join(sorted({r["Trait"] for r in rs})),
            "N_association_region_records":len(rs),
            "Associated_region_IDs":";".join(sorted({r["Region_ID"] for r in rs})),
            "Associated_lead_SNPs":";".join(sorted({r["Lead_SNP"] for r in rs})),
            "Minimum_distance_to_lead_bp":min(int(r["Distance_to_lead_bp"]) for r in rs)})
    return out

GENE_FIELDS=["Gene_ID","Chromosome","Gene_start","Gene_end","Gene_strand","Gene_length_bp","Gene_source","Present_in_longest_no_TE","Gene_attributes"]

def compare_genes(candidate):
    td=defaultdict(list); gd=defaultdict(list)
    for r in candidate:
        td[(r["Trait"],r["Dataset"],r["Gene_ID"])].append(r)
        gd[(r["Dataset"],r["Gene_ID"])].append(r)
    trait_rows=[]; trait_imp=[]
    for trait in ("BB","DBH","FD"):
        gids=sorted({gid for t,ds,gid in td if t==trait})
        for gid in gids:
            o=td.get((trait,"Original",gid),[]); i=td.get((trait,"Imputed",gid),[])
            g=(i or o)[0]
            cl="Shared" if o and i else ("Original-specific" if o else "Imputed-specific")
            row={"Trait":trait,"Gene_ID":gid,**{k:g[k] for k in GENE_FIELDS if k!="Gene_attributes"},
                 "Classification":cl,"Original_region_count":len(o),"Imputed_region_count":len(i),
                 "Original_lead_SNPs":";".join(sorted({r["Lead_SNP"] for r in o})),
                 "Imputed_lead_SNPs":";".join(sorted({r["Lead_SNP"] for r in i})),
                 "Original_region_IDs":";".join(sorted({r["Region_ID"] for r in o})),
                 "Imputed_region_IDs":";".join(sorted({r["Region_ID"] for r in i}))}
            trait_rows.append(row)
            if cl=="Imputed-specific": trait_imp.append(row)
    global_rows=[]; global_imp=[]
    gids=sorted({gid for ds,gid in gd})
    for gid in gids:
        o=gd.get(("Original",gid),[]); i=gd.get(("Imputed",gid),[])
        g=(i or o)[0]
        cl="Shared" if o and i else ("Original-specific" if o else "Imputed-specific")
        row={"Gene_ID":gid,**{k:g[k] for k in GENE_FIELDS if k!="Gene_attributes"},
             "Classification":cl,
             "Original_traits":";".join(sorted({r["Trait"] for r in o})),
             "Imputed_traits":";".join(sorted({r["Trait"] for r in i})),
             "Original_region_count":len(o),"Imputed_region_count":len(i),
             "Original_lead_SNPs":";".join(sorted({r["Lead_SNP"] for r in o})),
             "Imputed_lead_SNPs":";".join(sorted({r["Lead_SNP"] for r in i}))}
        global_rows.append(row)
        if cl=="Imputed-specific": global_imp.append(row)
    return trait_rows,trait_imp,global_rows,global_imp

def main():
    p=argparse.ArgumentParser(description="Fresh ±50-kb annotation of 53 GWAS index SNPs.")
    p.add_argument("--regions",required=True)
    p.add_argument("--gff",required=True)
    p.add_argument("--longest",required=True)
    p.add_argument("--outdir",required=True)
    p.add_argument("--window",type=int,default=50000)
    p.add_argument("--skip-count-check",action="store_true")
    a=p.parse_args()
    for fn in (a.regions,a.gff,a.longest):
        if not Path(fn).is_file(): raise FileNotFoundError(fn)
    if a.window<0: raise ValueError("--window must be non-negative")
    regions=read_regions(a.regions,a.skip_count_check)
    by_chr=read_gff(a.gff)
    longest=read_longest_ids(a.longest)
    missing=sorted({r["_chr"] for r in regions if r["_chr"] not in by_chr})
    if missing: raise ValueError("No canonical GFF genes on chromosome(s): "+",".join(map(str,missing)))
    idx=build_index(by_chr)
    candidates=[]; nearest_rows=[]; summaries=[]
    for r in regions:
        ch,bp=r["_chr"],r["_bp"]; chrom=r["_chrom"]
        lo=max(1,bp-a.window); hi=bp+a.window
        hits=overlap_hits(ch,lo,hi,idx)
        nearest=min(by_chr[ch],key=lambda g:(dist(bp,g),g["Gene_start"],g["Gene_ID"]))
        for g in hits:
            candidates.append({
                "Region_ID":r["Region_ID"],"Trait":r["Trait"],"Dataset":r["Dataset"],
                "Lead_SNP":r["Lead_SNP"],"Lead_P":r["P"],"Chromosome":chrom,
                "Lead_position":bp,"Window_start":lo,"Window_end":hi,
                **g,"Present_in_longest_no_TE":"YES" if g["Gene_ID"] in longest else "NO",
                "Distance_to_lead_bp":dist(bp,g),"Lead_gene_relationship":rel(bp,g)})
        nd=dist(bp,nearest); nin=nearest in hits
        near={"Nearest_gene_ID":nearest["Gene_ID"],"Nearest_gene_start":nearest["Gene_start"],
              "Nearest_gene_end":nearest["Gene_end"],"Nearest_gene_distance_bp":nd,
              "Nearest_gene_in_50kb":"YES" if nin else "NO",
              "Nearest_gene_source":nearest["Gene_source"],
              "Nearest_gene_present_in_longest_no_TE":"YES" if nearest["Gene_ID"] in longest else "NO"}
        summaries.append({"Region_ID":r["Region_ID"],"Trait":r["Trait"],"Dataset":r["Dataset"],
            "Lead_SNP":r["Lead_SNP"],"Lead_P":r["P"],"Chromosome":chrom,"Lead_position":bp,
            "Window_start":lo,"Window_end":hi,"Search_window_bp":a.window,
            "N_candidate_genes_within_window":len(hits),
            "N_candidate_genes_in_longest_no_TE":sum(g["Gene_ID"] in longest for g in hits),**near})
        nearest_rows.append({"Region_ID":r["Region_ID"],"Trait":r["Trait"],"Dataset":r["Dataset"],
            "Lead_SNP":r["Lead_SNP"],"Lead_P":r["P"],"Chromosome":chrom,"Lead_position":bp,
            "Search_window_bp":a.window,**near,
            "Nearest_gene_strand":nearest["Gene_strand"],
            "Nearest_gene_length_bp":nearest["Gene_length_bp"]})
        print(f"{r['Trait']:3s} {r['Dataset']:8s} {r['Lead_SNP']:>14s}: "
              f"{len(hits)} candidate gene(s); nearest={nearest['Gene_ID']} ({nd} bp)")
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    cand_cols=["Region_ID","Trait","Dataset","Lead_SNP","Lead_P","Chromosome","Lead_position",
        "Window_start","Window_end"]+GENE_FIELDS+["Distance_to_lead_bp","Lead_gene_relationship"]
    write_tsv(out/"region_gene_candidates_50kb.tsv",candidates,cand_cols)
    summary_cols=["Region_ID","Trait","Dataset","Lead_SNP","Lead_P","Chromosome","Lead_position",
        "Window_start","Window_end","Search_window_bp","N_candidate_genes_within_window",
        "N_candidate_genes_in_longest_no_TE","Nearest_gene_ID","Nearest_gene_start","Nearest_gene_end",
        "Nearest_gene_distance_bp","Nearest_gene_in_50kb","Nearest_gene_source",
        "Nearest_gene_present_in_longest_no_TE"]
    write_tsv(out/"region_annotation_summary.tsv",summaries,summary_cols)
    near_cols=["Region_ID","Trait","Dataset","Lead_SNP","Lead_P","Chromosome","Lead_position",
        "Search_window_bp","Nearest_gene_ID","Nearest_gene_start","Nearest_gene_end",
        "Nearest_gene_distance_bp","Nearest_gene_in_50kb","Nearest_gene_source",
        "Nearest_gene_present_in_longest_no_TE","Nearest_gene_strand","Nearest_gene_length_bp"]
    write_tsv(out/"region_nearest_genes.tsv",nearest_rows,near_cols)
    uniq_cols=GENE_FIELDS+["Traits","N_association_region_records","Associated_region_IDs",
        "Associated_lead_SNPs","Minimum_distance_to_lead_bp"]
    orig=aggregate(candidates,"Original"); imp=aggregate(candidates,"Imputed")
    write_tsv(out/"unique_genes_original.tsv",orig,uniq_cols)
    write_tsv(out/"unique_genes_imputed.tsv",imp,uniq_cols)
    tr,ti,gl,gi=compare_genes(candidates)
    tr_cols=["Trait","Gene_ID","Chromosome","Gene_start","Gene_end","Gene_strand","Gene_length_bp",
        "Gene_source","Present_in_longest_no_TE","Classification","Original_region_count",
        "Imputed_region_count","Original_lead_SNPs","Imputed_lead_SNPs",
        "Original_region_IDs","Imputed_region_IDs"]
    write_tsv(out/"gene_comparison_by_trait.tsv",tr,tr_cols)
    write_tsv(out/"imputed_specific_genes.tsv",ti,tr_cols)
    gl_cols=["Gene_ID","Chromosome","Gene_start","Gene_end","Gene_strand","Gene_length_bp",
        "Gene_source","Present_in_longest_no_TE","Classification","Original_traits","Imputed_traits",
        "Original_region_count","Imputed_region_count","Original_lead_SNPs","Imputed_lead_SNPs"]
    write_tsv(out/"gene_comparison_global.tsv",gl,gl_cols)
    write_tsv(out/"imputed_specific_genes_global.tsv",gi,gl_cols)
    count_rows=[]
    for scope,rows in [("BB",[x for x in tr if x["Trait"]=="BB"]),
                       ("DBH",[x for x in tr if x["Trait"]=="DBH"]),
                       ("FD",[x for x in tr if x["Trait"]=="FD"]),
                       ("GLOBAL",gl)]:
        count_rows.append({"Scope":scope,
            "Original_unique_genes":sum(x["Classification"] in ("Shared","Original-specific") for x in rows),
            "Imputed_unique_genes":sum(x["Classification"] in ("Shared","Imputed-specific") for x in rows),
            "Shared_genes":sum(x["Classification"]=="Shared" for x in rows),
            "Original_specific_genes":sum(x["Classification"]=="Original-specific" for x in rows),
            "Imputed_specific_genes":sum(x["Classification"]=="Imputed-specific" for x in rows)})
    write_tsv(out/"gene_comparison_counts.tsv",count_rows,
        ["Scope","Original_unique_genes","Imputed_unique_genes","Shared_genes",
         "Original_specific_genes","Imputed_specific_genes"])
    with open(out/"README_analysis.txt","w") as f:
        f.write(f"Regions input: {a.regions}\nFull GFF3: {a.gff}\nlongest_no_TE GFF3: {a.longest}\n")
        f.write(f"Window: +/-{a.window} bp\nCandidate definition: gene interval overlaps the window.\n")
        f.write("Full GFF3 canonical PA_chr01-PA_chr12 genes are included. longest_no_TE is a flag only, not a filter.\n")
        f.write("Nearest genes outside the window are supplementary and excluded from candidate gene sets.\n")
        f.write("Trait-specific comparison is primary; global comparison is supplementary.\n")
    print("\nCompleted.")
    print(f"Regions={len(regions)}; region-gene records={len(candidates)}")
    print(f"Original unique genes={len(orig)}; Imputed unique genes={len(imp)}")
    print(f"Trait-specific Imputed-specific records={len(ti)}; globally Imputed-specific genes={len(gi)}")
    print(f"Output: {out}")

if __name__=="__main__":
    try: main()
    except Exception as e:
        print(f"ERROR: {e}",file=sys.stderr); sys.exit(1)
