#!/bin/bash
#SBATCH -J 11a_BB01_chip_qtn_LD
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -o logs/11a_BB01_chip_qtn_LD_%j.out
#SBATCH -e logs/11a_BB01_chip_qtn_LD_%j.err

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR at line ${LINENO}: ${BASH_COMMAND} (exit ${rc})" >&2' ERR
export LC_ALL=C

module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
TRAIN="${WORK}/10_outer_test/02_imputed_inputs/outer_01/outer_train.keep"
CHIP="${WORK}/10_outer_test/01_chip_selection/outer_01/selected/BB_best_chip.snps"
QTN="${WORK}/10_outer_test/04_imputed_clumping_expanded_qtn/outer_01/BB/marker_sets/extra_20.snps"
OUT="${WORK}/10_outer_test/06_COJO_pilot/outer_01/BB"
mkdir -p "$OUT"

for f in "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam" "$TRAIN" "$CHIP" "$QTN"; do
    if [[ ! -s "$f" ]]; then echo "ERROR: missing/empty input: $f" >&2; exit 1; fi
done

[[ $(wc -l < "$CHIP") -eq 100 ]] || { echo "ERROR: expected 100 CHIP SNPs" >&2; exit 1; }
[[ $(wc -l < "$QTN") -eq 20 ]] || { echo "ERROR: expected 20 QTN SNPs" >&2; exit 1; }

# Check original CHIP and expanded-QTN sets are disjoint; do not overwrite either.
if [[ -n $(comm -12 <(sort -u "$CHIP") <(sort -u "$QTN")) ]]; then
    echo "ERROR: selected CHIP and QTN lists overlap" >&2
    exit 1
fi

cat "$CHIP" "$QTN" | LC_ALL=C sort -u > "$OUT/selected_120.snps"
[[ $(wc -l < "$OUT/selected_120.snps") -eq 120 ]] || { echo "ERROR: expected 120 distinct SNPs" >&2; exit 1; }

# Same outer-training individuals and variant-level QC as Step 10b1.
plink \
  --bfile "$MASTER" \
  --keep "$TRAIN" \
  --extract "$OUT/selected_120.snps" \
  --maf 0.01 --geno 0.05 \
  --keep-allele-order \
  --make-bed \
  --threads 2 --memory 6000 \
  --out "$OUT/selected_120_train_qc"

# Any marker excluded by QC would make the full 100+20 conditional analysis invalid.
awk '
    NR==FNR {missing[$1]=1; next}
    { delete missing[$2] }
    END {for (s in missing) print s}
' "$OUT/selected_120.snps" "$OUT/selected_120_train_qc.bim" \
    | LC_ALL=C sort > "$OUT/missing_after_qc.snps"

if [[ -s "$OUT/missing_after_qc.snps" ]]; then
    echo "ERROR: some selected CHIP/QTN SNPs do not pass training-fold QC:" >&2
    cat "$OUT/missing_after_qc.snps" >&2
    echo "STOP: inspect before constructing a full conditional model." >&2
    exit 1
fi

[[ $(wc -l < "$OUT/selected_120_train_qc.bim") -eq 120 ]] || {
  echo "ERROR: expected exactly 120 retained SNPs" >&2; exit 1;
}
[[ $(wc -l < "$OUT/selected_120_train_qc.fam") -eq 813 ]] || {
  echo "ERROR: expected 813 outer-training individuals" >&2; exit 1;
}

# Compute r^2 for all marker pairs on the same chromosome, irrespective of distance.
plink \
  --bfile "$OUT/selected_120_train_qc" \
  --r2 \
  --ld-window 999999 \
  --ld-window-kb 10000000 \
  --ld-window-r2 0 \
  --threads 2 --memory 6000 \
  --out "$OUT/selected_120_train_LD"

[[ -s "$OUT/selected_120_train_LD.ld" ]] || {
  echo "ERROR: missing LD output" >&2; exit 1;
}

# PLINK 1.9 --r2 .ld columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
# Report max within-chromosome r^2 for each QTN against the 100 selected CHIP SNPs.
awk 'BEGIN {OFS="\t"}
  FILENAME == ARGV[1] { qtn[$1]=1; order[++n]=$1; next }
  FILENAME == ARGV[2] { chip[$1]=1; next }
  FNR == 1 {
    if ($3!="SNP_A" || $6!="SNP_B" || $7!="R2") {
      print "ERROR: unexpected PLINK LD header" > "/dev/stderr"
      bad=1
    }
    next
  }
  {
    a=$3; b=$6; r=$7+0
    if ((a in chip) && (b in qtn)) {q=b; c=a}
    else if ((b in chip) && (a in qtn)) {q=a; c=b}
    else next
    if (!(q in max_r2) || r>max_r2[q]) {
      max_r2[q]=r; closest_chip[q]=c
    }
  }
  END {
    if (bad) exit 1
    print "QTN", "Max_r2_with_selected_CHIP", "Top_linked_CHIP"
    for (i=1; i<=n; i++) {
      q=order[i]
      if (q in max_r2) printf "%s\t%.7f\t%s\n",q,max_r2[q],closest_chip[q]
      else print q, "NA", "NA"
    }
  }
' "$QTN" "$CHIP" "$OUT/selected_120_train_LD.ld" \
    > "$OUT/BB_outer01_QTN_vs_CHIP_max_r2.tsv"

[[ $(wc -l < "$OUT/BB_outer01_QTN_vs_CHIP_max_r2.tsv") -eq 21 ]] || {
  echo "ERROR: expected 21 summary lines" >&2; exit 1;
}

printf '\n===== BB outer_01: QTN vs best CHIP LD =====\n'
column -t "$OUT/BB_outer01_QTN_vs_CHIP_max_r2.tsv" || cat "$OUT/BB_outer01_QTN_vs_CHIP_max_r2.tsv"
echo "DONE. Output: $OUT/BB_outer01_QTN_vs_CHIP_max_r2.tsv"
