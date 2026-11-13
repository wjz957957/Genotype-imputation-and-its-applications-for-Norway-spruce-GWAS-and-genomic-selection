#!/bin/bash
#SBATCH -J 11c_BB01_cojoLDpruned
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH -o logs/11c_BB01_COJO_%j.out
#SBATCH -e logs/11c_BB01_COJO_%j.err

# This is a diagnostic conditional analysis, NOT conditioning on all 100
# selected CHIP SNPs: high-LD CHIP predictors must first be de-redundified.
# Does NOT redo GCTA MLMA or modify the existing 11b output files.
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
export LC_ALL=C

module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
OUT="${WORK}/10_outer_test/06_COJO_pilot/outer_01/BB"
LD_BFILE="${OUT}/selected_120_train_qc"
MA="${OUT}/BB_outer01_genomewide_with_matched_chip.ma"

for f in "${LD_BFILE}.bed" "${LD_BFILE}.bim" "${LD_BFILE}.fam" "$MA" \
         "${OUT}/chr1_conditional_chip.snps" "${OUT}/chr12_conditional_chip.snps" \
         "${OUT}/chr1_selected.snps" "${OUT}/chr12_selected.snps"; do
    [[ -s "$f" ]] || { echo "Missing required input: $f" >&2; exit 1; }
done

# Across each whole chromosome, PLINK checks every selected CHIP SNP against
# all other selected CHIP SNPs (window 1000 SNPs > 73 or 27 SNPs).
# Threshold 0.2 is conservative: high pairwise LD is pruned before COJO.
# Note: pairwise pruning does not mathematically guarantee low multiple R^2.
for CHR in 1 12; do
    COND="${OUT}/chr${CHR}_conditional_chip.snps"
    PRUNE="${OUT}/chr${CHR}_cond_chip_LDprune_r2_0.2"
    echo "========== CHR${CHR}: LD-prune the selected CHIP conditions =========="

    plink \
      --bfile "$LD_BFILE" \
      --chr "$CHR" \
      --extract "$COND" \
      --indep-pairwise 1000 1 0.2 \
      --threads 4 --memory 9000 \
      --out "$PRUNE"

    [[ -s "${PRUNE}.prune.in" ]] || { echo "No conditional CHIP SNP retained: chr${CHR}" >&2; exit 1; }
    [[ -f "${PRUNE}.prune.out" ]] || { echo "Missing pruning report: chr${CHR}" >&2; exit 1; }

    N_ORIG=$(wc -l < "$COND")
    N_KEPT=$(wc -l < "${PRUNE}.prune.in")
    N_DROP=$(wc -l < "${PRUNE}.prune.out")
    [[ $((N_KEPT+N_DROP)) -eq "$N_ORIG" ]] || { echo "Pruning count mismatch: chr${CHR}" >&2; exit 1; }
    echo "[CONDITIONS] chr${CHR}: original=${N_ORIG}, retained=${N_KEPT}, removed=${N_DROP}"

    # COJO uses only the pruned CHIP-condition SNPs. The *tested* 15/5 QTN
    # remain unchanged and are still tested as in the original pilot.
    PREFIX="${OUT}/cojo_chr${CHR}_LDpruned"
    echo "[COJO] chr${CHR}: running conditional association with ${N_KEPT} CHIP SNPs"
    gcta64 \
      --bfile "$LD_BFILE" \
      --autosome-num 12 \
      --chr "$CHR" \
      --extract "${OUT}/chr${CHR}_selected.snps" \
      --cojo-file "$MA" \
      --cojo-cond "${PRUNE}.prune.in" \
      --cojo-wind 100000 \
      --thread-num 4 \
      --out "$PREFIX"

    if [[ ! -s "${PREFIX}.cma.cojo" && ! -s "${PREFIX}.cma" ]]; then
        echo "No GCTA cma result from chr${CHR}. Check ${PREFIX}.log" >&2
        exit 1
    fi
done

python3 "${WORK}/scripts/summarize_cojo_bb01_LDpruned.py" --work "$WORK" --out "$OUT"
REPORT="${OUT}/BB_outer01_20QTN_COJO_LDpruned_vs_marginal.tsv"
[[ -s "$REPORT" ]] || exit 1
[[ $(wc -l < "$REPORT") -eq 21 ]] || { echo 'Unexpected final report row count' >&2; exit 1; }
echo "DONE: $REPORT"
