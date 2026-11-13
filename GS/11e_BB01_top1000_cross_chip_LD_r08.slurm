#!/bin/bash
#SBATCH -J 11e_BB01_crossLD08
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/11e_BB01_crossLD08_%j.out
#SBATCH -e logs/11e_BB01_crossLD08_%j.err

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
export LC_ALL=C

module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
SCRIPT="${WORK}/scripts/11e_BB01_top1000_cross_chip_LD_r08.py"
MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
TRAIN="${WORK}/10_outer_test/02_imputed_inputs/outer_01/outer_train.keep"
CHIP="${WORK}/10_outer_test/01_chip_selection/outer_01/selected/BB_best_chip.snps"
RANKED="${WORK}/10_outer_test/04_imputed_clumping/outer_01/BB/ranked_extra_snps.tsv"
OUT="${WORK}/10_outer_test/07_cross_chip_qtn_LD_pilot/outer_01/BB/r2_0p8_top1000"
LD_BFILE="${OUT}/selected_1100_train_qc"
LD_PREFIX="${OUT}/selected_1100_train_LD"

# Reuse the existing outer_train genotypes, the 100 selected CHIP SNPs,
# and the existing MLMA-LOCO GWAS + PLINK clumping ranking (top 1000).
# No GWAS, clumping or phenotype-based selection is rerun.
for f in "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam" \
         "$TRAIN" "$CHIP" "$RANKED" "$SCRIPT"; do
    [[ -s "$f" ]] || { echo "ERROR: missing/empty input: $f" >&2; exit 1; }
done

mkdir -p "$OUT"
if [[ -s "${OUT}/complete.flag" && -s "${OUT}/filter_summary.tsv" && \
      -s "${OUT}/top1000_QTN_vs_selected_CHIP_LD.tsv" ]]; then
    echo "[SKIP] Completed LD diagnosis: $OUT"
    exit 0
fi
rm -f "${OUT}/complete.flag"

printf '=== BB outer_01 top1000 QTN / selected CHIP LD filter ===\n'
printf 'TRAIN=%s\nCHIP=%s\nRANKED=%s\nOUT=%s\n' "$TRAIN" "$CHIP" "$RANKED" "$OUT"
printf 'Remove QTN when max same-chromosome CHIP LD r2 >= 0.8\n'

python3 "$SCRIPT" prepare --ranked "$RANKED" --chip "$CHIP" --out "$OUT"

# Exactly the sample restriction and SNP QC used by Step 10b1:
# --keep outer_train, --maf 0.01, --geno 0.05; genotype source is
# the same all_904_genotypes as Step 10b1 and previous 11a pilot.
plink \
  --bfile "$MASTER" \
  --keep "$TRAIN" \
  --extract "${OUT}/all_1100_for_LD.snps" \
  --maf 0.01 --geno 0.05 --keep-allele-order \
  --make-bed \
  --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 24000 \
  --out "$LD_BFILE"

# Compute all within-chromosome SNP pairs (up to genome-wide distances);
# the Python report ONLY uses CHIP-QTN pairs. The 10,000,000 kb window
# comfortably exceeds the length of the longest spruce chromosome here.
plink \
  --bfile "$LD_BFILE" \
  --r2 \
  --ld-window 999999 \
  --ld-window-kb 10000000 \
  --ld-window-r2 0 \
  --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 24000 \
  --out "$LD_PREFIX"

python3 "$SCRIPT" report \
  --ranked "$RANKED" \
  --bim "${LD_BFILE}.bim" \
  --fam "${LD_BFILE}.fam" \
  --ld "${LD_PREFIX}.ld" \
  --out "$OUT" \
  --threshold 0.8

printf 'DONE BB outer_01 r2=0.8 top1000 cross-CHIP LD pilot\n' > "${OUT}/complete.flag"
echo '=========================================='
echo '[COMPLETE] Summary:'
cat "${OUT}/filter_summary.tsv"
echo '[COMPLETE] Candidate counts:'
cat "${OUT}/candidate_count_availability.tsv"
echo '[COMPLETE] Original top20 diagnostic:'
cat "${OUT}/original_top20_LD_diagnostic.tsv"
