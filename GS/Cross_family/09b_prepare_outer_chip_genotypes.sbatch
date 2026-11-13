#!/bin/bash
#SBATCH -J cf09_chip_prep
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=0-4%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09b_chip_prep_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09b_chip_prep_%A_%a.err
set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 09b line=${LINENO} exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9
PLINK_BIN=$(command -v plink)
PY_BIN=$(command -v python3)
VALIDATOR="$WORK/Cross_family/scripts/09b_validate_outer_chip_inputs.py"
TASK=${SLURM_ARRAY_TASK_ID:?Submit with sbatch array}
[[ "$TASK" =~ ^[0-4]$ ]] || { echo "Invalid task: $TASK" >&2; exit 1; }
OUTER=$(printf '%02d' "$((TASK+1))")
IN="$WORK/Cross_family/results/09_outer_inputs/outer_${OUTER}"
OUT="$WORK/Cross_family/results/09_outer_chip_genotypes/outer_${OUTER}"
MASTER="$WORK/08_imputed_GWAS/00_reference/all_904_genotypes"
PREFIX="$OUT/chip_ALL_outer_train"
for ext in bed bim fam; do [[ -s "${MASTER}.${ext}" ]] || { echo "missing master ${ext}" >&2; exit 1; }; done
[[ -s "$VALIDATOR" ]] || { echo "missing validator: $VALIDATOR" >&2; exit 1; }
mkdir -p "$OUT"
"$PY_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --out "$OUT" --prefix "$PREFIX" --phase pre
if [[ -s "$OUT/complete.flag" ]]; then
    grep -qx 'status=PASS' "$OUT/complete.flag" &&
    grep -qx "outer=$OUTER" "$OUT/complete.flag" || { echo "invalid prior flag" >&2; exit 1; }
    "$PY_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --out "$OUT" --prefix "$PREFIX" --phase post
    echo "[SKIP] 09b outer_${OUTER} complete and verified"
    exit 0
fi
# Retain the EXACT 43,267 chip roster. QC and association/clumping must be
# separately matched to the original Stage03 chip-selection method.
# Do not mix imputed QTN into the chip GWAS input.
TMP="$OUT/.chip_ALL_outer_train_${SLURM_JOB_ID:-manual}_${TASK}"
rm -f "${TMP}."{bed,bim,fam,log,nosex}
"$PLINK_BIN" --bfile "$MASTER" \
    --keep "$IN/outer_train.keep" \
    --extract "$IN/chip_ALL.snps" \
    --keep-allele-order --allow-extra-chr --make-bed \
    --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 18000 \
    --out "$TMP"
"$PY_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --out "$OUT" --prefix "$TMP" --phase post
# Publish only after the extracted SNP and sample identities pass validation.
for ext in bed bim fam log nosex; do
    if [[ -f "${TMP}.${ext}" ]]; then mv -f "${TMP}.${ext}" "${PREFIX}.${ext}"; fi
done
"$PY_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --out "$OUT" --prefix "$PREFIX" --phase post
N_TRAIN=$(wc -l < "${PREFIX}.fam")
N_CHIP=$(wc -l < "${PREFIX}.bim")
cat > "$OUT/input_summary.tsv" <<EOF
item\tvalue
outer\t$OUTER
training_samples\t$N_TRAIN
chip_snps_exact\t$N_CHIP
outer_test_phenotypes_read\tNO
qtn_genotype_included\tNO
chip_gwas_performed\tNO
chip_selection_performed\tNO
EOF
# Correct tabs for portable TSV parsing.
sed -i 's/\\t/\t/g' "$OUT/input_summary.tsv"
cat > "$OUT/complete.flag" <<EOF
status=PASS
outer=$OUTER
training_only=PASS
outer_test_phenotype_read=NO
chip_ALL=43267
chip_GWAS_not_run=YES
EOF
echo "PASS 09b: outer_${OUTER}, samples=$N_TRAIN, chip_SNPs=$N_CHIP"
