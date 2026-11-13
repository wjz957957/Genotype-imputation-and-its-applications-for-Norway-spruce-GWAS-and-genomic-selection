#!/bin/bash
#SBATCH -J 12a_BB01_innerLD
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-2%3
#SBATCH -o logs/12a_BB01_innerLD_%A_%a.out
#SBATCH -e logs/12a_BB01_innerLD_%A_%a.err

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
export LC_ALL=C

module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
PY="${WORK}/scripts/12a_filter_inner_qtn_ld.py"
MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
OUTER=01
INNER=$(printf '%02d' $(( ${SLURM_ARRAY_TASK_ID} + 1 )))
TRAIT=BB
R2=${R2:-0.8}
R2_TAG=${R2//./p}
FINAL="${WORK}/08_imputed_GWAS/08_crosschip_LD_pilot/outer_${OUTER}/inner_${INNER}/${TRAIT}/r2_${R2_TAG}"
RANKED="${WORK}/08_imputed_GWAS/07_imputed_marker_sets/outer_${OUTER}/inner_${INNER}/${TRAIT}/ranked_imputed_snps.tsv"
BEST_FILE="${WORK}/07_chip_RRBLUP/optimization/outer_${OUTER}/${TRAIT}_best_marker.txt"
CV="${WORK}/02_CV/outer_${OUTER}/inner_${INNER}"
TRAIN="${CV}/inner_train.samples"
VALID="${CV}/inner_valid.samples"
OUTER_TEST="${WORK}/02_CV/outer_${OUTER}/outer_test.samples"

for F in "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam" \
         "$PY" "$RANKED" "$BEST_FILE" "$TRAIN" "$VALID" "$OUTER_TEST"; do
  [[ -s "$F" ]] || { echo "ERROR: missing $F" >&2; exit 1; }
done

BEST=$(tr -d '\r\n[:space:]' < "$BEST_FILE")
[[ "$BEST" =~ ^chip_(50|100|200|300|400|500|1000|2000|3000|4000|5000|8000|10000|20000|30000|40000|LD_ALL|ALL)$ ]] || {
  echo "ERROR: unexpected best CHIP marker setting: $BEST" >&2; exit 1;
}
CHIP="${WORK}/06_chip_clumping/outer_${OUTER}/inner_${INNER}/${TRAIT}/marker_sets/${BEST}.snps"
[[ -s "$CHIP" ]] || { echo "ERROR: missing inner_train CHIP list $CHIP" >&2; exit 1; }

mkdir -p "$FINAL"
if [[ -s "$FINAL/complete.flag" && -s "$FINAL/top1000_extra.snps" && -s "$FINAL/filter_summary.tsv" ]]; then
  echo "[SKIP] Already completed: $FINAL"
  exit 0
fi

[[ "$R2" == "0.8" || "$R2" == "0.5" || "$R2" == "0.2" ]] || {
    echo 'ERROR: R2 must be one of 0.8, 0.5, 0.2 (first pilot use 0.8)' >&2; exit 1;
}
N_RANKED=$(( $(wc -l < "$RANKED") - 1 ))
(( N_RANKED >= 1000 )) || { echo "ERROR: too few ranked SNPs: $N_RANKED" >&2; exit 1; }

printf '===== BB outer_%s inner_%s threshold=%s =====\n' "$OUTER" "$INNER" "$R2"
printf 'Best CHIP: %s (%s markers)\n' "$BEST" "$(wc -l < "$CHIP")"
printf 'Total clumped candidates: %s\n' "$N_RANKED"
echo 'Only inner_train genotypes; no association rerun; no original output overwritten.'

# Expand only when needed. The final step uses the FIRST pool with >=1000 survivors.
POOL_CHOICES=(10000 20000 40000 80000 100000)
LAST_POOL=0
for REQUESTED in "${POOL_CHOICES[@]}"; do
  POOL=$REQUESTED
  (( POOL > N_RANKED )) && POOL=$N_RANKED
  (( POOL > LAST_POOL )) || continue
  LAST_POOL=$POOL
  STAGE="${FINAL}/stages/top${POOL}"
  mkdir -p "$STAGE"
  echo "===== Candidate pool $POOL of $N_RANKED ====="

  python3 "$PY" prepare \
      --ranked "$RANKED" --chip "$CHIP" \
      --train "$TRAIN" --valid "$VALID" --outer-test "$OUTER_TEST" \
      --pool "$POOL" --stage "$STAGE"

  LD_BFILE="${STAGE}/train_chip_plus_candidates"
  plink --bfile "$MASTER" --keep "${STAGE}/train.keep" \
      --extract "${STAGE}/union.snps" \
      --maf 0.01 --geno 0.05 --keep-allele-order \
      --make-bed --threads "${SLURM_CPUS_PER_TASK:-4}" \
      --memory 36000 --out "$LD_BFILE"

  LD_PREFIX="${STAGE}/selected_chip_to_qtn"
  # CRITICAL: sparse LD threshold reports only pairs r2>=R2; unreported
  # pairs are interpreted as BELOW R2, NOT as actual r2=0.
  # --ld-snp-list restricts one partner to selected CHIP. PLINK computes
  # same-chromosome pairs, without the earlier 1 Mb positional limit.
  plink --bfile "$LD_BFILE" --ld-snp-list "${STAGE}/selected_chip.snps" \
      --r2 gz --ld-window 999999 --ld-window-kb 10000000 \
      --ld-window-r2 "$R2" \
      --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 36000 \
      --out "$LD_PREFIX"

  python3 "$PY" report --ranked "$RANKED" --stage "$STAGE" \
      --bim "${LD_BFILE}.bim" --fam "${LD_BFILE}.fam" \
      --ld "${LD_PREFIX}.ld.gz" --train "$TRAIN" \
      --pool "$POOL" --threshold "$R2"

  N_RETAIN=$(awk -F '\t' '$1=="Retained"{print $2}' "$STAGE/filter_summary.tsv")
  echo "[PROGRESS] inner_${INNER}: pool=$POOL retained=$N_RETAIN"
  if (( N_RETAIN >= 1000 )); then
    python3 "$PY" publish --stage "$STAGE" --out "$FINAL"
    echo "[DONE] $FINAL"
    exit 0
  fi
  if (( POOL == N_RANKED )); then
    break
  fi
done

printf '%s\n' "[INSUFFICIENT] Only ${N_RETAIN} candidates after screening up to pool ${LAST_POOL}." \
    "No threshold was relaxed, and no original marker-set files were modified." \
    > "${FINAL}/insufficient_candidates.txt"
echo "ERROR: Cannot obtain 1000 QTN after max pool ${LAST_POOL}" >&2
exit 3
