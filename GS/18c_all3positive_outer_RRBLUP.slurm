#!/bin/bash
#SBATCH -J 18c_positive_RRBLUP
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-19%2
#SBATCH -o logs/18c_positive_%A_%a.out
#SBATCH -e logs/18c_positive_%A_%a.err
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
TASK=${SLURM_ARRAY_TASK_ID:?Submit as array}
(( TASK>=0 && TASK<20 )) || { echo 'Invalid task ID' >&2; exit 1; }
OUTER=$(printf '%02d' "$(( TASK/2+1 ))")
if (( TASK%2==0 )); then TRAIT=BB; else TRAIT=H6; fi
BASE="${WORK}/10_outer_test/08_final_genotypes_crosschip_r2_0p8_all3positive/outer_${OUTER}/${TRAIT}"
OUT="${WORK}/10_outer_test/09_outer_prediction_crosschip_r2_0p8_all3positive/outer_${OUTER}/${TRAIT}"
[[ -s "$BASE/complete.flag" && -s "$BASE/final_genotype_manifest.tsv" ]] || { echo 'Incomplete Step18b' >&2; exit 1; }
mkdir -p "$OUT"
if [[ -s "$OUT/complete.flag" && -s "$OUT/outer_metrics.tsv" && -s "$OUT/outer_test_predictions.tsv" ]]; then
  echo "[SKIP] already completed outer_${OUTER} $TRAIT"; exit 0
fi
rm -f "$OUT/complete.flag"
cd "$WORK"
Rscript "$WORK/scripts/18c_run_all3positive_outer_RRBLUP.R" "$WORK" "$OUTER" "$TRAIT"
[[ -s "$OUT/complete.flag" && -s "$OUT/outer_metrics.tsv" && -s "$OUT/outer_test_predictions.tsv" ]] || { echo 'Missing RR-BLUP outputs' >&2; exit 1; }
