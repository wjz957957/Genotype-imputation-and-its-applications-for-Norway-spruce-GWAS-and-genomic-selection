#!/bin/bash
#SBATCH -J cf11a_rrblup
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-9%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/11a_rrblup_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/11a_rrblup_%A_%a.err

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"

module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

TASK="${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array}"
(( TASK >= 0 && TASK < 10 )) || { echo "Invalid task ID: $TASK" >&2; exit 1; }
OUTER=$(printf '%02d' "$((TASK / 2 + 1))")
if (( TASK % 2 == 0 )); then TRAIT=BB; else TRAIT=H6; fi

SCRIPT="$WORK/Cross_family/scripts/11a_run_outer_rrblup.R"
OUT="$WORK/Cross_family/results/11_outer_rrblup/outer_${OUTER}/${TRAIT}"
[[ -s "$SCRIPT" ]] || { echo "Missing R script: $SCRIPT" >&2; exit 1; }
[[ -s "$WORK/Cross_family/results/10_outer_model_data/outer_${OUTER}/${TRAIT}/complete.flag" ]] || {
    echo "Missing/incomplete 10a model data for outer_${OUTER}/${TRAIT}" >&2; exit 1;
}
mkdir -p "$OUT"
Rscript "$SCRIPT" "$WORK" "$OUTER" "$TRAIT"
for F in complete.flag outer_metrics.tsv outer_test_predictions.tsv \
         outer_test_predictions.predictions_only.tsv model_manifest.tsv; do
    [[ -s "$OUT/$F" ]] || { echo "Missing RRBLUP output: $OUT/$F" >&2; exit 1; }
done
grep -qx 'status=PASS' "$OUT/complete.flag"
echo "PASS 11a outer_${OUTER}/${TRAIT}"
