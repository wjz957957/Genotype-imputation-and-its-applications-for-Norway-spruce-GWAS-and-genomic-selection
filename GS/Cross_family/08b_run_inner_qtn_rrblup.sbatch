#!/bin/bash
#SBATCH -J cf_qtn_rrblup
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 18:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --array=0-29%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/08b_qtn_rrblup_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/08b_qtn_rrblup_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03
# Protect rrBLUP::mixed.solve/LAPACK against observed multithreading segfaults.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
TASK="${SLURM_ARRAY_TASK_ID:?task required}"
[[ "$TASK" =~ ^[0-9]+$ ]] && ((TASK >= 0 && TASK < 30)) || { echo 'invalid task id' >&2; exit 2; }
OUTER=$(printf '%02d' "$((TASK/6+1))")
INNER=$(printf '%02d' "$((TASK%6/2+1))")
if ((TASK%2==0)); then TRAIT=BB; else TRAIT=H6; fi
INPUT="$WORK/Cross_family/results/08_inner_qtn_genotypes/outer_${OUTER}/inner_${INNER}/$TRAIT"
[[ -s "$INPUT/complete.flag" ]] && grep -qx 'status=PASS' "$INPUT/complete.flag" || { echo 'Stage08a input not PASS' >&2; exit 1; }
OUT="$WORK/Cross_family/results/08_inner_qtn_rrblup/outer_${OUTER}/inner_${INNER}/$TRAIT"
if [[ -s "$OUT/complete.flag" ]] && grep -qx 'status=PASS' "$OUT/complete.flag" &&
   [[ -s "$OUT/metrics.tsv" && -s "$OUT/validation_predictions.tsv" ]]; then
   echo "SKIP: already passed outer_${OUTER}/inner_${INNER}/$TRAIT"
   exit 0
fi
Rscript -e 'stopifnot(requireNamespace("rrBLUP",quietly=TRUE));cat("R:",as.character(getRversion())," rrBLUP:",as.character(packageVersion("rrBLUP")),"\n")'
Rscript "$WORK/Cross_family/scripts/08b_run_inner_qtn_rrblup.R" "$WORK" "$OUTER" "$INNER" "$TRAIT"
