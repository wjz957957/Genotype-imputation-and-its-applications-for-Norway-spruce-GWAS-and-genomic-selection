#!/bin/bash
#SBATCH -J cf_inner_rrblup
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-14%3
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/05_rrblup_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/05_rrblup_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "${WORK}"

# Use the PDC software stack specified for this project.
module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript is unavailable after loading the PDC R module" >&2
    exit 1
fi

Rscript -e '
cat("R version:", as.character(getRversion()), "\n")
if (!requireNamespace("rrBLUP", quietly=TRUE)) {
    stop("Package rrBLUP is unavailable in this R environment")
}
cat("rrBLUP version:", as.character(packageVersion("rrBLUP")), "\n")
'

export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if [[ ! "${TASK}" =~ ^[0-9]+$ ]] || (( TASK < 0 || TASK > 14 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 0 and 14" >&2
    exit 1
fi

OUTER=$(printf '%02d' "$((TASK / 3 + 1))")
INNER=$(printf '%02d' "$((TASK % 3 + 1))")

echo "========================================"
echo "Inner RRBLUP"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

Rscript "${WORK}/Cross_family/scripts/05_run_inner_rrblup.R" \
    "${WORK}" "${OUTER}" "${INNER}"
