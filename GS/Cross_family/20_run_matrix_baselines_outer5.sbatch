#!/bin/bash
#SBATCH -J cf20_matrix_base
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/20_matrix_baselines_%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/20_matrix_baselines_%j.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"
mkdir -p "$WORK/Cross_family/logs"

module load bioinfo-tools
if ! command -v Rscript >/dev/null 2>&1; then
    module load R
fi
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript is unavailable" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is unavailable" >&2; exit 1; }
Rscript -e 'if (!requireNamespace("rrBLUP", quietly=TRUE)) stop("rrBLUP package unavailable"); cat("rrBLUP", as.character(packageVersion("rrBLUP")), "\n")'

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
SCRIPT="$WORK/Cross_family/scripts/20_run_matrix_baselines_outer5.R"
VERIFY="$WORK/Cross_family/scripts/20b_verify_matrix_baselines_outer5.py"
[[ -s "$SCRIPT" ]] || { echo "ERROR: missing $SCRIPT" >&2; exit 1; }
[[ -s "$VERIFY" ]] || { echo "ERROR: missing $VERIFY" >&2; exit 1; }

echo "===== Validate current 5-fold inputs ====="
Rscript "$SCRIPT" --work "$WORK" --validate-only

echo "===== Fit ABLUP, RRBLUP-array, RRBLUP-imputed ====="
Rscript "$SCRIPT" --work "$WORK"

echo "===== Verify outputs ====="
python3 "$VERIFY" --work "$WORK"
echo "OVERALL PASS: 5-fold matrix baselines ready for seven-model summary."
