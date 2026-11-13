#!/bin/bash
#SBATCH -J 13c_BB01_LD08
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/13c_BB01_LD08_%j.out
#SBATCH -e logs/13c_BB01_LD08_%j.err

set -Eeuo pipefail

module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

SCRIPT="${WORK}/scripts/13c_BB01_crossLD08_summarize.R"

[[ -s "${SCRIPT}" ]] || {
    echo "ERROR: Missing R script" >&2
    exit 1
}

cd "${WORK}"

Rscript "${SCRIPT}" "${WORK}"

echo "DONE Step13c pilot summary"
