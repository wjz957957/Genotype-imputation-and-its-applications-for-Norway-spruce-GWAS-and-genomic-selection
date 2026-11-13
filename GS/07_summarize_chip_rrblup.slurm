#!/bin/bash
#SBATCH -J summarize_chip
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/07_summarize_chip_%j.out
#SBATCH -e logs/07_summarize_chip_%j.err

set -euo pipefail

module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

cd "${WORK}"

Rscript scripts/07_summarize_chip_rrblup.R "${WORK}"

echo "DONE summarize chip RRBLUP"
