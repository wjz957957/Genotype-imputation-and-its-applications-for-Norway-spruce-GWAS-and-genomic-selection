#!/bin/bash
#SBATCH -J 18d_positive_sum
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o logs/18d_positive_%j.out
#SBATCH -e logs/18d_positive_%j.err
set -Eeuo pipefail
module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
Rscript scripts/18d_summarize_all3positive_outer.R "$WORK"
