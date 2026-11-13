#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -J GWAS_QQ_lambda
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o GWAS_QQ_lambda_%j.out
#SBATCH -e GWAS_QQ_lambda_%j.err

set -euo pipefail

BASE_DIR="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GWAS"

cd "${BASE_DIR}"

echo "============================================================"
echo "GWAS QQ plot and genomic inflation analysis"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "============================================================"

# Load required software stack
module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

echo
echo "Loaded modules:"
module list

echo
echo "Rscript:"
which Rscript
Rscript --version

# Check required packages
Rscript -e '
packages <- c("data.table", "ggplot2", "cowplot")
missing <- packages[!sapply(packages, requireNamespace, quietly=TRUE)]

if (length(missing) > 0) {
    cat("ERROR: Missing R packages:", paste(missing, collapse=", "), "\n")
    quit(status=1)
} else {
    cat("All required R packages are available.\n")
}
'

# Run analysis
Rscript GWAS_QQ_lambda.R

echo
echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
