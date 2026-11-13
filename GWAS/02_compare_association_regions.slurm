#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH -t 00-02:30:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -J GWAS_region_cmp
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o GWAS_region_cmp_%j.out
#SBATCH -e GWAS_region_cmp_%j.err

set -euo pipefail

# ============================================================
# Step 2: Compare original and imputed association regions
#
# Run only AFTER Step 1 has completed successfully.
#
# Submit from the same GWAS root directory containing:
#   GWAS-Bud/
#   GWAS-DBH/
#   GWAS-FD/
#   compare_new_association_regions.R
#
# Expected main output:
#   association_region_summary.tsv
#
# Additional outputs:
#   BB_imputed_region_classification.tsv
#   DBH_imputed_region_classification.tsv
#   FD_imputed_region_classification.tsv
# ============================================================

module purge
module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

BASE_DIR="$(pwd)"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: ${BASE_DIR}"
echo "Rscript: $(which Rscript)"
echo

if [[ ! -f "${BASE_DIR}/compare_new_association_regions.R" ]]; then
    echo "ERROR: compare_new_association_regions.R not found in ${BASE_DIR}"
    exit 1
fi

# Confirm that Step 1 produced the six expected clumped files.
required_files=(
  "GWAS-Bud/gwas_clumps_revised_origin_Bud.clumped"
  "GWAS-Bud/gwas_clumps_revised_imputed_Bud.clumped"
  "GWAS-DBH/gwas_clumps_revised_origin_DBH.clumped"
  "GWAS-DBH/gwas_clumps_revised_imputed_DBH.clumped"
  "GWAS-FD/gwas_clumps_revised_origin_FD.clumped"
  "GWAS-FD/gwas_clumps_revised_imputed_FD.clumped"
)

for f in "${required_files[@]}"; do
    if [[ ! -f "${BASE_DIR}/${f}" ]]; then
        echo "ERROR: Missing Step-1 output: ${BASE_DIR}/${f}"
        exit 1
    fi
done

Rscript "${BASE_DIR}/compare_new_association_regions.R"

echo
echo "============================================================"
echo "Step 2 finished successfully: $(date)"
echo "============================================================"

echo
echo "Main summary:"
cat "${BASE_DIR}/association_region_summary.tsv"

