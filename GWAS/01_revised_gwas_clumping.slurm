#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH -t 00-06:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -J GWAS_clump_rev
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o GWAS_clump_rev_%j.out
#SBATCH -e GWAS_clump_rev_%j.err

set -euo pipefail

# ============================================================
# Step 1: Revised PLINK clumping for Reviewer #10
#
# Submit this job from the GWAS root directory containing:
#   GWAS-Bud/
#   GWAS-DBH/
#   GWAS-FD/
#
# Original clumping settings retained:
#   --clump-kb 50
#   --clump-r2 0.2
#   --clump-p2 1e-5
#
# Dataset-specific --clump-p1 thresholds:
#   BB original   9.1044830475e-07
#   BB imputed    9.6525655554e-08
#   DBH original  7.5243412439e-07
#   DBH imputed   8.8580079049e-08
#   FD original   6.7133017361e-07
#   FD imputed    8.5449489952e-08
#
# New output names are used, so old clumping files are not overwritten.
# ============================================================

module purge
module load PDC/26.03
module load bioinfo-tools plink/1.90b4.9

BASE_DIR="$(pwd)"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: ${BASE_DIR}"
echo "PLINK: $(which plink)"
echo

for d in GWAS-Bud GWAS-DBH GWAS-FD; do
    if [[ ! -d "${BASE_DIR}/${d}" ]]; then
        echo "ERROR: ${BASE_DIR}/${d} not found."
        echo "Submit this job from the GWAS root directory."
        exit 1
    fi
done

run_clump () {
    local trait="$1"
    local dataset="$2"
    local subdir="$3"
    local bfile="$4"
    local mlma="$5"
    local p1="$6"
    local outprefix="$7"

    echo
    echo "============================================================"
    echo "Trait=${trait}  Dataset=${dataset}"
    echo "p1=${p1}"
    echo "============================================================"

    cd "${BASE_DIR}/${subdir}"

    plink \
      --bfile "${bfile}" \
      --clump "${mlma}" \
      --clump-field p \
      --clump-snp-field SNP_ID \
      --clump-p1 "${p1}" \
      --clump-p2 1e-5 \
      --clump-r2 0.2 \
      --clump-kb 50 \
      --out "${outprefix}"

    cd "${BASE_DIR}"
}

# BB
run_clump \
  "BB" "Original" "GWAS-Bud" \
  "Bud-origin-plink" \
  "Bud-origin-GWAS.loco_fixed.mlma" \
  "9.1044830475e-07" \
  "gwas_clumps_revised_origin_Bud"

run_clump \
  "BB" "Imputed" "GWAS-Bud" \
  "Bud-imputed-plink" \
  "Bud-imputed-PC0-GWASloco_chrall_fixed.mlma" \
  "9.6525655554e-08" \
  "gwas_clumps_revised_imputed_Bud"

# DBH
run_clump \
  "DBH" "Original" "GWAS-DBH" \
  "DBH-origin-plink" \
  "DBH-origin-GWAS.loco_fixed.mlma" \
  "7.5243412439e-07" \
  "gwas_clumps_revised_origin_DBH"

run_clump \
  "DBH" "Imputed" "GWAS-DBH" \
  "DBH-imputed-plink" \
  "DBH-imputed-PC0-GWASloco_chrall_fixed.mlma" \
  "8.8580079049e-08" \
  "gwas_clumps_revised_imputed_DBH"

# FD
run_clump \
  "FD" "Original" "GWAS-FD" \
  "FD-origin-plink" \
  "FD-origin-GWAS.loco_fixed.mlma" \
  "6.7133017361e-07" \
  "gwas_clumps_revised_origin_FD"

run_clump \
  "FD" "Imputed" "GWAS-FD" \
  "FD-imputed-plink" \
  "FD-imputed-PC0-GWASloco_chrall_fixed.mlma" \
  "8.5449489952e-08" \
  "gwas_clumps_revised_imputed_FD"

echo
echo "============================================================"
echo "Step 1 finished successfully: $(date)"
echo "============================================================"

echo
echo "Revised .clumped files:"
find "${BASE_DIR}"/GWAS-Bud "${BASE_DIR}"/GWAS-DBH "${BASE_DIR}"/GWAS-FD \
    -maxdepth 1 -name 'gwas_clumps_revised_*.clumped' -print

