#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -J GWAS_LD_prune
#SBATCH -p shared
#SBATCH -t 12:00:00
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH --array=0-11%4
#SBATCH -o logs/LD_prune_%A_%a.out
#SBATCH -e logs/LD_prune_%A_%a.err

set -euo pipefail


# ============================================================
# 1. Working directory
# ============================================================

BASE_DIR="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GWAS"

cd "${BASE_DIR}"

mkdir -p logs
mkdir -p LD_threshold_results


# ============================================================
# 2. Load PLINK
#
# Modify this section if your cluster uses another PLINK module.
# If plink is already available in PATH, you can comment it out.
# ============================================================

module purge
module load plink

echo "============================================================"
echo "Job ID:        ${SLURM_JOB_ID}"
echo "Array task:    ${SLURM_ARRAY_TASK_ID}"
echo "Node:          $(hostname)"
echo "Start time:    $(date)"
echo "PLINK:         $(which plink)"
echo "============================================================"


# ============================================================
# 3. Define the six actual GWAS datasets
#
# Order:
# 0  Bud original
# 1  Bud imputed
# 2  DBH original
# 3  DBH imputed
# 4  FD original
# 5  FD imputed
# ============================================================

TRAITS=(
    "Bud"
    "Bud"
    "DBH"
    "DBH"
    "FD"
    "FD"
)

DATASETS=(
    "original"
    "imputed"
    "original"
    "imputed"
    "original"
    "imputed"
)

PREFIXES=(
    "GWAS-Bud/Bud-origin-plink"
    "GWAS-Bud/Bud-imputed-plink"
    "GWAS-DBH/DBH-origin-plink"
    "GWAS-DBH/DBH-imputed-plink"
    "GWAS-FD/FD-origin-plink"
    "GWAS-FD/FD-imputed-plink"
)


# ============================================================
# 4. Two LD-pruning thresholds
#
# Tasks 0-5  : r2 = 0.8
# Tasks 6-11 : r2 = 0.5
# ============================================================

N_DATASETS=6

DATA_INDEX=$(( SLURM_ARRAY_TASK_ID % N_DATASETS ))

if [ "${SLURM_ARRAY_TASK_ID}" -lt 6 ]; then
    R2="0.8"
    R2_TAG="08"
else
    R2="0.5"
    R2_TAG="05"
fi


TRAIT="${TRAITS[$DATA_INDEX]}"
DATASET="${DATASETS[$DATA_INDEX]}"
BFILE="${PREFIXES[$DATA_INDEX]}"


# ============================================================
# 5. Check PLINK input files
# ============================================================

for ext in bed bim fam; do

    if [ ! -f "${BFILE}.${ext}" ]; then
        echo "ERROR: missing file ${BFILE}.${ext}"
        exit 1
    fi

done


# ============================================================
# 6. Output prefix
# ============================================================

OUTDIR="LD_threshold_results/${TRAIT}_${DATASET}"

mkdir -p "${OUTDIR}"

OUTPREFIX="${OUTDIR}/${TRAIT}_${DATASET}_LD_r2_${R2_TAG}"


echo
echo "Trait:              ${TRAIT}"
echo "Dataset:            ${DATASET}"
echo "Input prefix:       ${BFILE}"
echo "LD window:          1000 kb"
echo "Step size:          1 variant"
echo "r2 threshold:       ${R2}"
echo "Output prefix:      ${OUTPREFIX}"
echo


# ============================================================
# 7. Count tested/input SNPs
# ============================================================

N_TESTED=$(wc -l < "${BFILE}.bim")

echo "Number of SNPs in BIM: ${N_TESTED}"


# ============================================================
# 8. LD pruning
#
# Same procedure as the original analysis:
#
# --indep-pairwise 1000kb 1 r2
# ============================================================

plink \
    --bfile "${BFILE}" \
    --indep-pairwise 1000kb 1 "${R2}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --out "${OUTPREFIX}"


# ============================================================
# 9. Count LD-pruned markers
# ============================================================

PRUNE_FILE="${OUTPREFIX}.prune.in"

if [ ! -s "${PRUNE_FILE}" ]; then
    echo "ERROR: ${PRUNE_FILE} does not exist or is empty."
    exit 1
fi

N_PRUNED=$(wc -l < "${PRUNE_FILE}")

echo "LD-pruned SNP number: ${N_PRUNED}"


# ============================================================
# 10. Calculate Bonferroni threshold
#
# threshold = 0.05 / N_LD-pruned
# ============================================================

P_THRESHOLD=$(awk -v n="${N_PRUNED}" \
    'BEGIN {printf "%.10e", 0.05/n}')


echo "Bonferroni threshold: ${P_THRESHOLD}"


# ============================================================
# 11. Write one summary file for this task
# ============================================================

SUMMARY="${OUTDIR}/${TRAIT}_${DATASET}_r2_${R2_TAG}_summary.tsv"

{
    echo -e "Trait\tDataset\tR2_cutoff\tTested_SNPs\tLD_pruned_SNPs\tBonferroni_threshold"
    echo -e "${TRAIT}\t${DATASET}\t${R2}\t${N_TESTED}\t${N_PRUNED}\t${P_THRESHOLD}"
} > "${SUMMARY}"


echo
echo "Summary:"
cat "${SUMMARY}"

echo
echo "Finished: $(date)"
echo "============================================================"
