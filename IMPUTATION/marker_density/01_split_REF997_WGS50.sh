#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=split_ref
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-12:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/split_ref_%A_%a.out
#SBATCH --error=logs/split_ref_%A_%a.err

set -euo pipefail

ml bcftools/1.20

# =========================
# Paths
# =========================

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

REFDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map

IDFILE=${WORKDIR}/WGS_50.ids

# Chromosomes included in this experiment
CHRS=(01 05 06 11)

chr=${CHRS[$SLURM_ARRAY_TASK_ID]}

INPUT=${REFDIR}/Chr${chr}_map_phased.bcf
OUTDIR=${WORKDIR}/Chr${chr}/base

mkdir -p "${OUTDIR}"

echo "================================================="
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome:   Chr${chr}"
echo "Input:        ${INPUT}"
echo "Output dir:   ${OUTDIR}"
echo "Start time:   $(date)"
echo "================================================="

# =========================
# 1. Extract 50 WGS validation individuals
# =========================

echo "[1/4] Extracting 50 validation individuals..."

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -S "${IDFILE}" \
    -Ob \
    -o "${OUTDIR}/Chr${chr}.WGS50.truth.bcf" \
    "${INPUT}"

# =========================
# 2. Index WGS50 truth
# =========================

echo "[2/4] Indexing WGS50 truth..."

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f \
    "${OUTDIR}/Chr${chr}.WGS50.truth.bcf"

# =========================
# 3. Remove the 50 validation individuals
#    to generate the 997-individual reference panel
# =========================

echo "[3/4] Generating REF997..."

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -S ^"${IDFILE}" \
    -Ob \
    -o "${OUTDIR}/Chr${chr}.REF997.bcf" \
    "${INPUT}"

# =========================
# 4. Index REF997
# =========================

echo "[4/4] Indexing REF997..."

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f \
    "${OUTDIR}/Chr${chr}.REF997.bcf"

# =========================
# QC
# =========================

echo
echo "========== Sample count QC =========="

N_INPUT=$(bcftools query -l "${INPUT}" | wc -l)
N_TRUTH=$(bcftools query -l "${OUTDIR}/Chr${chr}.WGS50.truth.bcf" | wc -l)
N_REF=$(bcftools query -l "${OUTDIR}/Chr${chr}.REF997.bcf" | wc -l)

echo "Original samples : ${N_INPUT}"
echo "WGS50 truth      : ${N_TRUTH}"
echo "REF997           : ${N_REF}"

if [[ "${N_INPUT}" -ne 1047 ]]; then
    echo "WARNING: Original sample number is not 1047!"
fi

if [[ "${N_TRUTH}" -ne 50 ]]; then
    echo "ERROR: WGS50 truth does not contain 50 samples!"
    exit 1
fi

if [[ "${N_REF}" -ne 997 ]]; then
    echo "ERROR: REF997 does not contain 997 samples!"
    exit 1
fi

if [[ $((N_TRUTH + N_REF)) -ne "${N_INPUT}" ]]; then
    echo "ERROR: REF + TRUTH sample counts do not sum to original!"
    exit 1
fi

echo
echo "Sample splitting passed."
echo "End time: $(date)"
echo "================================================="
