#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH -t 00-24:00:00
#SBATCH -p shared

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

#SBATCH --array=1-8%8

#SBATCH -J NRC-MAF

#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=FAIL,END

#SBATCH -o slurm-NRC-%A_%a.out


set -euo pipefail


# ============================================================
# Input
# ============================================================

TRUE="Chr11_MAF-347_with_MAF.vcf.gz"

IMP="impute5-output-347-maf-ac-an.vcf.gz"

REGION_FILE="NRC_regions.tsv"

OUTDIR="NRC_chunks"


mkdir -p "${OUTDIR}"


# ============================================================
# Obtain region for this array task
# ============================================================

REGION=$(
    awk \
        -v id="${SLURM_ARRAY_TASK_ID}" \
        '$1==id {print $2}' \
        "${REGION_FILE}"
)


if [[ -z "${REGION}" ]]
then

    echo \
        "ERROR: no region found for array task ${SLURM_ARRAY_TASK_ID}" \
        >&2

    exit 1

fi


echo "========================================"
echo "NRC MAF analysis"
echo "========================================"

echo "Task:   ${SLURM_ARRAY_TASK_ID}"
echo "Region: ${REGION}"
echo "Start:  $(date)"

echo "========================================"


# Avoid accidental multithreading by numerical libraries
export OMP_NUM_THREADS=1


# ============================================================
# Run NRC
# ============================================================

python3 -u maf_NRC_global_GT.py \
    -t "${TRUE}" \
    -i "${IMP}" \
    -r "${REGION}" \
    -o "${OUTDIR}/NRC_chunk_${SLURM_ARRAY_TASK_ID}.summary.tsv" \
    --matrix-output "${OUTDIR}/NRC_chunk_${SLURM_ARRAY_TASK_ID}.matrix.tsv" \
    --info-threshold 0.8


echo
echo "========================================"
echo "Task completed"
echo "End: $(date)"
echo "========================================"
