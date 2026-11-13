#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH -t 00-12:00:00
#SBATCH -p shared

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

#SBATCH --array=1-8%8

#SBATCH -J MAF-ACC

#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=FAIL,END

#SBATCH -o slurm-MAFACC-%A_%a.out


set -euo pipefail


TRUE="Chr11_MAF-347_with_MAF.vcf.gz"

IMPUTED="impute5-output-347-maf-ac-an.vcf.gz"

REGION_FILE="NRC_regions.tsv"

OUTDIR="MAF_accuracy_chunks"


mkdir -p "${OUTDIR}"


# ============================================================
# Region corresponding to array task
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
        "ERROR: no region found for task ${SLURM_ARRAY_TASK_ID}" \
        >&2

    exit 1

fi


echo "========================================"
echo "Unified MAF accuracy analysis"
echo "========================================"

echo "Task:   ${SLURM_ARRAY_TASK_ID}"
echo "Region: ${REGION}"
echo "Start:  $(date)"

echo "========================================"


export OMP_NUM_THREADS=1


python3 -u maf_accuracy_global_unified.py \
    -t "${TRUE}" \
    -i "${IMPUTED}" \
    -r "${REGION}" \
    -o "${OUTDIR}/accuracy_chunk_${SLURM_ARRAY_TASK_ID}.stats.tsv" \
    --info-threshold 0.8


echo
echo "========================================"
echo "Task completed"
echo "End: $(date)"
echo "========================================"
