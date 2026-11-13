#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=Chr11_D25_full
#SBATCH --array=0-9
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-02:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/Chr11_D25_full_%A_%a.out
#SBATCH --error=logs/Chr11_D25_full_%A_%a.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

IMPUTE5=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/impute5_v1.2.1_x86_64_static

REF=${WORKDIR}/Chr11/base/Chr11.REF997.bcf

TARGET=${WORKDIR}/Chr11/D25/rep1/Chr11.D25.rep1.target.phased.bcf

MAP=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map/LD.GenMap.chr11.txt

OUTDIR=${WORKDIR}/Chr11/D25/rep1/imputed_chunks

mkdir -p "${OUTDIR}"

# ==========================================
# Imputation regions
# ==========================================

REGIONS=(
"83036-100000000"
"100000001-200000000"
"200000001-300000000"
"300000001-400000000"
"400000001-500000000"
"500000001-600000000"
"600000001-700000000"
"700000001-800000000"
"800000001-900000000"
"900000001-954573755"
)

# Buffer regions
BUFFERS=(
"1-100050000"
"99950001-200050000"
"199950001-300050000"
"299950001-400050000"
"399950001-500050000"
"499950001-600050000"
"599950001-700050000"
"699950001-800050000"
"799950001-900050000"
"899950001-954579047"
)

REGION=${REGIONS[$SLURM_ARRAY_TASK_ID]}
BUFFER=${BUFFERS[$SLURM_ARRAY_TASK_ID]}

OUT=${OUTDIR}/Chr11.D25.rep1.chunk${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "==============================================="
echo "Task             : ${SLURM_ARRAY_TASK_ID}"
echo "Imputation region: 11:${REGION}"
echo "Buffer region    : 11:${BUFFER}"
echo "Reference        : ${REF}"
echo "Target           : ${TARGET}"
echo "Output           : ${OUT}"
echo "Start            : $(date)"
echo "==============================================="

"${IMPUTE5}" \
    --h "${REF}" \
    --g "${TARGET}" \
    --m "${MAP}" \
    --o "${OUT}" \
    --r 11:${REGION} \
    --buffer-region 11:${BUFFER} \
    --threads ${SLURM_CPUS_PER_TASK}

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${OUT}"

N_SAMPLE=$(bcftools query -l "${OUT}" | wc -l)
N_VARIANT=$(bcftools index -n "${OUT}")

echo
echo "Samples : ${N_SAMPLE}"
echo "Variants: ${N_VARIANT}"
echo "Finished: $(date)"

if [[ "${N_SAMPLE}" -ne 50 ]]; then
    echo "ERROR: output does not contain 50 samples."
    exit 1
fi
