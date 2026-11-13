#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=imp_r1493
#SBATCH -t 00-03:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/imp_r1493_%A_%a.out
#SBATCH --error=logs/imp_r1493_%A_%a.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

REP=$1

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
MAPDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map

IMPUTE5=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/impute5_v1.2.1_x86_64_static

REF=${WORKDIR}/Chr11/base/Chr11.REF997.bcf
TARGET=${WORKDIR}/Chr11/random1493/rep${REP}/Chr11.random1493.rep${REP}.target.phased.bcf

MAP=${MAPDIR}/LD.GenMap.chr11.txt
REGION_FILE=${WORKDIR}/Chr11/Chr11.impute_regions.tsv

OUTDIR=${WORKDIR}/Chr11/random1493/rep${REP}/imputed_chunks

mkdir -p "${OUTDIR}"
mkdir -p "${WORKDIR}/logs"

LINE_NO=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${LINE_NO}p" "${REGION_FILE}")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: no region for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

read -r CHUNK REGION BUFFER <<< "${LINE}"

OUT=${OUTDIR}/Chr11.random1493.rep${REP}.chunk${CHUNK}.vcf.gz

echo "============================================================"
echo "Chr        : 11"
echo "Markers    : 1493"
echo "Replicate  : ${REP}"
echo "Chunk      : ${CHUNK}"
echo "Region     : ${REGION}"
echo "Buffer     : ${BUFFER}"
echo "Reference  : ${REF}"
echo "Target     : ${TARGET}"
echo "============================================================"

for f in "${REF}" "${TARGET}" "${MAP}" "${REGION_FILE}"
do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: missing ${f}"
        exit 1
    fi
done

TARGET_N=$(bcftools index -n "${TARGET}")
TARGET_SAMPLES=$(bcftools query -l "${TARGET}" | wc -l)

if [[ "${TARGET_N}" -ne 1493 ]]; then
    echo "ERROR: target marker count is ${TARGET_N}, expected 1493."
    exit 1
fi

if [[ "${TARGET_SAMPLES}" -ne 50 ]]; then
    echo "ERROR: target does not contain 50 samples."
    exit 1
fi

"${IMPUTE5}" \
    --h "${REF}" \
    --g "${TARGET}" \
    --m "${MAP}" \
    --o "${OUT}" \
    --r "${REGION}" \
    --buffer-region "${BUFFER}" \
    --threads "${SLURM_CPUS_PER_TASK}"

bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f "${OUT}"

NSAMPLE=$(bcftools query -l "${OUT}" | wc -l)
NVAR=$(bcftools index -n "${OUT}")

echo
echo "Samples : ${NSAMPLE}"
echo "Variants: ${NVAR}"

if [[ "${NSAMPLE}" -ne 50 ]]; then
    echo "ERROR: imputed output does not contain 50 samples."
    exit 1
fi

echo "Finished: $(date)"
