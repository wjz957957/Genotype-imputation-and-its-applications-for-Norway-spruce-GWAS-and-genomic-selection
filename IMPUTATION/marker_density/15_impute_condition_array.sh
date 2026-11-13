#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=impute_condition
#SBATCH -t 00-03:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/impute_%A_%a.out
#SBATCH --error=logs/impute_%A_%a.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20


CHR2=$1
DENSITY=$2
REP=$3

CHRNUM=$((10#${CHR2}))


WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

MAPDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map

IMPUTE5=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/impute5_v1.2.1_x86_64_static


REF=${WORKDIR}/Chr${CHR2}/base/Chr${CHR2}.REF997.bcf

TARGET=${WORKDIR}/Chr${CHR2}/D${DENSITY}/rep${REP}/Chr${CHR2}.D${DENSITY}.rep${REP}.target.phased.bcf

MAP=${MAPDIR}/LD.GenMap.chr${CHRNUM}.txt

REGION_FILE=${WORKDIR}/Chr${CHR2}/Chr${CHR2}.impute_regions.tsv

OUTDIR=${WORKDIR}/Chr${CHR2}/D${DENSITY}/rep${REP}/imputed_chunks

mkdir -p "${OUTDIR}"


# ============================================================
# Read chromosome-specific region
# ============================================================

LINE_NO=$((SLURM_ARRAY_TASK_ID + 1))

LINE=$(sed -n "${LINE_NO}p" "${REGION_FILE}")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: no region for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

read -r CHUNK REGION BUFFER <<< "${LINE}"

OUT=${OUTDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.chunk${CHUNK}.vcf.gz


echo "============================================================"
echo "Chr        : ${CHR2}"
echo "Density    : ${DENSITY}"
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


"${IMPUTE5}" \
    --h "${REF}" \
    --g "${TARGET}" \
    --m "${MAP}" \
    --o "${OUT}" \
    --r "${REGION}" \
    --buffer-region "${BUFFER}" \
    --threads ${SLURM_CPUS_PER_TASK}


bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
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
