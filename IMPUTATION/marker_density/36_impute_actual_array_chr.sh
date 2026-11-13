#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=arr_imp
#SBATCH -t 00-04:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/arr_imp_%A_%a.out
#SBATCH --error=logs/arr_imp_%A_%a.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

CHR2=$1
CHRNUM=$((10#${CHR2}))

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
MAPDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map
IMPUTE5=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/impute5_v1.2.1_x86_64_static

REF=${WORKDIR}/Chr${CHR2}/base/Chr${CHR2}.REF997.bcf
TARGET=${WORKDIR}/Chr${CHR2}/array_real/Chr${CHR2}.ARRAY50.conformed.phased.bcf
MAP=${MAPDIR}/LD.GenMap.chr${CHRNUM}.txt
REGION_FILE=${WORKDIR}/Chr${CHR2}/Chr${CHR2}.impute_regions.tsv

OUTDIR=${WORKDIR}/Chr${CHR2}/array_real/imputed_chunks
mkdir -p "${OUTDIR}"

LINE_NO=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${LINE_NO}p" "${REGION_FILE}")

[[ -n "${LINE}" ]] || { echo "ERROR no region"; exit 1; }

read -r CHUNK REGION BUFFER <<< "${LINE}"

OUT=${OUTDIR}/Chr${CHR2}.ARRAY50.chunk${CHUNK}.vcf.gz

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

NS=$(bcftools query -l "${OUT}" | wc -l)

[[ "${NS}" -eq 50 ]] || { echo "ERROR sample count"; exit 1; }

echo "Chr${CHR2} chunk ${CHUNK} completed"
