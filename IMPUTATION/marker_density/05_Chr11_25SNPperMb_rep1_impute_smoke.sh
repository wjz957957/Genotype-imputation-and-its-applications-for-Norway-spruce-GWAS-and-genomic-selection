#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=Chr11_D25_impute_test
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-04:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/Chr11_D25_impute_test_%j.out
#SBATCH --error=logs/Chr11_D25_impute_test_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

IMPUTE5=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/impute5_v1.2.1_x86_64_static

REF=${WORKDIR}/Chr11/base/Chr11.REF997.bcf

TARGET=${WORKDIR}/Chr11/D25/rep1/Chr11.D25.rep1.target.phased.bcf

MAP=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map/LD.GenMap.chr11.txt

OUTDIR=${WORKDIR}/Chr11/D25/rep1/impute_test
mkdir -p "${OUTDIR}"

OUT=${OUTDIR}/Chr11.D25.rep1.83036-100000000.vcf.gz

echo "======================================"
echo "Reference : ${REF}"
echo "Target    : ${TARGET}"
echo "Map       : ${MAP}"
echo "Output    : ${OUT}"
echo "======================================"

if [[ ! -x "${IMPUTE5}" ]]; then
    echo "ERROR: IMPUTE5 executable not found:"
    echo "${IMPUTE5}"
    exit 1
fi

"${IMPUTE5}" \
    --h "${REF}" \
    --g "${TARGET}" \
    --m "${MAP}" \
    --o "${OUT}" \
    --r 11:83036-100000000 \
    --buffer-region 11:1-100050000 \
    --threads 4

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${OUT}"

echo
echo "========== IMPUTATION QC =========="

echo -n "Imputed samples: "
bcftools query -l "${OUT}" | wc -l

echo -n "Imputed variants: "
bcftools index -n "${OUT}"

echo
echo "Finished: $(date)"
