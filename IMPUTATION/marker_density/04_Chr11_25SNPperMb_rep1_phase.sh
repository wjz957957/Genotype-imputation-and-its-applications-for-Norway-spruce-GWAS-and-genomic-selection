#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=Chr11_D25_phase
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-02:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/Chr11_D25_phase_%j.out
#SBATCH --error=logs/Chr11_D25_phase_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20
ml SHAPEIT/v5.1.1

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

INPUT=${WORKDIR}/Chr11/D25/rep1/Chr11.D25.rep1.target.unphased.bcf

OUTPUT=${WORKDIR}/Chr11/D25/rep1/Chr11.D25.rep1.target.phased.bcf

MAP=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map/LD.GenMap.chr11.txt

echo "Input : ${INPUT}"
echo "Map   : ${MAP}"
echo "Output: ${OUTPUT}"

SHAPEIT5_phase_common \
    -T ${SLURM_CPUS_PER_TASK} \
    -I "${INPUT}" \
    -O "${OUTPUT}" \
    --region 11 \
    -M "${MAP}" \
    --seed 1001

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${OUTPUT}"

echo
echo "========== PHASING QC =========="

echo -n "Samples: "
bcftools query -l "${OUTPUT}" | wc -l

echo -n "Variants: "
bcftools index -n "${OUTPUT}"

echo
echo "Example phased GT:"
bcftools query -f '[%GT\t]\n' "${OUTPUT}" \
    | head -5 || true

echo
echo "Finished: $(date)"
