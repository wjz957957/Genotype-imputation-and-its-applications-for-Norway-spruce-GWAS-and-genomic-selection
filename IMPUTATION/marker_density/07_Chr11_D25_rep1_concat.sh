#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=Chr11_D25_concat
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-01:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/Chr11_D25_concat_%j.out
#SBATCH --error=logs/Chr11_D25_concat_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

INDIR=${WORKDIR}/Chr11/D25/rep1/imputed_chunks
OUTDIR=${WORKDIR}/Chr11/D25/rep1

OUT=${OUTDIR}/Chr11.D25.rep1.imputed.vcf.gz

FILES=(
"${INDIR}/Chr11.D25.rep1.chunk0.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk1.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk2.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk3.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk4.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk5.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk6.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk7.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk8.vcf.gz"
"${INDIR}/Chr11.D25.rep1.chunk9.vcf.gz"
)

echo "Checking input files..."

for f in "${FILES[@]}"
do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: missing ${f}"
        exit 1
    fi
done

echo "Concatenating..."

bcftools concat \
    --threads ${SLURM_CPUS_PER_TASK} \
    -Oz \
    -o "${OUT}" \
    "${FILES[@]}"

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${OUT}"

echo
echo "========== QC =========="

echo -n "Samples: "
bcftools query -l "${OUT}" | wc -l

echo -n "Variants: "
bcftools index -n "${OUT}"

echo
echo "First variants:"
bcftools query -f '%CHROM\t%POS\n' "${OUT}" | head -3 || true

echo
echo "Last variants:"
bcftools query -f '%CHROM\t%POS\n' "${OUT}" | tail -3

echo
echo "Finished: $(date)"
