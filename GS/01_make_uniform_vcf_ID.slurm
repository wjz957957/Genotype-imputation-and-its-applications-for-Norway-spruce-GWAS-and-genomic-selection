#!/bin/bash
#SBATCH -J make_vcf_ID
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o logs/01_make_vcf_ID_%j.out
#SBATCH -e logs/01_make_vcf_ID_%j.err


set -euo pipefail


echo "===================================="
echo "Step 2.1 Create uniform VCF IDs"
echo "Start: $(date)"
echo "===================================="


module load bioinfo-tools
ml bcftools/1.20


WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

INPUT=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/Imputed_array_combine.vcf.gz

OUTDIR=${WORK}/01_reference

OUTPUT=${OUTDIR}/Imputed_array_combine.ID.vcf.gz


mkdir -p ${OUTDIR}


echo "Input:"
echo ${INPUT}

echo "Output:"
echo ${OUTPUT}


echo "Creating IDs..."

bcftools annotate \
    --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' \
    -Oz \
    -o ${OUTPUT} \
    ${INPUT}


echo "Indexing..."

bcftools index \
    --csi \
    ${OUTPUT}


echo "Checking output..."

echo "Number of samples:"
bcftools query -l ${OUTPUT} | wc -l


echo "First variants:"
bcftools view -H ${OUTPUT} | head -5


echo "===================================="
echo "DONE"
echo "Finish: $(date)"
echo "===================================="
