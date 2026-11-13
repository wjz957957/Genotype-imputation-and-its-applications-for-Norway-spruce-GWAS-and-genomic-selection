#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=concat_eval
#SBATCH -t 00-16:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/final_%j.out
#SBATCH --error=logs/final_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20


CHR2=$1
DENSITY=$2
REP=$3


WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

CONDDIR=${WORKDIR}/Chr${CHR2}/D${DENSITY}/rep${REP}

REGION_FILE=${WORKDIR}/Chr${CHR2}/Chr${CHR2}.impute_regions.tsv

INDIR=${CONDDIR}/imputed_chunks

TRUTH=${WORKDIR}/Chr${CHR2}/base/Chr${CHR2}.WGS50.raw.vcf.gz

SCAFFOLD=${CONDDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.sites.tsv

MERGED=${CONDDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.imputed.vcf.gz

ACCURACY=${CONDDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.imputed_only.accuracy.tsv

EVAL=${WORKDIR}/11_evaluate_imputed_only.py


echo "============================================================"
echo "Chr       : ${CHR2}"
echo "Density   : ${DENSITY}"
echo "Replicate : ${REP}"
echo "============================================================"


# ============================================================
# Build ordered chunk list from region TSV
# ============================================================

FILES=()

while read -r CHUNK REGION BUFFER
do

    FILE=${INDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.chunk${CHUNK}.vcf.gz

    if [[ ! -f "${FILE}" ]]; then
        echo "ERROR: missing chunk:"
        echo "${FILE}"
        exit 1
    fi

    FILES+=("${FILE}")

done < "${REGION_FILE}"


echo "Number of chunks: ${#FILES[@]}"


# ============================================================
# CONCAT
# ============================================================

bcftools concat \
    --threads ${SLURM_CPUS_PER_TASK} \
    -Oz \
    -o "${MERGED}" \
    "${FILES[@]}"


bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${MERGED}"


NSAMPLE=$(bcftools query -l "${MERGED}" | wc -l)
NVAR=$(bcftools index -n "${MERGED}")


echo "Merged samples : ${NSAMPLE}"
echo "Merged variants: ${NVAR}"


if [[ "${NSAMPLE}" -ne 50 ]]; then
    echo "ERROR: merged file does not contain 50 samples."
    exit 1
fi


# ============================================================
# Accuracy — IMPUTED-ONLY
# ============================================================

python3 "${EVAL}" \
    --truth "${TRUTH}" \
    --imputed "${MERGED}" \
    --exclude-sites "${SCAFFOLD}" \
    --output "${ACCURACY}" \
    --progress 100000


# ============================================================
# Critical QC
# ============================================================

TOTAL_SCAFFOLD=$(awk \
    '$1=="Input_scaffold_SNPs_total" {print $2}' \
    "${ACCURACY}")

EXCLUDED_SCAFFOLD=$(awk \
    '$1=="Input_scaffold_SNPs_excluded" {print $2}' \
    "${ACCURACY}")

ALLELE_MISMATCH=$(awk \
    '$1=="Allele_mismatches" {print $2}' \
    "${ACCURACY}")


echo
echo "Input scaffold SNPs : ${TOTAL_SCAFFOLD}"
echo "Scaffold excluded   : ${EXCLUDED_SCAFFOLD}"
echo "Allele mismatches   : ${ALLELE_MISMATCH}"


if [[ "${TOTAL_SCAFFOLD}" -ne "${EXCLUDED_SCAFFOLD}" ]]; then
    echo "ERROR: not all scaffold SNPs were excluded."
    exit 1
fi


if [[ "${ALLELE_MISMATCH}" -ne 0 ]]; then
    echo "ERROR: allele mismatch detected."
    exit 1
fi


echo
echo "============================================"
echo "FINAL ACCURACY"
echo "============================================"

grep -E \
'Pooled_empirical_R2|Mean_per_variant_dosage_R2|Mean_GP_weighted_concordance|Hard_genotype_concordance|Non_reference_concordance' \
"${ACCURACY}"

echo
echo "Finished: $(date)"
