#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=eval_r1493
#SBATCH -t 00-16:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/eval_r1493_%j.out
#SBATCH --error=logs/eval_r1493_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

REP=$1

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

CONDDIR=${WORKDIR}/Chr11/random1493/rep${REP}
REGION_FILE=${WORKDIR}/Chr11/Chr11.impute_regions.tsv
INDIR=${CONDDIR}/imputed_chunks

TRUTH=${WORKDIR}/Chr11/base/Chr11.WGS50.raw.vcf.gz
SCAFFOLD=${CONDDIR}/Chr11.random1493.rep${REP}.sites.tsv

MERGED=${CONDDIR}/Chr11.random1493.rep${REP}.imputed.vcf.gz
ACCURACY=${CONDDIR}/Chr11.random1493.rep${REP}.imputed_only.accuracy.tsv

EVAL=${WORKDIR}/11_evaluate_imputed_only.py

echo "============================================================"
echo "Chr        : 11"
echo "Markers    : 1493"
echo "Replicate  : ${REP}"
echo "============================================================"

for f in "${REGION_FILE}" "${TRUTH}" "${SCAFFOLD}" "${EVAL}"
do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: missing ${f}"
        exit 1
    fi
done

N_SCAFFOLD=$(wc -l < "${SCAFFOLD}")

if [[ "${N_SCAFFOLD}" -ne 1493 ]]; then
    echo "ERROR: expected 1493 scaffold SNPs, found ${N_SCAFFOLD}."
    exit 1
fi

FILES=()

while read -r CHUNK REGION BUFFER
do
    FILE=${INDIR}/Chr11.random1493.rep${REP}.chunk${CHUNK}.vcf.gz

    if [[ ! -f "${FILE}" ]]; then
        echo "ERROR: missing chunk ${FILE}"
        exit 1
    fi

    FILES+=("${FILE}")

done < "${REGION_FILE}"

echo "Number of chunks: ${#FILES[@]}"

rm -f \
    "${MERGED}" \
    "${MERGED}.csi" \
    "${MERGED}.tbi"

bcftools concat \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "${MERGED}" \
    "${FILES[@]}"

bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f "${MERGED}"

NSAMPLE=$(bcftools query -l "${MERGED}" | wc -l)
NVAR=$(bcftools index -n "${MERGED}")

echo "Merged samples : ${NSAMPLE}"
echo "Merged variants: ${NVAR}"

if [[ "${NSAMPLE}" -ne 50 ]]; then
    echo "ERROR: merged file does not contain 50 samples."
    exit 1
fi

python3 "${EVAL}" \
    --truth "${TRUTH}" \
    --imputed "${MERGED}" \
    --exclude-sites "${SCAFFOLD}" \
    --output "${ACCURACY}" \
    --progress 100000

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

if [[ "${TOTAL_SCAFFOLD}" -ne 1493 ]]; then
    echo "ERROR: evaluator reports unexpected scaffold count."
    exit 1
fi

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
echo "FINAL ACCURACY — RANDOM1493 REP${REP}"
echo "============================================"

grep -E \
'Pooled_empirical_R2|Mean_per_variant_dosage_R2|Mean_GP_weighted_concordance|Hard_genotype_concordance|Non_reference_concordance' \
"${ACCURACY}"

echo
echo "Finished: $(date)"
