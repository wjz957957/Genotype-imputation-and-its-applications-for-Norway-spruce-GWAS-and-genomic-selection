#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=arr_eval
#SBATCH -t 00-16:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/arr_eval_%j.out
#SBATCH --error=logs/arr_eval_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

CHR2=$1

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

DIR=${WORKDIR}/Chr${CHR2}/array_real
REGION_FILE=${WORKDIR}/Chr${CHR2}/Chr${CHR2}.impute_regions.tsv
INDIR=${DIR}/imputed_chunks
TRUTH=${WORKDIR}/Chr${CHR2}/base/Chr${CHR2}.WGS50.raw.vcf.gz
SCAFFOLD=${DIR}/Chr${CHR2}.ARRAY50.sites.tsv
MERGED=${DIR}/Chr${CHR2}.ARRAY50.imputed.vcf.gz
ACCURACY=${DIR}/Chr${CHR2}.ARRAY50.imputed_only.accuracy.tsv
EVAL=${WORKDIR}/11_evaluate_imputed_only.py

FILES=()

while read -r CHUNK REGION BUFFER
do
    FILE=${INDIR}/Chr${CHR2}.ARRAY50.chunk${CHUNK}.vcf.gz
    [[ -f "${FILE}" ]] || { echo "ERROR missing ${FILE}"; exit 1; }
    FILES+=("${FILE}")
done < "${REGION_FILE}"

bcftools concat \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "${MERGED}" \
    "${FILES[@]}"

bcftools index -f "${MERGED}"

python3 "${EVAL}" \
    --truth "${TRUTH}" \
    --imputed "${MERGED}" \
    --exclude-sites "${SCAFFOLD}" \
    --output "${ACCURACY}" \
    --progress 100000

TOTAL=$(awk '$1=="Input_scaffold_SNPs_total"{print $2}' "${ACCURACY}")
EXCL=$(awk '$1=="Input_scaffold_SNPs_excluded"{print $2}' "${ACCURACY}")
MISMATCH=$(awk '$1=="Allele_mismatches"{print $2}' "${ACCURACY}")

if [[ "${TOTAL}" -ne "${EXCL}" || "${MISMATCH}" -ne 0 ]]; then
    echo "ERROR: evaluation QC failed"
    exit 1
fi

echo "===== Chr${CHR2} ACTUAL ARRAY ====="

grep -E \
'Pooled_empirical_R2|Mean_per_variant_dosage_R2|Mean_GP_weighted_concordance|Hard_genotype_concordance|Non_reference_concordance' \
"${ACCURACY}"
