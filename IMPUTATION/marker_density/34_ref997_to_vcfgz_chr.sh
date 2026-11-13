#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=ref_vcfgz
#SBATCH -t 00-12:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/ref_vcfgz_%j.out
#SBATCH --error=logs/ref_vcfgz_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

CHR2=$1

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
cd "${WORKDIR}"

mkdir -p logs
mkdir -p Chr${CHR2}/array_real

INPUT=${WORKDIR}/Chr${CHR2}/base/Chr${CHR2}.REF997.bcf
OUTPUT=${WORKDIR}/Chr${CHR2}/array_real/Chr${CHR2}.REF997.vcf.gz
QC=${WORKDIR}/Chr${CHR2}/array_real/Chr${CHR2}.REF997.vcfgz_QC.txt

IN_SAMPLES=$(bcftools query -l "${INPUT}" | wc -l)
IN_VAR=$(bcftools index -n "${INPUT}")

if [[ "${IN_SAMPLES}" -ne 997 ]]; then
    echo "ERROR: REF does not contain 997 samples"
    exit 1
fi

rm -f "${OUTPUT}" "${OUTPUT}.csi" "${OUTPUT}.tbi"

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "${OUTPUT}" \
    "${INPUT}"

bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -c \
    "${OUTPUT}"

OUT_SAMPLES=$(bcftools query -l "${OUTPUT}" | wc -l)
OUT_VAR=$(bcftools index -n "${OUTPUT}")

if [[ "${IN_SAMPLES}" -ne "${OUT_SAMPLES}" ]] || \
   [[ "${IN_VAR}" -ne "${OUT_VAR}" ]]
then
    echo "ERROR: BCF→VCF conversion QC failed"
    exit 1
fi

{
    echo -e "Metric\tValue"
    echo -e "Chromosome\t${CHR2}"
    echo -e "Input_samples\t${IN_SAMPLES}"
    echo -e "Output_samples\t${OUT_SAMPLES}"
    echo -e "Input_variants\t${IN_VAR}"
    echo -e "Output_variants\t${OUT_VAR}"
} > "${QC}"

echo "Chr${CHR2} REF997 VCF completed: $(date)"
