#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=WGS50_raw
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-02:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/WGS50_raw_%A_%a.out
#SBATCH --error=logs/WGS50_raw_%A_%a.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

DATADIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map

IDFILE=${WORKDIR}/WGS_50.ids

CHRS=(01 05 06 11)
chr=${CHRS[$SLURM_ARRAY_TASK_ID]}

INPUT=${DATADIR}/Chr${chr}_Miss_missing0.2_maf0.05_hwe2.5e-10_info0.8_plink2dosage.vcf.gz
OUTDIR=${WORKDIR}/Chr${chr}/base
OUTPUT=${OUTDIR}/Chr${chr}.WGS50.raw.vcf.gz

mkdir -p "${OUTDIR}"

echo "================================================="
echo "Chromosome : Chr${chr}"
echo "Input      : ${INPUT}"
echo "Output     : ${OUTPUT}"
echo "Start      : $(date)"
echo "================================================="

# ---------- Check input ----------
if [[ ! -f "${INPUT}" ]]; then
    echo "ERROR: input file not found:"
    echo "${INPUT}"
    exit 1
fi

# ---------- Check that all 50 IDs exist ----------
N_MATCH=$(bcftools query -l "${INPUT}" | grep -Fxf "${IDFILE}" | wc -l)

echo "Matched validation samples: ${N_MATCH}"

if [[ "${N_MATCH}" -ne 50 ]]; then
    echo "ERROR: expected 50 WGS validation individuals."
    exit 1
fi

# ---------- Extract 50 validation individuals ----------
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -S "${IDFILE}" \
    -Oz \
    -o "${OUTPUT}" \
    "${INPUT}"

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${OUTPUT}"

# ---------- QC ----------
N_SAMPLE=$(bcftools query -l "${OUTPUT}" | wc -l)
N_VARIANT=$(bcftools index -n "${OUTPUT}")

echo
echo "========== QC =========="
echo "Samples  : ${N_SAMPLE}"
echo "Variants : ${N_VARIANT}"

if [[ "${N_SAMPLE}" -ne 50 ]]; then
    echo "ERROR: output does not contain 50 samples."
    exit 1
fi

echo
echo "Example genotypes:"
bcftools query -f '[%GT\t]\n' "${OUTPUT}" \
    | head -5 || true

echo
echo "Finished: $(date)"
