#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -J Chr11_REF997_vcfgz
#SBATCH -o logs/Chr11_REF997_vcfgz_%j.out
#SBATCH -e logs/Chr11_REF997_vcfgz_%j.err

set -euo pipefail

module load bioinfo-tools
module load bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

cd "${WORKDIR}"

mkdir -p logs
mkdir -p Chr11/array_real

INPUT=${WORKDIR}/Chr11/base/Chr11.REF997.bcf
OUTPUT=${WORKDIR}/Chr11/array_real/Chr11.REF997.vcf.gz
QC=${WORKDIR}/Chr11/array_real/Chr11.REF997.vcfgz_QC.txt


echo "============================================================"
echo "Convert Chr11 REF997 BCF to VCF.GZ"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================================"


# ============================================================
# 1. Check input
# ============================================================

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: input file not found:"
    echo "${INPUT}"
    exit 1
fi

echo
echo "Input:"
ls -lh "${INPUT}"

INPUT_SAMPLES=$(bcftools query -l "${INPUT}" | wc -l)
INPUT_VARIANTS=$(bcftools index -n "${INPUT}")

echo "Input samples  : ${INPUT_SAMPLES}"
echo "Input variants : ${INPUT_VARIANTS}"

if [[ "${INPUT_SAMPLES}" -ne 997 ]]; then
    echo "ERROR: expected 997 reference samples, found ${INPUT_SAMPLES}"
    exit 1
fi


# ============================================================
# 2. Remove incomplete previous output if present
# ============================================================

echo
echo "Removing any previous incomplete output..."

rm -f \
    "${OUTPUT}" \
    "${OUTPUT}.csi" \
    "${OUTPUT}.tbi"


# ============================================================
# 3. Convert BCF -> VCF.GZ
# ============================================================

echo
echo "===== Converting BCF to VCF.GZ ====="

bcftools view \
    --threads 4 \
    -Oz \
    -o "${OUTPUT}" \
    "${INPUT}"


# ============================================================
# 4. Basic file check
# ============================================================

if [[ ! -s "${OUTPUT}" ]]; then
    echo "ERROR: output VCF.GZ was not created"
    exit 1
fi

echo
echo "Output created:"
ls -lh "${OUTPUT}"


# ============================================================
# 5. Test compressed VCF integrity
# ============================================================

echo
echo "===== Testing compressed file integrity ====="

bcftools view \
    -h \
    "${OUTPUT}" \
    >/dev/null

echo "Header read: PASS"


# ============================================================
# 6. Create CSI index
# ============================================================

echo
echo "===== Creating CSI index ====="

bcftools index \
    --threads 4 \
    -c \
    "${OUTPUT}"

if [[ ! -s "${OUTPUT}.csi" ]]; then
    echo "ERROR: CSI index was not created"
    exit 1
fi


# ============================================================
# 7. Compare sample count and variant count
# ============================================================

echo
echo "===== Output QC ====="

OUTPUT_SAMPLES=$(bcftools query -l "${OUTPUT}" | wc -l)
OUTPUT_VARIANTS=$(bcftools index -n "${OUTPUT}")

echo "Output samples  : ${OUTPUT_SAMPLES}"
echo "Output variants : ${OUTPUT_VARIANTS}"

if [[ "${INPUT_SAMPLES}" -ne "${OUTPUT_SAMPLES}" ]]; then
    echo "ERROR: sample count changed during conversion"
    exit 1
fi

if [[ "${INPUT_VARIANTS}" -ne "${OUTPUT_VARIANTS}" ]]; then
    echo "ERROR: variant count changed during conversion"
    echo "Input : ${INPUT_VARIANTS}"
    echo "Output: ${OUTPUT_VARIANTS}"
    exit 1
fi


# ============================================================
# 8. Check chromosome names
# ============================================================

CHR_VALUES=$(bcftools query \
    -f '%CHROM\n' \
    "${OUTPUT}" \
    | sort -u \
    | paste -sd, -)

echo "Chromosome values: ${CHR_VALUES}"

if [[ "${CHR_VALUES}" != "11" ]]; then
    echo "WARNING: unexpected chromosome value(s): ${CHR_VALUES}"
fi


# ============================================================
# 9. Write QC summary
# ============================================================

{
    echo -e "Metric\tValue"
    echo -e "Input_BCF\t${INPUT}"
    echo -e "Output_VCFGZ\t${OUTPUT}"
    echo -e "Input_samples\t${INPUT_SAMPLES}"
    echo -e "Output_samples\t${OUTPUT_SAMPLES}"
    echo -e "Input_variants\t${INPUT_VARIANTS}"
    echo -e "Output_variants\t${OUTPUT_VARIANTS}"
    echo -e "Chromosome\t${CHR_VALUES}"
    echo -e "CSI_index\t${OUTPUT}.csi"
} > "${QC}"


echo
echo "============================================================"
echo "Chr11 REF997 conversion completed successfully"
echo
echo "VCF.GZ:"
echo "${OUTPUT}"
echo
echo "CSI:"
echo "${OUTPUT}.csi"
echo
echo "QC:"
echo "${QC}"
echo
echo "Finished: $(date)"
echo "============================================================"
