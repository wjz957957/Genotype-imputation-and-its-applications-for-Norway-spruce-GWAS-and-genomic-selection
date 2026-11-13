#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=prep_array
#SBATCH -t 00-02:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/prep_array_%j.out
#SBATCH --error=logs/prep_array_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

CHR2=$1

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
ARRAYDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/SNParray

cd "${WORKDIR}"

mkdir -p logs
mkdir -p Chr${CHR2}/array_real

ARRAY=${ARRAYDIR}/mod.chr${CHR2}.consistent.phased.bcf

ARRAY_IDS=${WORKDIR}/ARRAY_50.frommap.ids
RENAME=${WORKDIR}/ARRAY_to_WGS.rename.tsv
WGS_IDS=${WORKDIR}/WGS_50.ids

ORIG=${WORKDIR}/Chr${CHR2}/array_real/Chr${CHR2}.ARRAY50.original_phased.vcf.gz
RENAMED=${WORKDIR}/Chr${CHR2}/array_real/Chr${CHR2}.ARRAY50.WGSid.phased.vcf.gz
UNPHASED=${WORKDIR}/Chr${CHR2}/array_real/Chr${CHR2}.ARRAY50.WGSid.unphased.vcf.gz

for f in "${ARRAY}" "${ARRAY_IDS}" "${RENAME}" "${WGS_IDS}"
do
    [[ -s "${f}" ]] || { echo "ERROR missing ${f}"; exit 1; }
done

echo "===== Chr${CHR2}: extract 50 array samples ====="

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -S "${ARRAY_IDS}" \
    -Oz \
    -o "${ORIG}" \
    "${ARRAY}"

bcftools index -c -f "${ORIG}"

NS=$(bcftools query -l "${ORIG}" | wc -l)

if [[ "${NS}" -ne 50 ]]; then
    echo "ERROR: expected 50 samples, found ${NS}"
    exit 1
fi

echo "===== Rename array IDs → WGS IDs ====="

bcftools reheader \
    -s "${RENAME}" \
    -o "${RENAMED}" \
    "${ORIG}"

bcftools index -c -f "${RENAMED}"

bcftools query -l "${RENAMED}" > Chr${CHR2}/array_real/sample_check.txt

if ! diff -q "${WGS_IDS}" Chr${CHR2}/array_real/sample_check.txt >/dev/null
then
    echo "ERROR: sample ID/order mismatch"
    exit 1
fi

echo "===== Check FORMAT ====="

N_FORMAT=$(bcftools view -h "${RENAMED}" | grep -c '^##FORMAT=')

if [[ "${N_FORMAT}" -ne 1 ]]; then
    echo "ERROR: expected GT-only FORMAT, found ${N_FORMAT} FORMAT definitions"
    bcftools view -h "${RENAMED}" | grep '^##FORMAT='
    exit 1
fi

if ! bcftools view -h "${RENAMED}" | grep '^##FORMAT=<ID=GT' >/dev/null
then
    echo "ERROR: sole FORMAT is not GT"
    exit 1
fi

echo "===== Remove previous phasing ====="

bcftools view -Ov "${RENAMED}" \
| awk '
BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
    for(i=10;i<=NF;i++) gsub(/\|/,"/",$i)
    print
}' \
| bgzip -c > "${UNPHASED}"

bcftools index -c -f "${UNPHASED}"

NVAR=$(bcftools index -n "${UNPHASED}")
NS=$(bcftools query -l "${UNPHASED}" | wc -l)

N_PHASED=$(bcftools query -f '[%GT\n]' "${UNPHASED}" | grep -c '|' || true)

echo "Samples           : ${NS}"
echo "Pre-conform SNPs  : ${NVAR}"
echo "Phased GT remains : ${N_PHASED}"

if [[ "${NS}" -ne 50 || "${N_PHASED}" -ne 0 ]]; then
    echo "ERROR: unphase QC failed"
    exit 1
fi

echo "Finished Chr${CHR2}: $(date)"
