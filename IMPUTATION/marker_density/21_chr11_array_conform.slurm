#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J Chr11_array_conform
#SBATCH -o logs/Chr11_array_conform_%j.out
#SBATCH -e logs/Chr11_array_conform_%j.err

set -euo pipefail

module load bioinfo-tools
module load bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
cd "${WORKDIR}"

mkdir -p logs
mkdir -p Chr11/array_real

ARRAY=${WORKDIR}/Chr11/array_real/Chr11.ARRAY50.WGSid.unphased.vcf.gz
REF=${WORKDIR}/Chr11/array_real/Chr11.REF997.vcf.gz

CONFORM=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/conform-gt.24May16.cee.jar

OUTPREFIX=${WORKDIR}/Chr11/array_real/Chr11.ARRAY50.conformed
OUT=${OUTPREFIX}.vcf.gz

WGS_IDS=${WORKDIR}/WGS_50.ids
QC=${WORKDIR}/Chr11/array_real/Chr11.array_conform_QC.txt


echo "============================================================"
echo "Chr11 real-array conform-gt"
echo "Job ID : ${SLURM_JOB_ID}"
echo "Start  : $(date)"
echo "============================================================"


# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

for F in "${ARRAY}" "${REF}" "${CONFORM}" "${WGS_IDS}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: missing input: ${F}" >&2
        exit 1
    fi
done

ARRAY_SAMPLES=$(bcftools query -l "${ARRAY}" | wc -l)
ARRAY_VARIANTS=$(bcftools index -n "${ARRAY}")

REF_SAMPLES=$(bcftools query -l "${REF}" | wc -l)
REF_VARIANTS=$(bcftools index -n "${REF}")

echo
echo "Array samples  : ${ARRAY_SAMPLES}"
echo "Array variants : ${ARRAY_VARIANTS}"
echo "REF samples    : ${REF_SAMPLES}"
echo "REF variants   : ${REF_VARIANTS}"

if [[ "${ARRAY_SAMPLES}" -ne 50 ]]; then
    echo "ERROR: expected 50 array samples" >&2
    exit 1
fi

if [[ "${REF_SAMPLES}" -ne 997 ]]; then
    echo "ERROR: expected 997 REF samples" >&2
    exit 1
fi


# ------------------------------------------------------------
# 2. Sample-name check
# ------------------------------------------------------------

bcftools query -l "${ARRAY}" \
    > Chr11/array_real/Chr11.array50.preconform.samples.txt

if ! diff -q \
    "${WGS_IDS}" \
    Chr11/array_real/Chr11.array50.preconform.samples.txt \
    >/dev/null
then
    echo "ERROR: array sample IDs/order do not match WGS_50.ids" >&2
    diff -u \
        "${WGS_IDS}" \
        Chr11/array_real/Chr11.array50.preconform.samples.txt || true
    exit 1
fi

echo "Array sample ID/order check: PASS"


# ------------------------------------------------------------
# 3. Confirm array GT is unphased
# ------------------------------------------------------------

N_PHASED_BEFORE=$(bcftools query \
    -f '[%GT\n]' \
    "${ARRAY}" \
    | grep -c '|' || true)

echo "Phased GT before conform: ${N_PHASED_BEFORE}"

if [[ "${N_PHASED_BEFORE}" -ne 0 ]]; then
    echo "ERROR: array target still contains phased GT" >&2
    exit 1
fi


# ------------------------------------------------------------
# 4. Remove previous conform output
# ------------------------------------------------------------

rm -f \
    "${OUT}" \
    "${OUT}.csi" \
    "${OUT}.tbi"


# ------------------------------------------------------------
# 5. Run conform-gt
# ------------------------------------------------------------

echo
echo "===== Running conform-gt ====="
echo "Start conform: $(date)"

java -Xmx28g \
    -jar "${CONFORM}" \
    gt="${ARRAY}" \
    ref="${REF}" \
    match=POS \
    chrom=11 \
    out="${OUTPREFIX}"

echo "End conform: $(date)"


# ------------------------------------------------------------
# 6. Check output
# ------------------------------------------------------------

echo
echo "===== Checking conform output ====="

if [[ ! -s "${OUT}" ]]; then
    echo "ERROR: expected conform output not found:"
    echo "${OUT}"
    echo
    echo "Files matching output prefix:"
    ls -lh "${OUTPREFIX}"* || true
    exit 1
fi

bcftools index \
    --threads 2 \
    -c \
    "${OUT}"

CONFORM_SAMPLES=$(bcftools query -l "${OUT}" | wc -l)
CONFORM_VARIANTS=$(bcftools index -n "${OUT}")

echo "Conformed samples  : ${CONFORM_SAMPLES}"
echo "Conformed variants : ${CONFORM_VARIANTS}"


# ------------------------------------------------------------
# 7. Sample-name check after conform
# ------------------------------------------------------------

bcftools query -l "${OUT}" \
    > Chr11/array_real/Chr11.array50.postconform.samples.txt

if ! diff -q \
    "${WGS_IDS}" \
    Chr11/array_real/Chr11.array50.postconform.samples.txt \
    >/dev/null
then
    echo "ERROR: sample IDs/order changed after conform" >&2
    diff -u \
        "${WGS_IDS}" \
        Chr11/array_real/Chr11.array50.postconform.samples.txt || true
    exit 1
fi

echo "Post-conform sample ID/order check: PASS"


# ------------------------------------------------------------
# 8. Chromosome check
# ------------------------------------------------------------

CHR_VALUES=$(bcftools query \
    -f '%CHROM\n' \
    "${OUT}" \
    | sort -u \
    | paste -sd, -)

echo "Chromosome values: ${CHR_VALUES}"

if [[ "${CHR_VALUES}" != "11" ]]; then
    echo "ERROR: unexpected chromosome value(s)" >&2
    exit 1
fi


# ------------------------------------------------------------
# 9. Phase-state check after conform
# ------------------------------------------------------------

N_PHASED_AFTER=$(bcftools query \
    -f '[%GT\n]' \
    "${OUT}" \
    | grep -c '|' || true)

echo "Phased GT after conform: ${N_PHASED_AFTER}"

if [[ "${N_PHASED_AFTER}" -ne 0 ]]; then
    echo "WARNING: conform output contains phased GT"
fi


# ------------------------------------------------------------
# 10. Marker-retention statistics
# ------------------------------------------------------------

REMOVED=$((ARRAY_VARIANTS - CONFORM_VARIANTS))

RETENTION=$(awk \
    -v a="${ARRAY_VARIANTS}" \
    -v c="${CONFORM_VARIANTS}" \
    'BEGIN {
        if (a>0) printf "%.4f", 100*c/a;
        else print "NA"
    }')

echo
echo "Pre-conform markers : ${ARRAY_VARIANTS}"
echo "Post-conform markers: ${CONFORM_VARIANTS}"
echo "Removed markers     : ${REMOVED}"
echo "Retention (%)       : ${RETENTION}"


# ------------------------------------------------------------
# 11. Write QC summary
# ------------------------------------------------------------

{
    echo -e "Metric\tValue"
    echo -e "Chromosome\t11"

    echo -e "Array_samples_before_conform\t${ARRAY_SAMPLES}"
    echo -e "Array_markers_before_conform\t${ARRAY_VARIANTS}"

    echo -e "REF997_samples\t${REF_SAMPLES}"
    echo -e "REF997_variants\t${REF_VARIANTS}"

    echo -e "Array_samples_after_conform\t${CONFORM_SAMPLES}"
    echo -e "Usable_array_markers_after_conform\t${CONFORM_VARIANTS}"

    echo -e "Markers_removed_by_conform\t${REMOVED}"
    echo -e "Marker_retention_percent\t${RETENTION}"

    echo -e "Phased_GT_before_conform\t${N_PHASED_BEFORE}"
    echo -e "Phased_GT_after_conform\t${N_PHASED_AFTER}"

    echo -e "Conformed_VCF\t${OUT}"
} > "${QC}"


echo
echo "============================================================"
echo "Chr11 conform-gt completed successfully"
echo
echo "Before conform : ${ARRAY_VARIANTS}"
echo "After conform  : ${CONFORM_VARIANTS}"
echo "Retention (%)  : ${RETENTION}"
echo
echo "QC:"
echo "${QC}"
echo
echo "Finished: $(date)"
echo "============================================================"
