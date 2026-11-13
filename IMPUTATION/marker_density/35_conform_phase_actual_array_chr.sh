#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=conf_phase
#SBATCH -t 00-06:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/conf_phase_%j.out
#SBATCH --error=logs/conf_phase_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20
ml SHAPEIT/v5.1.1

CHR2=$1
CHRNUM=$((10#${CHR2}))

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
MAPDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map

CONFORM=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/conform-gt.24May16.cee.jar

cd "${WORKDIR}"

DIR=${WORKDIR}/Chr${CHR2}/array_real

ARRAY=${DIR}/Chr${CHR2}.ARRAY50.WGSid.unphased.vcf.gz
REF=${DIR}/Chr${CHR2}.REF997.vcf.gz
MAP=${MAPDIR}/LD.GenMap.chr${CHRNUM}.txt

PREFIX=${DIR}/Chr${CHR2}.ARRAY50.conformed
CONFORMED=${PREFIX}.vcf.gz
TAGGED=${DIR}/Chr${CHR2}.ARRAY50.conformed.tags.vcf.gz
PHASED=${DIR}/Chr${CHR2}.ARRAY50.conformed.phased.bcf
SITES=${DIR}/Chr${CHR2}.ARRAY50.sites.tsv
QC=${DIR}/Chr${CHR2}.array_conform_QC.txt

PRE_N=$(bcftools index -n "${ARRAY}")

echo "Pre-conform markers: ${PRE_N}"

rm -f "${CONFORMED}" "${CONFORMED}.csi" "${TAGGED}" "${TAGGED}.csi" \
      "${PHASED}" "${PHASED}.csi"

java -Xmx28g \
    -jar "${CONFORM}" \
    gt="${ARRAY}" \
    ref="${REF}" \
    match=POS \
    chrom="${CHRNUM}" \
    out="${PREFIX}"

bcftools index -c -f "${CONFORMED}"

POST_N=$(bcftools index -n "${CONFORMED}")
POST_NS=$(bcftools query -l "${CONFORMED}" | wc -l)

echo "Post-conform markers: ${POST_N}"

if [[ "${POST_NS}" -ne 50 || "${POST_N}" -le 0 ]]; then
    echo "ERROR: conform QC failed"
    exit 1
fi

# Scaffold positions used later for imputed-only evaluation
bcftools query \
    -f '%CHROM\t%POS\n' \
    "${CONFORMED}" \
    > "${SITES}"

if [[ $(wc -l < "${SITES}") -ne "${POST_N}" ]]; then
    echo "ERROR: site-list count mismatch"
    exit 1
fi

# SHAPEIT5 requires AC
bcftools +fill-tags \
    "${CONFORMED}" \
    -Oz \
    -o "${TAGGED}" \
    -- -t AC,AN,AF

bcftools index -c -f "${TAGGED}"

SHAPEIT5_phase_common \
    -T "${SLURM_CPUS_PER_TASK}" \
    -I "${TAGGED}" \
    -O "${PHASED}" \
    --region "${CHRNUM}" \
    -M "${MAP}" \
    --seed "$((20260000 + CHRNUM))"

bcftools index -f "${PHASED}"

PHASED_N=$(bcftools index -n "${PHASED}")
PHASED_NS=$(bcftools query -l "${PHASED}" | wc -l)

if [[ "${PHASED_N}" -ne "${POST_N}" || "${PHASED_NS}" -ne 50 ]]; then
    echo "ERROR: SHAPEIT output QC failed"
    exit 1
fi

{
    echo -e "Metric\tValue"
    echo -e "Chromosome\t${CHR2}"
    echo -e "Array_markers_before_conform\t${PRE_N}"
    echo -e "Usable_array_markers_after_conform\t${POST_N}"
    echo -e "Markers_removed_by_conform\t$((PRE_N-POST_N))"
    echo -e "Samples\t${PHASED_NS}"
} > "${QC}"

echo "Chr${CHR2}: usable actual-array markers = ${POST_N}"
echo "Finished: $(date)"
