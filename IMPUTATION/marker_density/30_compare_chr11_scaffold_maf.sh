#!/bin/bash

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
cd "${WORKDIR}"

OUTDIR=${WORKDIR}/Chr11/maf_comparison
mkdir -p "${OUTDIR}"

declare -A INPUTS

INPUTS[Actual_array]=${WORKDIR}/Chr11/array_real/Chr11.ARRAY50.conformed.vcf.gz
INPUTS[Random_rep1]=${WORKDIR}/Chr11/random1493/rep1/Chr11.random1493.rep1.target.unphased.bcf
INPUTS[Random_rep2]=${WORKDIR}/Chr11/random1493/rep2/Chr11.random1493.rep2.target.unphased.bcf
INPUTS[Random_rep3]=${WORKDIR}/Chr11/random1493/rep3/Chr11.random1493.rep3.target.unphased.bcf

for NAME in Actual_array Random_rep1 Random_rep2 Random_rep3
do
    IN=${INPUTS[$NAME]}
    TAGGED=${OUTDIR}/${NAME}.tags.bcf
    TSV=${OUTDIR}/${NAME}.maf.tsv

    echo "===== ${NAME} ====="

    if [[ ! -f "${IN}" ]]; then
        echo "ERROR: missing ${IN}" >&2
        exit 1
    fi

    NS=$(bcftools query -l "${IN}" | wc -l)
    NV=$(bcftools index -n "${IN}")

    echo "samples=${NS}"
    echo "variants=${NV}"

    if [[ "${NS}" -ne 50 ]]; then
        echo "ERROR: ${NAME} does not contain 50 samples" >&2
        exit 1
    fi

    if [[ "${NV}" -ne 1493 ]]; then
        echo "ERROR: ${NAME} does not contain 1493 markers" >&2
        exit 1
    fi

    # Recalculate AC, AN, AF from the 50 validation individuals
    bcftools +fill-tags \
        "${IN}" \
        -Ob \
        -o "${TAGGED}" \
        -- -t AC,AN,AF

    bcftools index -f "${TAGGED}"

    # Export AF and convert to MAF = min(AF, 1-AF)
    bcftools query \
        -f '%CHROM\t%POS\t%AC\t%AN\t%AF\n' \
        "${TAGGED}" \
    | awk 'BEGIN{OFS="\t"; print "CHROM","POS","AC","AN","AF","MAF"}
           {
             af=$5;
             if (af=="." || $4==0) {
                 maf="NA"
             } else {
                 maf=(af<=0.5 ? af : 1-af)
             }
             print $1,$2,$3,$4,$5,maf
           }' \
    > "${TSV}"

done

echo "Finished."
