#!/bin/bash

set -euo pipefail


TRUE="Chr11_MAF-347_with_MAF.vcf.gz"

NCHUNKS="${1:-8}"

OUT="NRC_regions.tsv"


# ------------------------------------------------------------
# Check bcftools
# ------------------------------------------------------------

if ! command -v bcftools >/dev/null 2>&1
then
    echo "ERROR: bcftools not found in PATH" >&2
    exit 1
fi


# ------------------------------------------------------------
# Check VCF index
# ------------------------------------------------------------

if [[ ! -f "${TRUE}.csi" && ! -f "${TRUE}.tbi" ]]
then
    echo "Index not found. Creating index for ${TRUE}..."
    bcftools index "${TRUE}"
fi


# ------------------------------------------------------------
# Get chromosome
# ------------------------------------------------------------

CHR=$(
    bcftools index -s "${TRUE}" |
    head -n 1 |
    cut -f1
)


if [[ -z "${CHR}" ]]
then
    echo "ERROR: could not determine chromosome name." >&2
    exit 1
fi


# ------------------------------------------------------------
# Obtain maximum observed position
#
# This is robust even when the index does not contain
# chromosome length information.
# ------------------------------------------------------------

echo "Determining maximum position..."

MAXPOS=$(
    bcftools query \
        -r "${CHR}" \
        -f '%POS\n' \
        "${TRUE}" |
    tail -n 1
)


if ! [[ "${MAXPOS}" =~ ^[0-9]+$ ]]
then
    echo "ERROR: invalid maximum position: ${MAXPOS}" >&2
    exit 1
fi


echo "Chromosome: ${CHR}"
echo "Maximum position: ${MAXPOS}"
echo "Number of chunks: ${NCHUNKS}"


# ------------------------------------------------------------
# Calculate chunk size
# Ceiling division:
#
# STEP = ceil(MAXPOS / NCHUNKS)
# ------------------------------------------------------------

STEP=$(( (MAXPOS + NCHUNKS - 1) / NCHUNKS ))


echo "Chunk size: ${STEP}"


# ------------------------------------------------------------
# Create output file
# ------------------------------------------------------------

: > "${OUT}"


for IDX in $(seq 1 "${NCHUNKS}")
do

    START=$(( (IDX - 1) * STEP + 1 ))

    END=$(( IDX * STEP ))


    if [[ "${END}" -gt "${MAXPOS}" ]]
    then
        END="${MAXPOS}"
    fi


    # Stop if start exceeds chromosome endpoint
    if [[ "${START}" -gt "${MAXPOS}" ]]
    then
        break
    fi


    printf "%d\t%s:%d-%d\n" \
        "${IDX}" \
        "${CHR}" \
        "${START}" \
        "${END}" \
        >> "${OUT}"

done


echo
echo "Created ${OUT}:"
cat "${OUT}"
