#!/bin/bash

set -euo pipefail

# ============================================================
# Usage:
#
# bash generate_impute_regions.sh CHR REF_BCF CHR_LENGTH
#
# Example:
#
# bash generate_impute_regions.sh \
# 11 \
# Chr11/base/Chr11.REF997.bcf \
# 954579047
# ============================================================

CHR=$1
REF=$2
CHR_LENGTH=$3

CHUNK=100000000
BUFFER=50000
SEARCH_WINDOW=10000000


# ============================================================
# Check input
# ============================================================

if [[ ! -f "${REF}" ]]; then
    echo "ERROR: reference file not found:"
    echo "${REF}"
    exit 1
fi

if [[ ! -f "${REF}.csi" && ! -f "${REF}.tbi" ]]; then
    echo "ERROR: reference file is not indexed:"
    echo "${REF}"
    exit 1
fi


# ============================================================
# Find first panel SNP
# ============================================================

FIRST_POS=""
SEARCH_START=1

while [[ "${SEARCH_START}" -le "${CHR_LENGTH}" ]]
do

    SEARCH_END=$((SEARCH_START + SEARCH_WINDOW - 1))

    if [[ "${SEARCH_END}" -gt "${CHR_LENGTH}" ]]; then
        SEARCH_END=${CHR_LENGTH}
    fi

    # head closes pipe early, so temporarily disable pipefail
    set +o pipefail

    FIRST_POS=$(
        bcftools query \
        -r "${CHR}:${SEARCH_START}-${SEARCH_END}" \
        -f '%POS\n' \
        "${REF}" 2>/dev/null \
        | head -1
    )

    set -o pipefail

    if [[ -n "${FIRST_POS}" ]]; then
        break
    fi

    SEARCH_START=$((SEARCH_END + 1))
done


if [[ -z "${FIRST_POS}" ]]; then
    echo "ERROR: no SNP found on chromosome ${CHR}"
    exit 1
fi


# ============================================================
# Find last panel SNP
# ============================================================

LAST_POS=""
SEARCH_END=${CHR_LENGTH}

while [[ "${SEARCH_END}" -ge 1 ]]
do

    SEARCH_START=$((SEARCH_END - SEARCH_WINDOW + 1))

    if [[ "${SEARCH_START}" -lt 1 ]]; then
        SEARCH_START=1
    fi

    LAST_POS=$(
        bcftools query \
        -r "${CHR}:${SEARCH_START}-${SEARCH_END}" \
        -f '%POS\n' \
        "${REF}" 2>/dev/null \
        | tail -1
    )

    if [[ -n "${LAST_POS}" ]]; then
        break
    fi

    SEARCH_END=$((SEARCH_START - 1))
done


if [[ -z "${LAST_POS}" ]]; then
    echo "ERROR: no SNP found on chromosome ${CHR}"
    exit 1
fi


# ============================================================
# General information
# ============================================================

echo "============================================"
echo "Chromosome       : ${CHR}"
echo "Physical length  : ${CHR_LENGTH}"
echo "First panel SNP  : ${FIRST_POS}"
echo "Last panel SNP   : ${LAST_POS}"
echo "Chunk size       : ${CHUNK}"
echo "Buffer size      : ${BUFFER}"
echo "============================================"
echo

printf "Chunk\tRegion\tBuffer_region\n"


# ============================================================
# Generate regions
# ============================================================

REGION_START=${FIRST_POS}
INDEX=0

while [[ "${REGION_START}" -le "${LAST_POS}" ]]
do

    # --------------------------------------------------------
    # Find next 100-Mb boundary
    #
    # e.g.
    # 42,607       -> 100,000,000
    # 100,000,001  -> 200,000,000
    # --------------------------------------------------------

    REGION_END=$(( ((REGION_START - 1) / CHUNK + 1) * CHUNK ))


    # Last formal region ends at last reference-panel SNP
    if [[ "${REGION_END}" -gt "${LAST_POS}" ]]; then
        REGION_END=${LAST_POS}
    fi


    # ========================================================
    # Buffer start
    # ========================================================

    if [[ "${INDEX}" -eq 0 ]]; then

        # First buffer begins at chromosome coordinate 1
        BUFFER_START=1

    else

        BUFFER_START=$((REGION_START - BUFFER))

        if [[ "${BUFFER_START}" -lt 1 ]]; then
            BUFFER_START=1
        fi

    fi


    # ========================================================
    # Buffer end
    # ========================================================

    if [[ "${REGION_END}" -eq "${LAST_POS}" ]]; then

        # Last buffer reaches physical chromosome end
        BUFFER_END=${CHR_LENGTH}

    else

        BUFFER_END=$((REGION_END + BUFFER))

        if [[ "${BUFFER_END}" -gt "${CHR_LENGTH}" ]]; then
            BUFFER_END=${CHR_LENGTH}
        fi

    fi


    # ========================================================
    # Print
    # ========================================================

    printf "%02d\t%s:%d-%d\t%s:%d-%d\n" \
        "${INDEX}" \
        "${CHR}" "${REGION_START}" "${REGION_END}" \
        "${CHR}" "${BUFFER_START}" "${BUFFER_END}"


    # Next region
    REGION_START=$((REGION_END + 1))
    INDEX=$((INDEX + 1))

done
