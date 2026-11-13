#!/bin/bash
#SBATCH -J cf_rrblup_geno
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-14%3
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/04_rrblup_geno_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/04_rrblup_geno_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

########################################
# 1. Paths and configuration
########################################

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"

SAMPLE_ROOT="${WORK}/Cross_family/cv_samples_v1"
MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
CHIP_MARKERS="${WORK}/01_reference/chip.marker"
OUTROOT="${WORK}/Cross_family/results/04_inner_rrblup_genotypes"

N_EXPECTED_SNP=43267
SCRIPT_VERSION=2

########################################
# 2. Software
########################################

module load bioinfo-tools
ml plink/1.90b4.9

########################################
# 3. Map task to outer/inner fold
########################################

TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

if [[ ! "${TASK}" =~ ^[0-9]+$ ]] || (( TASK < 0 || TASK > 14 )); then
    echo "ERROR: Task must be between 0 and 14" >&2
    exit 1
fi

OUTER=$(printf '%02d' "$((TASK / 3 + 1))")
INNER=$(printf '%02d' "$((TASK % 3 + 1))")

FOLD="${SAMPLE_ROOT}/outer_${OUTER}"

TRAIN="${FOLD}/inner_${INNER}/inner_train.samples"
VALID="${FOLD}/inner_${INNER}/inner_valid.samples"

OUTER_TRAIN="${FOLD}/outer_train.samples"
OUTER_TEST="${FOLD}/outer_test.samples"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}"

mkdir -p "${OUT}"

echo "========================================"
echo "Prepare RRBLUP genotype inputs"
echo "Script version: ${SCRIPT_VERSION}"
echo "Task:  ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

########################################
# 4. Validate input files
########################################

for EXT in bed bim fam; do
    [[ -s "${MASTER}.${EXT}" ]] || {
        echo "ERROR: Missing master ${EXT}" >&2
        exit 1
    }
done

for F in "${TRAIN}" "${VALID}" "${OUTER_TRAIN}" \
         "${OUTER_TEST}" "${CHIP_MARKERS}"; do
    [[ -s "${F}" ]] || {
        echo "ERROR: Missing input ${F}" >&2
        exit 1
    }
done

N_MASTER=$(wc -l < "${MASTER}.fam")

if (( N_MASTER != 904 )); then
    echo "ERROR: Master FAM must contain 904 samples" >&2
    exit 1
fi

N_MARKER_LINES=$(wc -l < "${CHIP_MARKERS}")
N_MARKER_UNIQUE=$(sort -u "${CHIP_MARKERS}" | wc -l)

if (( N_MARKER_LINES != N_EXPECTED_SNP ||
      N_MARKER_UNIQUE != N_EXPECTED_SNP )); then
    echo "ERROR: Chip marker count/uniqueness mismatch" >&2
    echo "Lines=${N_MARKER_LINES}" >&2
    echo "Unique=${N_MARKER_UNIQUE}" >&2
    exit 1
fi

########################################
# 5. Validate sample partitions
########################################

for F in "${TRAIN}" "${VALID}" "${OUTER_TRAIN}" "${OUTER_TEST}"; do
    awk '
    {
        sub(/\r$/, "", $0)

        if (NF != 1 || seen[$1]++) {
            print "ERROR: Invalid or duplicate sample ID:", $0 \
                > "/dev/stderr"
            bad=1
        }
    }
    END {
        if (NR == 0 || bad) exit 1
    }
    ' "${F}"
done

sort "${TRAIN}" > "${OUT}/train_ids.sorted"
sort "${VALID}" > "${OUT}/valid_ids.sorted"
sort "${OUTER_TRAIN}" > "${OUT}/outer_train_ids.sorted"
sort "${OUTER_TEST}" > "${OUT}/outer_test_ids.sorted"

if comm -12 "${OUT}/train_ids.sorted" \
            "${OUT}/valid_ids.sorted" | grep -q .; then
    echo "ERROR: Inner train/validation overlap" >&2
    exit 1
fi

cat "${OUT}/train_ids.sorted" "${OUT}/valid_ids.sorted" \
    | sort > "${OUT}/inner_union_ids.sorted"

if ! cmp -s "${OUT}/inner_union_ids.sorted" \
             "${OUT}/outer_train_ids.sorted"; then
    echo "ERROR: Inner train + validation do not equal outer train" >&2
    exit 1
fi

if comm -12 "${OUT}/train_ids.sorted" \
            "${OUT}/outer_test_ids.sorted" | grep -q .; then
    echo "ERROR: Outer-test leakage into training" >&2
    exit 1
fi

if comm -12 "${OUT}/valid_ids.sorted" \
            "${OUT}/outer_test_ids.sorted" | grep -q .; then
    echo "ERROR: Outer-test leakage into validation" >&2
    exit 1
fi

N_TRAIN=$(wc -l < "${TRAIN}")
N_VALID=$(wc -l < "${VALID}")

echo "Master samples:          ${N_MASTER}"
echo "Inner training samples:  ${N_TRAIN}"
echo "Inner validation samples:${N_VALID}"
echo "Expected chip SNPs:      ${N_EXPECTED_SNP}"
echo "PASS: Sample partitions and outer-test exclusion"

########################################
# 6. Prepare expected SNP IDs
########################################

EXPECTED_SNPS="${OUT}/expected_chip_snps.sorted"

sort -u "${CHIP_MARKERS}" > "${EXPECTED_SNPS}"

########################################
# 7. Validate an exported PLINK raw file
########################################

validate_raw() {
    local RAW="$1"
    local SAMPLE_FILE="$2"
    local LABEL="$3"

    [[ -s "${RAW}" ]] || {
        echo "ERROR: Missing ${LABEL}: ${RAW}" >&2
        return 1
    }

    local N_EXPECTED
    local N_OBSERVED
    local N_SNPS

    N_EXPECTED=$(wc -l < "${SAMPLE_FILE}")
    N_OBSERVED=$(( $(wc -l < "${RAW}") - 1 ))
    N_SNPS=$(awk 'NR == 1 {print NF-6; exit}' "${RAW}")

    if (( N_OBSERVED != N_EXPECTED )); then
        echo "ERROR: ${LABEL} sample count mismatch" >&2
        echo "Expected=${N_EXPECTED}; observed=${N_OBSERVED}" >&2
        return 1
    fi

    if (( N_SNPS != N_EXPECTED_SNP )); then
        echo "ERROR: ${LABEL} SNP count mismatch" >&2
        echo "Expected=${N_EXPECTED_SNP}; observed=${N_SNPS}" >&2
        return 1
    fi

    awk '
    NR > 1 && $1 != $2 {
        print "ERROR: FID/IID mismatch at row", NR > "/dev/stderr"
        bad=1
    }
    END {
        if (bad) exit 1
    }
    ' "${RAW}"

    local OBSERVED_IDS="${OUT}/${LABEL}_observed_ids.sorted"
    local OBSERVED_SNPS="${OUT}/${LABEL}_observed_snps.sorted"
    local EXPECTED_IDS="${OUT}/${LABEL}_expected_ids.sorted"

    awk 'NR > 1 {print $2}' "${RAW}" \
        | sort > "${OBSERVED_IDS}"

    sort "${SAMPLE_FILE}" > "${EXPECTED_IDS}"

    if ! cmp -s "${EXPECTED_IDS}" "${OBSERVED_IDS}"; then
        echo "ERROR: ${LABEL} sample IDs differ from expected" >&2
        diff -u "${EXPECTED_IDS}" "${OBSERVED_IDS}" | head -n 30 || true
        return 1
    fi

    # Recover SNP IDs from PLINK --recode A headers.
    # Example: SNP_ID_A1 -> SNP_ID.
    awk '
    NR == 1 {
        for (i=7; i<=NF; i++) {
            id=$i
            sub(/_[^_]+$/, "", id)
            print id
        }
    }
    ' "${RAW}" | sort > "${OBSERVED_SNPS}"

    if ! cmp -s "${EXPECTED_SNPS}" "${OBSERVED_SNPS}"; then
        echo "ERROR: ${LABEL} SNP IDs differ from chip.marker" >&2
        diff -u "${EXPECTED_SNPS}" "${OBSERVED_SNPS}" | head -n 30 || true
        return 1
    fi

    echo "PASS: ${LABEL}: samples=${N_OBSERVED}; SNPs=${N_SNPS}"
}

########################################
# 8. Validate train/validation SNP headers
########################################

validate_raw_headers() {
    local TRAIN_RAW="$1"
    local VALID_RAW="$2"

    local TRAIN_HEADERS="${OUT}/train_snp_headers.txt"
    local VALID_HEADERS="${OUT}/valid_snp_headers.txt"

    awk '
    NR == 1 {
        for (i=7; i<=NF; i++) print $i
        exit
    }
    ' "${TRAIN_RAW}" > "${TRAIN_HEADERS}"

    awk '
    NR == 1 {
        for (i=7; i<=NF; i++) print $i
        exit
    }
    ' "${VALID_RAW}" > "${VALID_HEADERS}"

    local N_TRAIN_HEADERS
    local N_VALID_HEADERS

    N_TRAIN_HEADERS=$(wc -l < "${TRAIN_HEADERS}")
    N_VALID_HEADERS=$(wc -l < "${VALID_HEADERS}")

    if (( N_TRAIN_HEADERS != N_EXPECTED_SNP ||
          N_VALID_HEADERS != N_EXPECTED_SNP )); then
        echo "ERROR: SNP header count mismatch" >&2
        echo "Training headers: ${N_TRAIN_HEADERS}" >&2
        echo "Validation headers: ${N_VALID_HEADERS}" >&2
        return 1
    fi

    if ! cmp -s "${TRAIN_HEADERS}" "${VALID_HEADERS}"; then
        echo "ERROR: Training and validation SNP headers differ" >&2
        echo "SNP order or counted alleles may differ." >&2

        diff -u "${TRAIN_HEADERS}" "${VALID_HEADERS}" \
            | head -n 40 || true

        return 1
    fi

    echo "PASS: All ${N_EXPECTED_SNP} SNP headers match exactly"
    echo "PASS: SNP order and counted-allele suffixes match"
}

########################################
# 9. Output prefixes
########################################

PREFIX_TRAIN="${OUT}/chip_all_train"
PREFIX_VALID="${OUT}/chip_all_valid"

# Temporary prefixes prevent unvalidated output from replacing
# the existing final .raw files.
TMP_TRAIN="${OUT}/chip_all_train.v2tmp"
TMP_VALID="${OUT}/chip_all_valid.v2tmp"

FLAG="${OUT}/complete.flag"

########################################
# 10. Reuse only verified version-2 results
########################################

if [[ -s "${FLAG}" &&
      -s "${PREFIX_TRAIN}.raw" &&
      -s "${PREFIX_VALID}.raw" ]] &&
   grep -qx "version=${SCRIPT_VERSION}" "${FLAG}"; then

    if validate_raw "${PREFIX_TRAIN}.raw" "${TRAIN}" train &&
       validate_raw "${PREFIX_VALID}.raw" "${VALID}" valid &&
       validate_raw_headers "${PREFIX_TRAIN}.raw" "${PREFIX_VALID}.raw"; then

        grep -qx "train_valid_snp_header_match=PASS" "${FLAG}" || {
            echo "ERROR: Completion flag lacks header-match status" >&2
            exit 1
        }

        echo "Previously completed version-2 files passed revalidation."
        exit 0
    fi

    echo "Existing output failed validation; regenerating."
fi

########################################
# 11. Prepare KEEP files
########################################

KEEP_TRAIN="${OUT}/train.keep"
KEEP_VALID="${OUT}/valid.keep"

awk '{print $1, $1}' "${TRAIN}" > "${KEEP_TRAIN}"
awk '{print $1, $1}' "${VALID}" > "${KEEP_VALID}"

########################################
# 12. Remove stale temporary outputs
########################################

rm -f "${FLAG}"

for PREFIX in "${TMP_TRAIN}" "${TMP_VALID}"; do
    rm -f \
        "${PREFIX}.raw" \
        "${PREFIX}.log" \
        "${PREFIX}.nosex"
done

########################################
# 13. Export training genotypes
########################################

echo "========================================"
echo "Exporting training genotypes"
echo "Preserving master BIM allele order"
echo "========================================"

plink \
    --bfile "${MASTER}" \
    --keep-allele-order \
    --keep "${KEEP_TRAIN}" \
    --extract "${CHIP_MARKERS}" \
    --recode A \
    --allow-extra-chr \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --memory 24000 \
    --out "${TMP_TRAIN}"

########################################
# 14. Export validation genotypes
########################################

echo "========================================"
echo "Exporting validation genotypes"
echo "Preserving master BIM allele order"
echo "========================================"

plink \
    --bfile "${MASTER}" \
    --keep-allele-order \
    --keep "${KEEP_VALID}" \
    --extract "${CHIP_MARKERS}" \
    --recode A \
    --allow-extra-chr \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --memory 24000 \
    --out "${TMP_VALID}"

########################################
# 15. Validate temporary raw files
########################################

validate_raw "${TMP_TRAIN}.raw" "${TRAIN}" train
validate_raw "${TMP_VALID}.raw" "${VALID}" valid

########################################
# 16. Exact train/validation header comparison
########################################

validate_raw_headers "${TMP_TRAIN}.raw" "${TMP_VALID}.raw"

########################################
# 17. Publish validated outputs
########################################

# Only publish after both raw files and all SNP headers pass.
# If validation fails, the previous final raw files are retained.

mv -f "${TMP_TRAIN}.raw" "${PREFIX_TRAIN}.raw"
mv -f "${TMP_VALID}.raw" "${PREFIX_VALID}.raw"

for EXT in log nosex; do
    if [[ -e "${TMP_TRAIN}.${EXT}" ]]; then
        mv -f "${TMP_TRAIN}.${EXT}" "${PREFIX_TRAIN}.${EXT}"
    fi

    if [[ -e "${TMP_VALID}.${EXT}" ]]; then
        mv -f "${TMP_VALID}.${EXT}" "${PREFIX_VALID}.${EXT}"
    fi
done

########################################
# 18. Validate published outputs again
########################################

validate_raw "${PREFIX_TRAIN}.raw" "${TRAIN}" train
validate_raw "${PREFIX_VALID}.raw" "${VALID}" valid

validate_raw_headers "${PREFIX_TRAIN}.raw" "${PREFIX_VALID}.raw"

########################################
# 19. Write completion flag
########################################

{
    echo "version=${SCRIPT_VERSION}"
    echo "outer=${OUTER}"
    echo "inner=${INNER}"
    echo "n_train=${N_TRAIN}"
    echo "n_valid=${N_VALID}"
    echo "n_chip_snps=${N_EXPECTED_SNP}"
    echo "allele_order=preserved_from_master_bim"
    echo "train_sample_ID_check=PASS"
    echo "valid_sample_ID_check=PASS"
    echo "train_SNP_ID_check=PASS"
    echo "valid_SNP_ID_check=PASS"
    echo "train_valid_snp_header_match=PASS"
    echo "train_valid_snp_order_match=PASS"
    echo "outer_test_excluded=PASS"
} > "${FLAG}.tmp"

mv -f "${FLAG}.tmp" "${FLAG}"

echo "========================================"
echo "PASS outer_${OUTER}/inner_${INNER}"
echo "Script version: ${SCRIPT_VERSION}"
echo "Training samples:   ${N_TRAIN}"
echo "Validation samples: ${N_VALID}"
echo "Chip SNPs:          ${N_EXPECTED_SNP}"
echo "Allele order:       preserved"
echo "SNP headers:        exact match"
echo "Outer-test leakage: none"
echo "========================================"

