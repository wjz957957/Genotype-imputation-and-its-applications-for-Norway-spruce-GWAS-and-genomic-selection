#!/bin/bash
#SBATCH -J 10d_outer_final
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-19%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/10d_outer_final_%A_%a.out
#SBATCH -e logs/10d_outer_final_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

########################################
# 1. Modules
########################################

module load bioinfo-tools
ml plink/1.90b4.9

########################################
# 2. Paths
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"

CV_ROOT="${WORK}/02_CV"
CHIP_ROOT="${WORK}/10_outer_test/01_chip_selection"
OPT_ROOT="${WORK}/09_chip_plus_imputed_RRBLUP/optimization"
QTN_ROOT="${WORK}/10_outer_test/04_imputed_clumping"

OUT_ROOT="${WORK}/10_outer_test/05_final_genotypes"

TASK=${SLURM_ARRAY_TASK_ID:?ERROR: submit as Slurm array}

if (( TASK < 0 || TASK > 19 )); then
    echo "ERROR: task must be 0..19" >&2
    exit 1
fi

########################################
# 3. Map array task
########################################

OUTER=$(printf '%02d' "$((TASK / 2 + 1))")

if (( TASK % 2 == 0 )); then
    TRAIT=BB
else
    TRAIT=H6
fi

CV="${CV_ROOT}/outer_${OUTER}"

CHIP_DIR="${CHIP_ROOT}/outer_${OUTER}/selected"
QTN_DIR="${QTN_ROOT}/outer_${OUTER}/${TRAIT}"
OPT_DIR="${OPT_ROOT}/outer_${OUTER}"

TRAIN="${CV}/outer_train.samples"
TEST="${CV}/outer_test.samples"

CHIP_LIST="${CHIP_DIR}/${TRAIT}_best_chip.snps"
CHIP_SETTING="${CHIP_DIR}/${TRAIT}_best_marker.txt"

BEST_EXTRA_FILE="${OPT_DIR}/${TRAIT}_best_extra.txt"

OUT="${OUT_ROOT}/outer_${OUTER}/${TRAIT}"

FLAG="${OUT}/complete.flag"

mkdir -p "${OUT}"

rm -f "${FLAG}"

echo "========================================"
echo "Step 10d: final outer-train/test genotype"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Trait: ${TRAIT}"
echo "========================================"

########################################
# 4. Check required inputs
########################################

for F in \
    "${MASTER}.bed" \
    "${MASTER}.bim" \
    "${MASTER}.fam" \
    "${TRAIN}" \
    "${TEST}" \
    "${CHIP_LIST}" \
    "${CHIP_SETTING}" \
    "${BEST_EXTRA_FILE}"
do
    [[ -s "${F}" ]] || {
        echo "ERROR: missing/empty input: ${F}" >&2
        exit 1
    }
done

########################################
# 5. Read selected SNP numbers
########################################

BEST_CHIP=$(tr -d '[:space:]' < "${CHIP_SETTING}")
BEST_EXTRA=$(tr -d '[:space:]' < "${BEST_EXTRA_FILE}")

case "${BEST_EXTRA}" in
    0|1|5|10|20|50|100|200|500|1000)
        ;;
    *)
        echo "ERROR: invalid Best_Extra=${BEST_EXTRA}" >&2
        exit 1
        ;;
esac

########################################
# 6. Locate selected extra SNP set
########################################

EXTRA_LIST="${QTN_DIR}/marker_sets/extra_${BEST_EXTRA}.snps"

[[ -f "${EXTRA_LIST}" ]] || {
    echo "ERROR: missing extra SNP list: ${EXTRA_LIST}" >&2
    exit 1
}

########################################
# 7. Count selected markers
########################################

N_CHIP=$(wc -l < "${CHIP_LIST}")
N_EXTRA=$(wc -l < "${EXTRA_LIST}")

if (( N_CHIP < 1 )); then
    echo "ERROR: empty selected chip set" >&2
    exit 1
fi

########################################
# 8. Special handling of zero-QTN case
########################################

if (( BEST_EXTRA == 0 )); then

    if (( N_EXTRA != 0 )); then
        echo "ERROR: Best_Extra=0 but extra_0.snps is not empty" >&2
        exit 1
    fi

    echo "INFO: Best_Extra=0"
    echo "INFO: no additional QTN will be included"

else

    if (( N_EXTRA != BEST_EXTRA )); then
        echo "ERROR: extra count mismatch: ${N_EXTRA} vs ${BEST_EXTRA}" >&2
        exit 1
    fi

fi

########################################
# 9. Check marker uniqueness
########################################

awk '
{
    sub(/\r$/, "", $0)

    if (NF != 1 || seen[$1]++) {
        print "ERROR: invalid/duplicate SNP: " $0 > "/dev/stderr"
        bad=1
    }
}
END {
    if (bad)
        exit 1
}
' "${CHIP_LIST}"

if (( N_EXTRA > 0 )); then

    awk '
    {
        sub(/\r$/, "", $0)

        if (NF != 1 || seen[$1]++) {
            print "ERROR: invalid/duplicate extra SNP: " $0 > "/dev/stderr"
            bad=1
        }
    }
    END {
        if (bad)
            exit 1
    }
    ' "${EXTRA_LIST}"

fi

########################################
# 10. Check chip / extra separation
########################################

if (( N_EXTRA > 0 )); then

    OVERLAP=$(comm -12 \
        <(sort "${CHIP_LIST}") \
        <(sort "${EXTRA_LIST}") \
        | wc -l)

    if (( OVERLAP != 0 )); then
        echo "ERROR: chip/extra SNP overlap = ${OVERLAP}" >&2
        exit 1
    fi

    echo "PASS: chip/extra overlap = 0"

fi

########################################
# 11. Check selected markers exist
#     in master BIM
########################################

check_markers_in_bim() {

    local LIST="$1"

    awk -v list_file="${LIST}" '
    FILENAME == list_file {
        wanted[$1]=1
        total++
        next
    }

    ($2 in wanted) {
        found[$2]++
        next
    }

    END {

        missing=0

        for (id in wanted) {
            if (!(id in found)) {
                print id
                missing++
            }
        }

        if (missing > 0 || length(wanted) != total)
            exit 1
    }
    ' "${LIST}" "${MASTER}.bim" > "${OUT}/.missing.tmp" || {

        echo "ERROR: selected marker absent from master BIM: ${LIST}" >&2

        cat "${OUT}/.missing.tmp" >&2

        rm -f "${OUT}/.missing.tmp"

        exit 1
    }

    rm -f "${OUT}/.missing.tmp"

    echo "PASS: all markers from ${LIST} exist in master BIM"

}

check_markers_in_bim "${CHIP_LIST}"

if (( N_EXTRA > 0 )); then
    check_markers_in_bim "${EXTRA_LIST}"
else
    echo "INFO: skipping master-BIM check for zero extra QTN"
fi

########################################
# 12. Prepare train / test keep files
########################################

awk 'NF==1 {print $1, $1}' \
    "${TRAIN}" > "${OUT}/train.keep"

awk 'NF==1 {print $1, $1}' \
    "${TEST}" > "${OUT}/test.keep"

N_TRAIN=$(wc -l < "${OUT}/train.keep")
N_TEST=$(wc -l < "${OUT}/test.keep")

if (( N_TRAIN + N_TEST != 904 )); then
    echo "ERROR: outer train + test != 904" >&2
    exit 1
fi

########################################
# 13. Extract genotype data
#
# No second MAF/missingness filtering.
########################################

extract_raw() {

    local list="$1"
    local keep="$2"
    local prefix="$3"
    local expected="$4"

    if (( expected == 0 )); then
        return 0
    fi

    rm -f "${prefix}.raw"

    plink \
        --bfile "${MASTER}" \
        --extract "${list}" \
        --keep "${keep}" \
        --keep-allele-order \
        --allow-extra-chr \
        --recode A \
        --out "${prefix}" \
        --threads 4 \
        --memory 24000

    [[ -s "${prefix}.raw" ]] || {
        echo "ERROR: missing ${prefix}.raw" >&2
        exit 1
    }

    ACTUAL=$(awk 'NR==1 {print NF-6; exit}' "${prefix}.raw")

    if (( ACTUAL != expected )); then
        echo "ERROR: ${prefix}: expected ${expected} SNPs, got ${ACTUAL}" >&2
        exit 1
    fi
}

########################################
# 14. Extract chip SNPs
########################################

extract_raw \
    "${CHIP_LIST}" \
    "${OUT}/train.keep" \
    "${OUT}/chip_train" \
    "${N_CHIP}"

extract_raw \
    "${CHIP_LIST}" \
    "${OUT}/test.keep" \
    "${OUT}/chip_test" \
    "${N_CHIP}"

########################################
# 15. Extract extra/QTN SNPs
########################################

if (( N_EXTRA > 0 )); then

    extract_raw \
        "${EXTRA_LIST}" \
        "${OUT}/train.keep" \
        "${OUT}/extra_train" \
        "${N_EXTRA}"

    extract_raw \
        "${EXTRA_LIST}" \
        "${OUT}/test.keep" \
        "${OUT}/extra_test" \
        "${N_EXTRA}"

else

    echo "INFO: no extra_train.raw / extra_test.raw because Best_Extra=0"

fi

########################################
# 16. Validate row counts
########################################

for F in \
    "${OUT}/chip_train.raw" \
    "${OUT}/chip_test.raw"
do

    rows=$(( $(wc -l < "${F}") - 1 ))

    if [[ "${F}" == *chip_train.raw ]]; then
        expected_rows=${N_TRAIN}
    else
        expected_rows=${N_TEST}
    fi

    (( rows == expected_rows )) || {
        echo "ERROR: row count mismatch in ${F}: ${rows} vs ${expected_rows}" >&2
        exit 1
    }

done

if (( N_EXTRA > 0 )); then

    for F in \
        "${OUT}/extra_train.raw" \
        "${OUT}/extra_test.raw"
    do

        rows=$(( $(wc -l < "${F}") - 1 ))

        if [[ "${F}" == *extra_train.raw ]]; then
            expected_rows=${N_TRAIN}
        else
            expected_rows=${N_TEST}
        fi

        (( rows == expected_rows )) || {
            echo "ERROR: row count mismatch in ${F}: ${rows} vs ${expected_rows}" >&2
            exit 1
        }

    done

fi

########################################
# 17. Save manifest
########################################

{
    printf 'Outer\tTrait\tBest_Chip\tN_chip\tBest_Extra\tN_total\tN_train\tN_test\n'

    printf '%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n' \
        "${OUTER}" \
        "${TRAIT}" \
        "${BEST_CHIP}" \
        "${N_CHIP}" \
        "${BEST_EXTRA}" \
        "$((N_CHIP + N_EXTRA))" \
        "${N_TRAIN}" \
        "${N_TEST}"

} > "${OUT}/final_genotype_manifest.tsv"

########################################
# 18. Completion flag
########################################

printf \
    'DONE Step10d outer_%s trait_%s\n' \
    "${OUTER}" \
    "${TRAIT}" \
    > "${FLAG}"

echo "PASS: final genotype inputs prepared"

cat "${OUT}/final_genotype_manifest.tsv"
