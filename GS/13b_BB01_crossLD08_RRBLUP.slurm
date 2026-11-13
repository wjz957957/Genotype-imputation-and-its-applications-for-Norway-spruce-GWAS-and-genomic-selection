#!/bin/bash
#SBATCH -J 13b_BB01_LD08
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-2%3
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/13b_BB01_LD08_%A_%a.out
#SBATCH -e logs/13b_BB01_LD08_%A_%a.err

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

########################################
# Software
########################################

module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

########################################
# Paths and fold mapping
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

SCRIPT="${WORK}/scripts/13b_BB01_crossLD08_RRBLUP.R"

TASK=${SLURM_ARRAY_TASK_ID}

OUTER=01
INNER=$(printf "%02d" $((TASK + 1)))
[[ "$TASK" =~ ^[012]$ ]] || { echo "ERROR: pilot array index must be 0,1,2" >&2; exit 1; }

# Step 13a and Step 13b use the same isolated pilot fold directory.
INPUT_OUT="${WORK}/09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_pilot/outer_${OUTER}/inner_${INNER}"
OUT="${INPUT_OUT}"

FLAG="${OUT}/step9b.complete.flag"

echo "========================================"
echo "Step 13b: BB cross-chip LD filtered RRBLUP"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

########################################
# Input checks
########################################

[[ -s "${SCRIPT}" ]] || {
    echo "ERROR: Missing R script" >&2
    exit 1
}

[[ -s "${INPUT_OUT}/step9a.complete.flag" ]] || {
    echo "ERROR: Step 9a incomplete" >&2
    exit 1
}

# Skip only if the BB result contains header + 20 candidate settings.
RESULTS_READY=1

for TRAIT in BB
do
    FILE="${OUT}/${TRAIT}_chip_plus_imputed_accuracy.tsv"

    if [[ ! -s "${FILE}" ]]; then
        RESULTS_READY=0
        continue
    fi

    N_LINES=$(wc -l < "${FILE}")
    if (( N_LINES != 21 )); then
        RESULTS_READY=0
    fi
done

if [[ -s "${FLAG}" ]] && (( RESULTS_READY == 1 )); then
    echo "SKIP: Cross-CHIP LD pilot Step13b already completed"
    exit 0
fi

mkdir -p "${OUT}"
rm -f "${FLAG}"

########################################
# Run RRBLUP
########################################

cd "${WORK}"

Rscript "${SCRIPT}" \
    "${WORK}" \
    "${OUTER}" \
    "${INNER}"

########################################
# Result validation
########################################

for TRAIT in BB
do
    FILE="${OUT}/${TRAIT}_chip_plus_imputed_accuracy.tsv"

    [[ -s "${FILE}" ]] || {
        echo "ERROR: Missing result ${FILE}" >&2
        exit 1
    }

    awk '
    END {
        if (NR != 21) {
            print "ERROR: Expected header + 20 RRBLUP results" > "/dev/stderr"
            exit 1
        }
    }
    ' "${FILE}"

    echo "PASS: ${TRAIT} has 20 RRBLUP results"
done

########################################
# Completion
########################################

printf "DONE Step13b outer_%s inner_%s\n" \
    "${OUTER}" "${INNER}" > "${FLAG}"

echo "========================================"
echo "DONE Step13b"
echo "========================================"
