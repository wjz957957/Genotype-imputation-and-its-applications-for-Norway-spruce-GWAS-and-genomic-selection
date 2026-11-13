#!/bin/bash
#SBATCH -J 10b1_outer_inputs
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-9%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/10b1_outer_inputs_%A_%a.out
#SBATCH -e logs/10b1_outer_inputs_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

########################################
# 1. Modules
########################################

module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1

if command -v gcta64 >/dev/null 2>&1; then
    GCTA=gcta64
elif command -v gcta >/dev/null 2>&1; then
    GCTA=gcta
else
    echo "ERROR: GCTA executable unavailable" >&2
    exit 1
fi

########################################
# 2. Paths
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"

CHIP_LIST="${WORK}/01_reference/chip.marker"
IMP_LIST="${WORK}/01_reference/imputed.marker"

TASK=${SLURM_ARRAY_TASK_ID}
OUTER=$(printf "%02d" $((TASK + 1)))

CV="${WORK}/02_CV/outer_${OUTER}"

TRAIN="${CV}/outer_train.samples"
TEST="${CV}/outer_test.samples"

CHIP_STEP10A="${WORK}/10_outer_test/01_chip_selection/outer_${OUTER}"

OUT="${WORK}/10_outer_test/02_imputed_inputs/outer_${OUTER}"

PLINK_DIR="${OUT}/plink"
GRM_DIR="${OUT}/grm"
PHENO_DIR="${OUT}/pheno"

CHIP="${PLINK_DIR}/chip_train_qc"
IMP="${PLINK_DIR}/imputed_train_qc"

GRM_ALL="${GRM_DIR}/chip_all"

KEEP="${OUT}/outer_train.keep"
FLAG="${OUT}/complete.flag"

mkdir -p "${PLINK_DIR}" "${GRM_DIR}" "${PHENO_DIR}"

echo "========================================"
echo "Step 10b-1: Outer-train imputed inputs"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "========================================"

########################################
# 3. Check completion
########################################

if [[ -s "${FLAG}" ]]; then
    READY=1

    for PREFIX in "${CHIP}" "${IMP}"
    do
        for EXT in bed bim fam
        do
            [[ -s "${PREFIX}.${EXT}" ]] || READY=0
        done
    done

    for PREFIX in "${GRM_ALL}" \
        "${GRM_DIR}"/chip_chr{01..12}
    do
        for EXT in grm.bin grm.N.bin grm.id
        do
            [[ -s "${PREFIX}.${EXT}" ]] || READY=0
        done
    done

    for TRAIT in BB H6
    do
        [[ -s "${PHENO_DIR}/${TRAIT}.pheno" ]] || READY=0
    done

    [[ -s "${OUT}/input_summary.tsv" ]] || READY=0

    if (( READY == 1 )); then
        echo "SKIP: Previously completed"
        exit 0
    fi
fi

rm -f "${FLAG}"

########################################
# 4. Validate inputs
########################################

for EXT in bed bim fam
do
    [[ -s "${MASTER}.${EXT}" ]] || {
        echo "ERROR: Missing master ${EXT}" >&2
        exit 1
    }
done

for FILE in \
    "${CHIP_LIST}" \
    "${IMP_LIST}" \
    "${TRAIN}" \
    "${TEST}" \
    "${CHIP_STEP10A}/complete.flag" \
    "${CHIP_STEP10A}/pheno/BB.pheno" \
    "${CHIP_STEP10A}/pheno/H6.pheno"
do
    [[ -s "${FILE}" ]] || {
        echo "ERROR: Missing prerequisite ${FILE}" >&2
        exit 1
    }
done

if [[ $(wc -l < "${CHIP_LIST}") -ne 43267 ]]; then
    echo "ERROR: Unexpected chip marker count" >&2
    exit 1
fi

########################################
# 5. Prepare outer-train keep file
########################################

awk '
{
    sub(/\r$/, "", $0)

    if (NF != 1 || seen[$1]++) {
        print "ERROR: Invalid training ID" > "/dev/stderr"
        bad=1
        exit 1
    }

    print $1, $1
}
END {
    if (bad || NR == 0) exit 1
}
' "${TRAIN}" > "${KEEP}"

N_TRAIN=$(wc -l < "${KEEP}")
N_TEST=$(wc -l < "${TEST}")

if (( N_TRAIN < 813 || N_TRAIN > 814 )); then
    echo "ERROR: Unexpected outer-train sample count" >&2
    exit 1
fi

if (( N_TRAIN + N_TEST != 904 )); then
    echo "ERROR: Outer sample total is not 904" >&2
    exit 1
fi

awk '
NR == FNR {
    sub(/\r$/, "", $0)
    train[$1]=1
    next
}
{
    sub(/\r$/, "", $0)

    if (NF != 1 || test_seen[$1]++) {
        print "ERROR: Invalid outer-test ID" > "/dev/stderr"
        bad=1
    }

    if ($1 in train) {
        print "ERROR: Outer train/test overlap" > "/dev/stderr"
        bad=1
    }
}
END {
    if (bad) exit 1
}
' "${TRAIN}" "${TEST}"

awk '{print $1}' "${KEEP}" \
    | sort > "${OUT}/expected.samples.sorted"

echo "Outer train: ${N_TRAIN}"
echo "Outer test: ${N_TEST}"

########################################
# 6. Extract QC chip SNPs
#
# EXACT same QC as Step 8.2a:
# --maf 0.01
# --geno 0.05
########################################

echo "========================================"
echo "Preparing QC chip genotype"
echo "========================================"

plink \
    --bfile "${MASTER}" \
    --keep "${KEEP}" \
    --extract "${CHIP_LIST}" \
    --maf 0.01 \
    --geno 0.05 \
    --keep-allele-order \
    --allow-extra-chr \
    --make-bed \
    --threads 4 \
    --memory 24000 \
    --out "${CHIP}"

########################################
# 7. Extract QC imputed SNPs
#
# EXACT same QC as Step 8.2b.
########################################

echo "========================================"
echo "Preparing QC imputed genotype"
echo "========================================"

plink \
    --bfile "${MASTER}" \
    --keep "${KEEP}" \
    --extract "${IMP_LIST}" \
    --maf 0.01 \
    --geno 0.05 \
    --keep-allele-order \
    --allow-extra-chr \
    --make-bed \
    --threads 4 \
    --memory 24000 \
    --out "${IMP}"

########################################
# 8. Validate genotype dimensions
########################################

check_bfile () {
    local PREFIX="$1"
    local MIN_SNPS="$2"

    for EXT in bed bim fam
    do
        [[ -s "${PREFIX}.${EXT}" ]] || {
            echo "ERROR: Missing ${PREFIX}.${EXT}" >&2
            exit 1
        }
    done

    local N_SAMPLES
    local N_SNPS
    local BYTES_PER_SNP
    local EXPECTED_BYTES
    local ACTUAL_BYTES

    N_SAMPLES=$(wc -l < "${PREFIX}.fam")
    N_SNPS=$(wc -l < "${PREFIX}.bim")

    if (( N_SAMPLES != N_TRAIN || N_SNPS < MIN_SNPS )); then
        echo "ERROR: Invalid genotype dimensions for ${PREFIX}" >&2
        exit 1
    fi

    BYTES_PER_SNP=$(( (N_SAMPLES + 3) / 4 ))
    EXPECTED_BYTES=$(( N_SNPS * BYTES_PER_SNP + 3 ))
    ACTUAL_BYTES=$(stat -c '%s' "${PREFIX}.bed")

    if (( ACTUAL_BYTES != EXPECTED_BYTES )); then
        echo "ERROR: BED size mismatch for ${PREFIX}" >&2
        exit 1
    fi

    awk '
    {
        if (NF < 2 || $1 != $2 || seen[$2]++) {
            print "ERROR: Invalid genotype sample ID" > "/dev/stderr"
            bad=1
            exit 1
        }

        print $2
    }
    END {
        if (bad) exit 1
    }
    ' "${PREFIX}.fam" \
        | sort > "${PREFIX}.samples.sorted"

    if ! cmp -s \
        "${PREFIX}.samples.sorted" \
        "${OUT}/expected.samples.sorted"
    then
        echo "ERROR: Genotype sample identities mismatch" >&2
        exit 1
    fi

    echo "PASS genotype: ${PREFIX}"
    echo "Samples: ${N_SAMPLES}"
    echo "Markers: ${N_SNPS}"
}

check_bfile "${CHIP}" 10000
check_bfile "${IMP}" 100000

N_CHIP=$(wc -l < "${CHIP}.bim")
N_IMP=$(wc -l < "${IMP}.bim")

########################################
# 9. Validate marker membership
########################################

for KIND in CHIP IMP
do
    if [[ "${KIND}" == "CHIP" ]]; then
        LIST="${CHIP_LIST}"
        BIM="${CHIP}.bim"
    else
        LIST="${IMP_LIST}"
        BIM="${IMP}.bim"
    fi

    awk '
    NR == FNR {
        if (NF != 1 || allowed[$1]++) {
            print "ERROR: Invalid source SNP list" > "/dev/stderr"
            bad=1
        }
        next
    }

    {
        if (!($2 in allowed) || observed[$2]++) {
            print "ERROR: Incorrect or duplicate SNP in BIM" > "/dev/stderr"
            bad=1
        }
    }

    END {
        if (bad) exit 1
    }
    ' "${LIST}" "${BIM}"

    echo "PASS: ${KIND} SNP membership"
done

########################################
# 10. Chromosome distribution
########################################

for KIND in CHIP IMP
do
    if [[ "${KIND}" == "CHIP" ]]; then
        BIM="${CHIP}.bim"
    else
        BIM="${IMP}.bim"
    fi

    awk '
    {
        if ($1 !~ /^[0-9]+$/ || $1+0 < 1 || $1+0 > 12) {
            print "ERROR: Unexpected chromosome:", $1 > "/dev/stderr"
            bad=1
        }

        count[$1]++
        total++
    }

    END {
        for (chr=1; chr<=12; chr++) {
            print chr "\t" count[chr]+0

            if (count[chr]+0 < 1000) {
                print "ERROR: Too few markers on chromosome", chr > "/dev/stderr"
                bad=1
            }
        }

        if (bad || total < 10000) exit 1
    }
    ' "${BIM}" \
        > "${OUT}/${KIND,,}_chr_counts.tsv"

    echo "PASS: ${KIND} chromosome distribution"
done

########################################
# 11. Reuse Step 10a training phenotypes
#
# Step 10a already normalized CRLF and
# checked phenotype values.
########################################

echo "========================================"
echo "Preparing training-only phenotypes"
echo "========================================"

for TRAIT in BB H6
do
    SOURCE="${CHIP_STEP10A}/pheno/${TRAIT}.pheno"
    TARGET="${PHENO_DIR}/${TRAIT}.pheno"

    cp "${SOURCE}" "${TARGET}"

    awk -v expected="${N_TRAIN}" '
    {
        if (NF != 3 || $1 != $2 || seen[$2]++) {
            print "ERROR: Invalid phenotype row" > "/dev/stderr"
            bad=1
            exit 1
        }

        if ($3 !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
            print "ERROR: Invalid numeric phenotype" > "/dev/stderr"
            bad=1
            exit 1
        }

        print $2
        n++
    }

    END {
        if (bad || n != expected) exit 1
    }
    ' "${TARGET}" \
        | sort > "${PHENO_DIR}/${TRAIT}.samples.sorted"

    if ! cmp -s \
        "${PHENO_DIR}/${TRAIT}.samples.sorted" \
        "${OUT}/expected.samples.sorted"
    then
        echo "ERROR: ${TRAIT} phenotype sample mismatch" >&2
        exit 1
    fi

    echo "PASS: ${TRAIT} phenotype ${N_TRAIN} samples"
done

########################################
# 12. Build full-chip QC GRM
########################################

echo "========================================"
echo "Building QC full-chip GRM"
echo "========================================"

"${GCTA}" \
    --bfile "${CHIP}" \
    --autosome-num 12 \
    --make-grm \
    --thread-num 4 \
    --out "${GRM_ALL}"

for EXT in grm.bin grm.N.bin grm.id
do
    [[ -s "${GRM_ALL}.${EXT}" ]] || {
        echo "ERROR: Missing full-chip GRM ${EXT}" >&2
        exit 1
    }
done

if [[ $(wc -l < "${GRM_ALL}.grm.id") -ne "${N_TRAIN}" ]]; then
    echo "ERROR: Full-chip GRM sample count mismatch" >&2
    exit 1
fi

########################################
# 13. Build 12 chromosome-specific GRMs
########################################

for CHR in $(seq 1 12)
do
    CHR_PAD=$(printf "%02d" "${CHR}")
    GRM_CHR="${GRM_DIR}/chip_chr${CHR_PAD}"

    echo "========================================"
    echo "Building chip GRM chromosome ${CHR}"
    echo "========================================"

    "${GCTA}" \
        --bfile "${CHIP}" \
        --autosome-num 12 \
        --chr "${CHR}" \
        --make-grm \
        --thread-num 4 \
        --out "${GRM_CHR}"

    for EXT in grm.bin grm.N.bin grm.id
    do
        [[ -s "${GRM_CHR}.${EXT}" ]] || {
            echo "ERROR: Missing chr${CHR} GRM ${EXT}" >&2
            exit 1
        }
    done

    if [[ $(wc -l < "${GRM_CHR}.grm.id") -ne "${N_TRAIN}" ]]; then
        echo "ERROR: Chromosome GRM sample count mismatch" >&2
        exit 1
    fi

    if ! cmp -s \
        "${GRM_ALL}.grm.id" \
        "${GRM_CHR}.grm.id"
    then
        echo "ERROR: Chromosome GRM IDs/order mismatch" >&2
        exit 1
    fi

    echo "PASS: Chromosome ${CHR} GRM"
done

########################################
# 14. GRM/genotype identity check
########################################

awk '
{
    if (NF != 2 || $1 != $2 || seen[$2]++) {
        print "ERROR: Invalid GRM sample ID" > "/dev/stderr"
        bad=1
        exit 1
    }

    print $2
}
END {
    if (bad) exit 1
}
' "${GRM_ALL}.grm.id" \
    | sort > "${OUT}/grm.samples.sorted"

if ! cmp -s \
    "${OUT}/grm.samples.sorted" \
    "${OUT}/expected.samples.sorted"
then
    echo "ERROR: GRM and training sample mismatch" >&2
    exit 1
fi

if ! cmp -s \
    "${CHIP}.samples.sorted" \
    "${IMP}.samples.sorted"
then
    echo "ERROR: Chip/imputed sample mismatch" >&2
    exit 1
fi

########################################
# 15. Save summary
########################################

SUMMARY="${OUT}/input_summary.tsv"

{
    printf 'Outer\tN_outer_train\tN_outer_test\tN_chip_QC\tN_imputed_QC\tN_GRMs\n'

    printf '%s\t%s\t%s\t%s\t%s\t13\n' \
        "${OUTER}" \
        "${N_TRAIN}" \
        "${N_TEST}" \
        "${N_CHIP}" \
        "${N_IMP}"
} > "${SUMMARY}"

########################################
# 16. Completion flag
########################################

printf 'DONE Step10b-1 outer_%s\n' "${OUTER}" \
    > "${FLAG}"

echo "========================================"
echo "DONE Step10b-1"
echo "Outer: ${OUTER}"
echo "Train samples: ${N_TRAIN}"
echo "QC chip SNPs: ${N_CHIP}"
echo "QC imputed SNPs: ${N_IMP}"
echo "GRMs built: 13"
echo "========================================"

