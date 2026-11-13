#!/bin/bash
#SBATCH -J 08_imp_GWAS
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-359%6
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_run_imputed_gwas_%A_%a.out
#SBATCH -e logs/08_run_imputed_gwas_%A_%a.err

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR: line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

########################################
# 1. Modules
########################################

module load bioinfo-tools
ml gcta/1.94.1

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1

if command -v gcta64 >/dev/null 2>&1; then
    GCTA=gcta64
elif command -v gcta >/dev/null 2>&1; then
    GCTA=gcta
else
    echo "ERROR: GCTA executable not found" >&2
    exit 1
fi

########################################
# 2. Project paths
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

INPUT_ROOT="${WORK}/08_imputed_GWAS/02_imputed_inputs"

GRM_ROOT="${WORK}/08_imputed_GWAS/01_chip_GRM"

OUTROOT="${WORK}/08_imputed_GWAS/03_imputed_GWAS"

########################################
# 3. Task mapping
#
# 30 inner folds x 12 chromosomes
# = 360 array tasks.
#
# Each task runs both BB and H6.
########################################

TASK=${SLURM_ARRAY_TASK_ID}

FOLD_INDEX=$(( TASK / 12 ))
CHR=$(( TASK % 12 + 1 ))

OUTER=$(printf "%02d" $(( FOLD_INDEX / 3 + 1 )))
INNER=$(printf "%02d" $(( FOLD_INDEX % 3 + 1 )))

CHR_PAD=$(printf "%02d" "${CHR}")

echo "========================================"
echo "Step 8.2c: Imputed SNP GWAS"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "Chromosome: ${CHR}"
echo "========================================"

########################################
# 4. Input/output paths
########################################

INPUT_DIR="${INPUT_ROOT}/outer_${OUTER}/inner_${INNER}"

BFILE="${INPUT_DIR}/imputed_train_qc"

GRM_DIR="${GRM_ROOT}/outer_${OUTER}/inner_${INNER}/grm"

GRM_ALL="${GRM_DIR}/chip_all"

GRM_CHR="${GRM_DIR}/chip_chr${CHR_PAD}"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}/chr_${CHR_PAD}"

mkdir -p "${OUT}"

########################################
# 5. Check Step 8.2b completion
########################################

if [[ ! -s "${INPUT_DIR}/complete.flag" ]]; then
    echo "ERROR: Step 8.2b incomplete for this fold" >&2
    exit 1
fi

########################################
# 6. Check genotype files
########################################

for EXT in bed bim fam
do
    if [[ ! -s "${BFILE}.${EXT}" ]]; then
        echo "ERROR: Missing genotype ${EXT}" >&2
        exit 1
    fi
done

N_TRAIN=$(wc -l < "${BFILE}.fam")

if (( N_TRAIN < 500 || N_TRAIN > 600 )); then
    echo "ERROR: Unexpected sample count" >&2
    exit 1
fi

########################################
# 7. Check GRMs
########################################

for PREFIX in "${GRM_ALL}" "${GRM_CHR}"
do

    for EXT in grm.bin grm.N.bin grm.id
    do
        if [[ ! -s "${PREFIX}.${EXT}" ]]; then
            echo "ERROR: Missing GRM: ${PREFIX}.${EXT}" >&2
            exit 1
        fi
    done

    if [[ $(wc -l < "${PREFIX}.grm.id") -ne "${N_TRAIN}" ]]; then
        echo "ERROR: GRM sample count mismatch" >&2
        exit 1
    fi

done

if ! cmp -s \
    "${GRM_ALL}.grm.id" \
    "${GRM_CHR}.grm.id"
then
    echo "ERROR: GRM sample IDs/order mismatch" >&2
    exit 1
fi

########################################
# 8. Confirm expected chromosome SNPs
########################################

CHR_COUNT_FILE="${INPUT_DIR}/imputed_chr_counts.tsv"

if [[ ! -s "${CHR_COUNT_FILE}" ]]; then
    echo "ERROR: Missing chromosome SNP counts" >&2
    exit 1
fi

N_CHR_SNP=$(awk -v chr="${CHR}" '
$1 == chr {
    print $2
    found = 1
}
END {
    if (!found) exit 1
}
' "${CHR_COUNT_FILE}")

if (( N_CHR_SNP < 1000 )); then
    echo "ERROR: Too few SNPs on chromosome ${CHR}" >&2
    exit 1
fi

echo "Training individuals: ${N_TRAIN}"
echo "Expected chromosome SNPs: ${N_CHR_SNP}"

########################################
# 9. Confirm phenotype sample sets
########################################

GENO_IDS="${OUT}/genotype.samples.sorted"

awk '
{
    if (NF < 2 || $1 != $2) {
        exit 1
    }
    print $2
}
' "${BFILE}.fam" \
    | LC_ALL=C sort > "${GENO_IDS}"

for TRAIT in BB H6
do

    PHENO="${INPUT_DIR}/${TRAIT}.pheno"

    if [[ ! -s "${PHENO}" ]]; then
        echo "ERROR: Missing ${TRAIT} phenotype" >&2
        exit 1
    fi

    if [[ $(wc -l < "${PHENO}") -ne "${N_TRAIN}" ]]; then
        echo "ERROR: ${TRAIT} phenotype count mismatch" >&2
        exit 1
    fi

    PHENO_IDS="${OUT}/${TRAIT}.samples.sorted"

    awk '
    {
        if (NF != 3 || $1 != $2) {
            exit 1
        }
        print $2
    }
    ' "${PHENO}" \
        | LC_ALL=C sort > "${PHENO_IDS}"

    if ! cmp -s "${GENO_IDS}" "${PHENO_IDS}"; then
        echo "ERROR: ${TRAIT} phenotype ID mismatch" >&2
        exit 1
    fi

done

echo "PASS: Genotype, GRM and phenotype inputs"

########################################
# 10. Reuse completed analyses
########################################

BB_RESULT="${OUT}/BB_chr${CHR_PAD}.mlma"

H6_RESULT="${OUT}/H6_chr${CHR_PAD}.mlma"

if [[ -s "${OUT}/complete.flag" &&
      -s "${BB_RESULT}" &&
      -s "${H6_RESULT}" &&
      -s "${OUT}/gwas_qc.tsv" ]]; then

    echo "All outputs already completed; skipping."
    echo "DONE Step8.2c (previously completed)"
    exit 0
fi

rm -f "${OUT}/complete.flag"

########################################
# 11. Run LOCO GWAS
########################################

QC_FILE="${OUT}/gwas_qc.tsv"

printf "Outer\tInner\tTrait\tChromosome\tExpected_SNPs\tTested_SNPs\tValid_P\tInvalid_P\n" \
    > "${QC_FILE}"

for TRAIT in BB H6
do

    PHENO="${INPUT_DIR}/${TRAIT}.pheno"

    PREFIX="${OUT}/${TRAIT}_chr${CHR_PAD}"

    echo "========================================"
    echo "Running ${TRAIT}, chromosome ${CHR}"
    echo "========================================"

    echo "GCTA: ${GCTA}"
    echo "BFILE: ${BFILE}"
    echo "Full GRM: ${GRM_ALL}"
    echo "Subtracted GRM: ${GRM_CHR}"
    echo "Phenotype: ${PHENO}"

    ####################################
    # GWAS:
    #
    # --bfile:
    #   Imputed markers only.
    #
    # --grm:
    #   Whole-chip GRM.
    #
    # --mlma-subtract-grm:
    #   Remove chip GRM component
    #   from the tested chromosome.
    #
    # --chr:
    #   Test only that chromosome.
    ####################################

    "${GCTA}" \
        --mlma \
        --bfile "${BFILE}" \
        --autosome-num 12 \
        --chr "${CHR}" \
        --grm "${GRM_ALL}" \
        --mlma-subtract-grm "${GRM_CHR}" \
        --pheno "${PHENO}" \
        --thread-num 4 \
        --out "${PREFIX}"

    ####################################
    # Verify GCTA output
    ####################################

    RESULT="${PREFIX}.mlma"

    if [[ ! -s "${RESULT}" ]]; then
        echo "ERROR: Missing GWAS result ${RESULT}" >&2
        exit 1
    fi

    # The standard GCTA MLMA output uses:
    #
    # Chr SNP bp A1 A2 Freq b se p
    #
    # Some problematic markers can have
    # missing/non-finite P-values.
    # Record those rather than silently
    # treating them as valid associations.

    awk \
        -v outer="${OUTER}" \
        -v inner="${INNER}" \
        -v trait="${TRAIT}" \
        -v chr="${CHR}" \
        -v expected="${N_CHR_SNP}" '
    function is_probability(x) {
        if (x !~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ &&
            x !~ /^[.][0-9]+([eE][-+]?[0-9]+)?$/) {
            return 0
        }

        return (x + 0 >= 0 && x + 0 <= 1)
    }

    NR == 1 {
        if (NF < 9 || $1 != "Chr" || $2 != "SNP" || $9 != "p") {
            print "ERROR: Unexpected GWAS header" > "/dev/stderr"
            bad = 1
            exit 1
        }
        next
    }

    {
        if (NF < 9 || $1 != chr) {
            print "ERROR: Incorrect GWAS chromosome/columns" > "/dev/stderr"
            bad = 1
            exit 1
        }

        tested++

        if (is_probability($9)) {
            valid++
        } else {
            invalid++
        }
    }

    END {
        if (bad) exit 1

        if (tested < 1000 || valid < 1000 || tested > expected) {
            print "ERROR: Unexpected number of GWAS results" > "/dev/stderr"
            exit 1
        }

        printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n", \
            outer, inner, trait, chr, expected, tested, valid, invalid
    }
    ' "${RESULT}" >> "${QC_FILE}"

    echo "PASS: ${TRAIT} chr${CHR} GWAS output"

done

########################################
# 12. Save completion status
########################################

printf "DONE Step8.2c outer_%s inner_%s chr_%s\n" \
    "${OUTER}" "${INNER}" "${CHR_PAD}" \
    > "${OUT}/complete.flag"

echo "========================================"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "Chromosome: ${CHR}"
echo "Results: ${BB_RESULT}"
echo "Results: ${H6_RESULT}"
echo "DONE Step8.2c"
echo "========================================"

