#!/bin/bash
#SBATCH -J 08_imp_input
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-29%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_prepare_imputed_input_%A_%a.out
#SBATCH -e logs/08_prepare_imputed_input_%A_%a.err

set -Eeuo pipefail

trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

############################################
# 1. Modules
############################################

module load bioinfo-tools
ml plink/1.90b4.9

############################################
# 2. Paths
############################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"

IMP_LIST="${WORK}/01_reference/imputed.marker"

CV_DIR="${WORK}/02_CV"

GRM_ROOT="${WORK}/08_imputed_GWAS/01_chip_GRM"

OUTROOT="${WORK}/08_imputed_GWAS/02_imputed_inputs"

############################################
# 3. Task mapping
############################################

TASK=${SLURM_ARRAY_TASK_ID}

OUTER=$(printf "%02d" $((TASK / 3 + 1)))
INNER=$(printf "%02d" $((TASK % 3 + 1)))

echo "========================================"
echo "Step 8.2b"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

############################################
# 4. Fold-specific paths
############################################

TRAIN_SAMPLES="${CV_DIR}/outer_${OUTER}/inner_${INNER}/inner_train.samples"

GRM_DIR="${GRM_ROOT}/outer_${OUTER}/inner_${INNER}"

GRM_PREFIX="${GRM_DIR}/grm/chip_all"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}"

BFILE="${OUT}/imputed_train_qc"

KEEP="${OUT}/train.keep"

mkdir -p "${OUT}"

rm -f "${OUT}/complete.flag"

############################################
# 5. Input checks
############################################

echo "Checking inputs"

for EXT in bed bim fam
do
    test -s "${MASTER}.${EXT}" || {
        echo "ERROR: Missing master ${EXT}" >&2
        exit 1
    }
done

test -s "${IMP_LIST}" || {
    echo "ERROR: Missing imputed.marker" >&2
    exit 1
}

test -s "${TRAIN_SAMPLES}" || {
    echo "ERROR: Missing inner_train.samples" >&2
    exit 1
}

# Previous GRM completion flags are zero bytes.
# Therefore use -f, not -s.

test -f "${GRM_DIR}/complete.flag" || {
    echo "ERROR: GRM step incomplete" >&2
    exit 1
}

for EXT in grm.bin grm.N.bin grm.id
do
    test -s "${GRM_PREFIX}.${EXT}" || {
        echo "ERROR: Missing chip GRM ${EXT}" >&2
        exit 1
    }
done

echo "PASS: Input and GRM files"

############################################
# 6. Prepare keep file
############################################

N_EXPECTED=$(wc -l < "${TRAIN_SAMPLES}")

echo "Expected training samples: ${N_EXPECTED}"

if (( N_EXPECTED < 500 || N_EXPECTED > 600 )); then
    echo "ERROR: Unexpected training sample count" >&2
    exit 1
fi

awk '
{
    if (NF != 1 || seen[$1]++) {
        print "ERROR: Invalid sample ID file" > "/dev/stderr"
        exit 1
    }

    print $1, $1
}
' "${TRAIN_SAMPLES}" > "${KEEP}"

if [[ $(wc -l < "${KEEP}") -ne "${N_EXPECTED}" ]]; then
    echo "ERROR: Incorrect keep count" >&2
    exit 1
fi

############################################
# 7. Reuse existing genotype if valid
############################################

REUSE=0

if [[ -s "${BFILE}.bed" &&
      -s "${BFILE}.bim" &&
      -s "${BFILE}.fam" ]]; then

    OLD_N_SAMPLE=$(wc -l < "${BFILE}.fam")
    OLD_N_SNP=$(wc -l < "${BFILE}.bim")

    BYTES_PER_SNP=$(( (OLD_N_SAMPLE + 3) / 4 ))
    EXPECTED_BYTES=$(( OLD_N_SNP * BYTES_PER_SNP + 3 ))

    ACTUAL_BYTES=$(stat -c '%s' "${BFILE}.bed")

    if [[ "${OLD_N_SAMPLE}" -eq "${N_EXPECTED}" &&
          "${OLD_N_SNP}" -ge 100000 &&
          "${ACTUAL_BYTES}" -eq "${EXPECTED_BYTES}" ]]; then

        awk '
        {
            if ($1 != $2) exit 1
            print $2
        }
        ' "${BFILE}.fam" > "${OUT}/reuse.samples"

        LC_ALL=C sort "${OUT}/reuse.samples" \
            > "${OUT}/reuse.samples.sorted"

        LC_ALL=C sort "${TRAIN_SAMPLES}" \
            > "${OUT}/expected.samples.sorted"

        if cmp -s \
            "${OUT}/reuse.samples.sorted" \
            "${OUT}/expected.samples.sorted"
        then
            REUSE=1
        fi
    fi
fi

if [[ "${REUSE}" -eq 1 ]]; then

    echo "========================================"
    echo "Reusing existing PLINK genotype"
    echo "========================================"

    echo "Samples: ${OLD_N_SAMPLE}"
    echo "SNPs: ${OLD_N_SNP}"

else

    echo "========================================"
    echo "Preparing missing/invalid PLINK genotype"
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
        --out "${BFILE}"

fi

############################################
# 8. Final genotype checks
############################################

for EXT in bed bim fam
do
    test -s "${BFILE}.${EXT}" || {
        echo "ERROR: Missing ${BFILE}.${EXT}" >&2
        exit 1
    }
done

N_PEOPLE=$(wc -l < "${BFILE}.fam")
N_SNPS=$(wc -l < "${BFILE}.bim")

if [[ "${N_PEOPLE}" -ne "${N_EXPECTED}" ]]; then
    echo "ERROR: Genotype sample count mismatch" >&2
    exit 1
fi

if [[ "${N_SNPS}" -lt 100000 ]]; then
    echo "ERROR: Too few imputed SNPs" >&2
    exit 1
fi

BYTES_PER_SNP=$(( (N_PEOPLE + 3) / 4 ))
EXPECTED_BYTES=$(( N_SNPS * BYTES_PER_SNP + 3 ))
ACTUAL_BYTES=$(stat -c '%s' "${BFILE}.bed")

if [[ "${ACTUAL_BYTES}" -ne "${EXPECTED_BYTES}" ]]; then
    echo "ERROR: BED file size mismatch" >&2
    exit 1
fi

awk '
{
    if ($1 != $2) exit 1
    print $2
}
' "${BFILE}.fam" \
    | LC_ALL=C sort \
    > "${OUT}/genotype.samples.sorted"

LC_ALL=C sort "${TRAIN_SAMPLES}" \
    > "${OUT}/expected.samples.sorted"

if ! cmp -s \
    "${OUT}/genotype.samples.sorted" \
    "${OUT}/expected.samples.sorted"
then
    echo "ERROR: Genotype sample identities mismatch" >&2
    exit 1
fi

awk '{print $2}' "${GRM_PREFIX}.grm.id" \
    | LC_ALL=C sort \
    > "${OUT}/grm.samples.sorted"

if ! cmp -s \
    "${OUT}/genotype.samples.sorted" \
    "${OUT}/grm.samples.sorted"
then
    echo "ERROR: GRM/genotype sample identities mismatch" >&2
    exit 1
fi

echo "PASS: Genotype and GRM sample checks"

############################################
# 9. Prepare training-only phenotypes
############################################

for TRAIT in BB H6
do

    echo "========================================"
    echo "Preparing phenotype ${TRAIT}"
    echo "========================================"

    if [[ "${TRAIT}" == "BB" ]]; then

        PHENO_SOURCE="${WORK}/00_input/Trait-Adj_Budburst.txt"

    else

        PHENO_SOURCE="${WORK}/00_input/Trait-Adj_H6.txt"

    fi

    PHENO_OUT="${OUT}/${TRAIT}.pheno"

    PHENO_TMP="${PHENO_OUT}.tmp.$$"

    test -s "${PHENO_SOURCE}" || {
        echo "ERROR: Missing phenotype source" >&2
        exit 1
    }

    ########################################
    # IMPORTANT:
    #
    # 1. Remove CR (carriage return)
    #    if phenotype file has CRLF endings.
    #
    # 2. Avoid the previous fragile
    #    multiline regular-expression test.
    #
    # 3. Match phenotype by IID.
    #
    # 4. Output only inner_train individuals.
    ########################################

    awk '
    NR == FNR {

        if (NF != 1 || seen_train[$1]++) {
            print "ERROR: Invalid training IDs" > "/dev/stderr"
            bad = 1
            exit 1
        }

        wanted[$1] = 1
        next
    }

    {
        # Remove Windows carriage returns
        # before checking fields and values.
        sub(/\015+$/, "", $0)

        if (FNR == 1) {

            if ($1 != "FID" || $2 != "IID" || NF != 3) {
                print "ERROR: Invalid phenotype header" > "/dev/stderr"
                bad = 1
                exit 1
            }

            next
        }

        if (NF != 3) {
            print "ERROR: Invalid phenotype columns" > "/dev/stderr"
            bad = 1
            exit 1
        }

        if (!($2 in wanted)) {
            next
        }

        if ($1 != $2) {
            print "ERROR: Phenotype FID/IID mismatch", $2 > "/dev/stderr"
            bad = 1
            exit 1
        }

        if (seen_pheno[$2]++) {
            print "ERROR: Duplicate phenotype", $2 > "/dev/stderr"
            bad = 1
            exit 1
        }

        # Check missing phenotype codes.

        if ($3 == "" || $3 == "NA" || $3 == "NaN" || $3 == "-9" || $3 == "Inf" || $3 == "-Inf") {
            print "ERROR: Missing phenotype", $2, $3 > "/dev/stderr"
            bad = 1
            exit 1
        }

        # Validate numeric conversion without
        # relying on the previous regex.
        #
        # Normal numeric values, including
        # decimal and scientific notation,
        # should compare equal numerically.

        numeric_value = $3 + 0

        if ($3 != numeric_value) {
            print "ERROR: Non-numeric phenotype", $2, $3 > "/dev/stderr"
            bad = 1
            exit 1
        }

        print $1, $2, $3
    }

    END {

        if (bad) {
            exit 1
        }

        for (id in wanted) {

            if (seen_pheno[id] != 1) {
                print "ERROR: Missing phenotype for", id > "/dev/stderr"
                bad = 1
            }
        }

        if (bad) {
            exit 1
        }
    }
    ' "${TRAIN_SAMPLES}" "${PHENO_SOURCE}" \
        > "${PHENO_TMP}"

    N_PHENO=$(wc -l < "${PHENO_TMP}")

    if [[ "${N_PHENO}" -ne "${N_EXPECTED}" ]]; then
        echo "ERROR: Wrong ${TRAIT} phenotype count" >&2
        exit 1
    fi

    mv "${PHENO_TMP}" "${PHENO_OUT}"

    echo "PASS: ${TRAIT}, ${N_PHENO} samples"

done

############################################
# 10. Chromosome distribution
############################################

echo "Checking chromosome distribution"

awk '
{
    count[$1]++
}
END {
    for (chr = 1; chr <= 12; chr++) {
        print chr, count[chr] + 0
    }
}
' "${BFILE}.bim" \
    > "${OUT}/imputed_chr_counts.tsv"

N_AUTOSOME=$(awk '
$1 >= 1 && $1 <= 12 {
    n++
}
END {
    print n+0
}
' "${BFILE}.bim")

if [[ "${N_AUTOSOME}" -ne "${N_SNPS}" ]]; then
    echo "ERROR: Unexpected chromosome labels" >&2
    exit 1
fi

if awk '$2 == 0 {bad=1} END {exit !bad}' \
    "${OUT}/imputed_chr_counts.tsv"
then
    echo "ERROR: Empty chromosome" >&2
    exit 1
fi

############################################
# 11. Save summary
############################################

{
    printf "Outer\tInner\tN_train\tN_imputed_QC\n"

    printf "%s\t%s\t%s\t%s\n" \
        "${OUTER}" \
        "${INNER}" \
        "${N_PEOPLE}" \
        "${N_SNPS}"

} > "${OUT}/input_summary.tsv"

############################################
# 12. Completion marker
############################################

printf "DONE Step8.2b outer_%s inner_%s\n" \
    "${OUTER}" "${INNER}" \
    > "${OUT}/complete.flag"

echo "========================================"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "Training samples: ${N_PEOPLE}"
echo "Imputed SNPs after QC: ${N_SNPS}"
echo "DONE Step8.2b"
echo "========================================"

