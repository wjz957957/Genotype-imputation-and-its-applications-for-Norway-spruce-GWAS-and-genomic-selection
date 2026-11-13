#!/bin/bash
#SBATCH -J 08_chip_GRM
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-29%5
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_build_chip_grm_%A_%a.out
#SBATCH -e logs/08_build_chip_grm_%A_%a.err

set -euo pipefail

########################################
# 1. Modules
########################################

module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

if command -v gcta64 >/dev/null 2>&1; then
    GCTA=gcta64
elif command -v gcta >/dev/null 2>&1; then
    GCTA=gcta
else
    echo "ERROR: GCTA executable not found" >&2
    exit 1
fi

########################################
# 2. Paths
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"

CHIP_LIST="${WORK}/01_reference/chip.marker"

CV_DIR="${WORK}/02_CV"

OUTROOT="${WORK}/08_imputed_GWAS/01_chip_GRM"

########################################
# 3. Task mapping
########################################

TASK=${SLURM_ARRAY_TASK_ID}

OUTER=$(printf "%02d" $((TASK / 3 + 1)))
INNER=$(printf "%02d" $((TASK % 3 + 1)))

echo "========================================"
echo "Step 8.2a: Chip GRM"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

########################################
# 4. Inputs and outputs
########################################

TRAIN_SAMPLES="${CV_DIR}/outer_${OUTER}/inner_${INNER}/inner_train.samples"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}"

PLINK_DIR="${OUT}/plink"
GRM_DIR="${OUT}/grm"

mkdir -p "${PLINK_DIR}" "${GRM_DIR}"

CHIP_BFILE="${PLINK_DIR}/chip_train_qc"

GRM_ALL="${GRM_DIR}/chip_all"

rm -f "${OUT}/complete.flag"

########################################
# 5. Validate inputs
########################################

for EXT in bed bim fam
do
    test -s "${MASTER}.${EXT}" || {
        echo "ERROR: Missing master .${EXT}" >&2
        exit 1
    }
done

test -s "${CHIP_LIST}"
test -s "${TRAIN_SAMPLES}"

N_EXPECTED=$(wc -l < "${TRAIN_SAMPLES}")

echo "Expected training samples: ${N_EXPECTED}"

if [[ "${N_EXPECTED}" -lt 500 ||
      "${N_EXPECTED}" -gt 600 ]]; then
    echo "ERROR: Unexpected training sample count" >&2
    exit 1
fi

########################################
# 6. Create PLINK --keep
########################################

KEEP="${OUT}/inner_train.keep"

awk '
NF != 1 {
    print "ERROR: sample list must have one column" > "/dev/stderr"
    exit 1
}
{
    if (seen[$1]++) {
        print "ERROR: duplicate sample ID" > "/dev/stderr"
        exit 1
    }
    print $1, $1
}
' "${TRAIN_SAMPLES}" > "${KEEP}"

########################################
# 7. Extract chip SNPs and inner_train
########################################

echo "Extracting chip genotypes"

# QC is applied after restricting to inner_train.
# The output contains only training samples.
#
# --maf 0.01 : retain MAF >= 0.01
# --geno 0.05: retain SNP missingness <= 5%
#
# These thresholds are fixed in advance.

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
    --memory 16000 \
    --out "${CHIP_BFILE}"

########################################
# 8. Validate samples
########################################

N_TRAIN=$(wc -l < "${CHIP_BFILE}.fam")
N_CHIP=$(wc -l < "${CHIP_BFILE}.bim")

echo "Retained samples: ${N_TRAIN}"
echo "Retained chip SNPs after QC: ${N_CHIP}"

if [[ "${N_TRAIN}" -ne "${N_EXPECTED}" ]]; then
    echo "ERROR: Incorrect sample count" >&2
    exit 1
fi

if [[ "${N_CHIP}" -lt 10000 ]]; then
    echo "ERROR: Too few chip SNPs after QC" >&2
    exit 1
fi

awk '{print $2}' "${CHIP_BFILE}.fam" \
    | LC_ALL=C sort > "${OUT}/observed.samples.sorted"

LC_ALL=C sort "${TRAIN_SAMPLES}" \
    > "${OUT}/expected.samples.sorted"

if ! cmp -s \
    "${OUT}/observed.samples.sorted" \
    "${OUT}/expected.samples.sorted"
then
    echo "ERROR: Training sample identities mismatch" >&2
    exit 1
fi

echo "PASS: Training samples match CV partition"

########################################
# 9. Count SNPs by chromosome
########################################

awk '
{
    count[$1]++
}
END {
    for (chr = 1; chr <= 12; chr++) {
        print chr, count[chr] + 0
    }
}
' "${CHIP_BFILE}.bim" \
    > "${OUT}/chip_chr_counts.tsv"

if awk '$2 == 0 {bad=1} END {exit !bad}' \
    "${OUT}/chip_chr_counts.tsv"
then
    echo "ERROR: At least one chromosome has zero SNPs" >&2
    exit 1
fi

N_AUTOSOME=$(awk '
$1 >= 1 && $1 <= 12 {n++}
END {print n+0}
' "${CHIP_BFILE}.bim")

if [[ "${N_AUTOSOME}" -ne "${N_CHIP}" ]]; then
    echo "ERROR: Non-autosomal SNPs detected" >&2
    exit 1
fi

########################################
# 10. Whole-chip GRM
########################################

echo "========================================"
echo "Building GRM from all chip SNPs"
echo "========================================"

"${GCTA}" \
    --bfile "${CHIP_BFILE}" \
    --autosome-num 12 \
    --make-grm \
    --thread-num 4 \
    --out "${GRM_ALL}"

for EXT in grm.bin grm.N.bin grm.id
do
    test -s "${GRM_ALL}.${EXT}" || {
        echo "ERROR: Missing full GRM ${EXT}" >&2
        exit 1
    }
done

if [[ $(wc -l < "${GRM_ALL}.grm.id") -ne "${N_TRAIN}" ]]; then
    echo "ERROR: Full GRM sample count mismatch" >&2
    exit 1
fi

########################################
# 11. Chromosome-specific GRMs
########################################

for CHR in $(seq 1 12)
do

    CHR_PAD=$(printf "%02d" "${CHR}")

    GRM_CHR="${GRM_DIR}/chip_chr${CHR_PAD}"

    echo "========================================"
    echo "Building chip GRM: chromosome ${CHR}"
    echo "========================================"

    "${GCTA}" \
        --bfile "${CHIP_BFILE}" \
        --autosome-num 12 \
        --chr "${CHR}" \
        --make-grm \
        --thread-num 4 \
        --out "${GRM_CHR}"

    for EXT in grm.bin grm.N.bin grm.id
    do
        test -s "${GRM_CHR}.${EXT}" || {
            echo "ERROR: Missing chr${CHR} GRM ${EXT}" >&2
            exit 1
        }
    done

    if [[ $(wc -l < "${GRM_CHR}.grm.id") -ne "${N_TRAIN}" ]]; then
        echo "ERROR: Chromosome GRM sample mismatch" >&2
        exit 1
    fi

    # All GRMs must contain the identical individuals,
    # in the identical order.

    if ! cmp -s \
        "${GRM_ALL}.grm.id" \
        "${GRM_CHR}.grm.id"
    then
        echo "ERROR: Chromosome GRM IDs/order mismatch" >&2
        exit 1
    fi

    echo "PASS chromosome ${CHR}"

done

########################################
# 12. Save summary
########################################

SUMMARY="${OUT}/grm_summary.tsv"

{
    printf "Outer\tInner\tN_train\tN_chip_QC\tN_chromosomes\n"
    printf "%s\t%s\t%s\t%s\t12\n" \
        "${OUTER}" "${INNER}" "${N_TRAIN}" "${N_CHIP}"
} > "${SUMMARY}"

touch "${OUT}/complete.flag"

echo "========================================"
echo "Outer ${OUTER} Inner ${INNER}"
echo "Training samples: ${N_TRAIN}"
echo "Chip SNPs after QC: ${N_CHIP}"
echo "GRMs produced: 13"
echo "DONE Step8.2a"
echo "========================================"

