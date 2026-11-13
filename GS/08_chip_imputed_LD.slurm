#!/bin/bash
#SBATCH -J 08_chip_imp_LD
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-359%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_chip_imputed_LD_%A_%a.out
#SBATCH -e logs/08_chip_imputed_LD_%A_%a.err

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

############################################
# 1. Software
############################################

module load bioinfo-tools
ml plink/1.90b4.9

############################################
# 2. Paths
############################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"

CHIP_ROOT="${WORK}/08_imputed_GWAS/01_chip_GRM"

IMP_ROOT="${WORK}/08_imputed_GWAS/02_imputed_inputs"

CV_ROOT="${WORK}/02_CV"

OUTROOT="${WORK}/08_imputed_GWAS/05_chip_imputed_LD"

############################################
# 3. LD parameters
############################################

TAG_R2=0.1
TAG_KB=1000

############################################
# 4. Array mapping
############################################

TASK=${SLURM_ARRAY_TASK_ID}

FOLD_INDEX=$((TASK / 12))
CHR=$((TASK % 12 + 1))

OUTER=$(printf "%02d" $((FOLD_INDEX / 3 + 1)))
INNER=$(printf "%02d" $((FOLD_INDEX % 3 + 1)))
CHR_PAD=$(printf "%02d" "${CHR}")

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}/chr_${CHR_PAD}"

mkdir -p "${OUT}"

echo "========================================"
echo "Step 8.3a: Full chip vs imputed LD"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "Chromosome: ${CHR}"
echo "R2 threshold: ${TAG_R2}"
echo "Window: ${TAG_KB} kb"
echo "========================================"

############################################
# 5. Files
############################################

TRAIN="${CV_ROOT}/outer_${OUTER}/inner_${INNER}/inner_train.samples"

CHIP_BIM="${CHIP_ROOT}/outer_${OUTER}/inner_${INNER}/plink/chip_train_qc.bim"

IMP_DIR="${IMP_ROOT}/outer_${OUTER}/inner_${INNER}"

IMP_BIM="${IMP_DIR}/imputed_train_qc.bim"

CHIP_LIST="${OUT}/chip.snps"
IMP_LIST="${OUT}/imputed.snps"
UNION_LIST="${OUT}/chip_plus_imputed.snps"

KEEP="${OUT}/train.keep"

PREFIX="${OUT}/chip_vs_imputed"

EXCLUDED="${OUT}/chip_ld_excluded.snps"
INDEPENDENT="${OUT}/chip_independent_imputed.snps"

SUMMARY="${OUT}/ld_summary.tsv"
FLAG="${OUT}/complete.flag"

############################################
# 6. Skip successfully completed tasks
############################################

if [[ -s "${FLAG}" && -s "${SUMMARY}" && -f "${INDEPENDENT}" ]]; then
    echo "SKIP: Already completed"
    echo "DONE Step8.3a (existing)"
    exit 0
fi

rm -f "${FLAG}"

############################################
# 7. Input validation
############################################

for EXT in bed bim fam
do
    [[ -s "${MASTER}.${EXT}" ]] || {
        echo "ERROR: Missing master ${EXT}" >&2
        exit 1
    }
done

for F in "${TRAIN}" "${CHIP_BIM}" "${IMP_BIM}"
do
    [[ -s "${F}" ]] || {
        echo "ERROR: Missing input ${F}" >&2
        exit 1
    }
done

[[ -s "${IMP_DIR}/complete.flag" ]] || {
    echo "ERROR: Step 8.2b incomplete" >&2
    exit 1
}

############################################
# 8. Training samples
############################################

# Master genotype uses FID = IID.
# The CV sample files contain one sample ID
# per line.

awk '
NF != 1 {
    print "ERROR: Invalid sample list" > "/dev/stderr"
    exit 1
}
{
    if (seen[$1]++) {
        print "ERROR: Duplicate sample ID" > "/dev/stderr"
        exit 1
    }
    print $1, $1
}
' "${TRAIN}" > "${KEEP}"

N_TRAIN=$(wc -l < "${KEEP}")

if (( N_TRAIN < 500 || N_TRAIN > 600 )); then
    echo "ERROR: Unexpected training sample count ${N_TRAIN}" >&2
    exit 1
fi

echo "Training samples: ${N_TRAIN}"

############################################
# 9. All QC-passing chip SNPs
############################################

awk -v chr="${CHR}" '
$1 == chr {
    print $2
}
' "${CHIP_BIM}" > "${CHIP_LIST}"

N_CHIP=$(wc -l < "${CHIP_LIST}")

if (( N_CHIP < 100 )); then
    echo "ERROR: Too few chip SNPs" >&2
    exit 1
fi

echo "Chip SNPs: ${N_CHIP}"

############################################
# 10. All QC-passing imputed SNPs
############################################

awk -v chr="${CHR}" '
$1 == chr {
    print $2
}
' "${IMP_BIM}" > "${IMP_LIST}"

N_IMP=$(wc -l < "${IMP_LIST}")

if (( N_IMP < 1000 )); then
    echo "ERROR: Too few imputed SNPs" >&2
    exit 1
fi

echo "Imputed SNPs: ${N_IMP}"

############################################
# 11. Load BOTH SNP sets
############################################

# CRITICAL FIX:
#
# Previous failed script loaded only imputed
# SNPs with --extract.
#
# Here we explicitly create the union of
# all chip and all imputed SNPs.
#
# All chip SNPs remain available as
# --show-tags reference variants.

cat "${CHIP_LIST}" "${IMP_LIST}" \
    | LC_ALL=C sort -u \
    > "${UNION_LIST}"

N_UNION=$(wc -l < "${UNION_LIST}")

if (( N_UNION != N_CHIP + N_IMP )); then
    echo "ERROR: Chip/imputed SNP overlap or duplicate IDs" >&2
    exit 1
fi

echo "Total SNPs loaded for LD: ${N_UNION}"

############################################
# 12. Full chip vs full imputed LD
############################################

echo "========================================"
echo "Running PLINK --show-tags"
echo "========================================"

plink \
    --bfile "${MASTER}" \
    --keep "${KEEP}" \
    --chr "${CHR}" \
    --extract "${UNION_LIST}" \
    --show-tags "${CHIP_LIST}" \
    --tag-r2 "${TAG_R2}" \
    --tag-kb "${TAG_KB}" \
    --threads 4 \
    --memory 24000 \
    --out "${PREFIX}"

############################################
# 13. Validate PLINK tag output
############################################

TAG_FILE="${PREFIX}.tags"

if [[ ! -s "${TAG_FILE}" ]]; then
    echo "ERROR: No valid PLINK .tags output" >&2
    exit 1
fi

# Default --show-tags output is a
# one-column list of SNP IDs.

awk '
NF != 1 {
    print "ERROR: Invalid .tags format" > "/dev/stderr"
    exit 1
}
' "${TAG_FILE}"

############################################
# 14. Verify all chip references
############################################

# PLINK includes each reference SNP
# as its own tag.
#
# Therefore every chip SNP must appear
# in the output.

awk '
NR == FNR {
    tagged[$1] = 1
    next
}
!($1 in tagged) {
    print "ERROR: Missing chip reference:", $1 > "/dev/stderr"
    bad = 1
}
END {
    if (bad) exit 1
}
' "${TAG_FILE}" "${CHIP_LIST}"

echo "PASS: All chip SNPs used as LD references"

############################################
# 15. Imputed SNPs excluded by chip LD
############################################

# Retain ONLY imputed SNPs appearing
# in the chip tag list.
#
# These have r2 >= 0.1 with at least
# one chip SNP within the specified window.

awk '
NR == FNR {
    tagged[$1] = 1
    next
}
($1 in tagged) {
    print $1
}
' "${TAG_FILE}" "${IMP_LIST}" \
    > "${EXCLUDED}"

############################################
# 16. Chip-independent imputed SNPs
############################################

awk '
NR == FNR {
    tagged[$1] = 1
    next
}
!($1 in tagged) {
    print $1
}
' "${TAG_FILE}" "${IMP_LIST}" \
    > "${INDEPENDENT}"

############################################
# 17. Check SNP counts
############################################

N_EXCLUDED=$(wc -l < "${EXCLUDED}")

N_INDEPENDENT=$(wc -l < "${INDEPENDENT}")

echo "========================================"
echo "LD filtering summary"
echo "========================================"

echo "Chip SNPs: ${N_CHIP}"
echo "Imputed SNPs: ${N_IMP}"
echo "Excluded imputed SNPs: ${N_EXCLUDED}"
echo "Remaining imputed SNPs: ${N_INDEPENDENT}"

if (( N_EXCLUDED + N_INDEPENDENT != N_IMP )); then
    echo "ERROR: SNP count mismatch" >&2
    exit 1
fi

if (( N_INDEPENDENT == 0 )); then
    echo "ERROR: No remaining imputed SNPs" >&2
    exit 1
fi

############################################
# 18. Verify exact partition
############################################

LC_ALL=C sort "${IMP_LIST}" \
    > "${OUT}/imputed.original.sorted"

cat "${EXCLUDED}" "${INDEPENDENT}" \
    | LC_ALL=C sort \
    > "${OUT}/imputed.filtered.sorted"

if ! cmp -s \
    "${OUT}/imputed.original.sorted" \
    "${OUT}/imputed.filtered.sorted"
then
    echo "ERROR: SNP partition verification failed" >&2
    exit 1
fi

echo "PASS: SNP partition verification"

############################################
# 19. Save QC summary
############################################

{
    printf "Outer\tInner\tChr\tN_train\tChip_SNPs\tImputed_SNPs\tExcluded_by_chip_LD\tIndependent_imputed\tR2\tWindow_kb\n"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${OUTER}" \
        "${INNER}" \
        "${CHR}" \
        "${N_TRAIN}" \
        "${N_CHIP}" \
        "${N_IMP}" \
        "${N_EXCLUDED}" \
        "${N_INDEPENDENT}" \
        "${TAG_R2}" \
        "${TAG_KB}"

} > "${SUMMARY}"

############################################
# 20. Completion flag
############################################

printf "DONE Step8.3a outer_%s inner_%s chr_%s\n" \
    "${OUTER}" "${INNER}" "${CHR_PAD}" \
    > "${FLAG}"

echo "========================================"
echo "DONE Step8.3a"
echo "========================================"

