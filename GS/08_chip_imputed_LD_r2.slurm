#!/bin/bash
#SBATCH -J 08_chip_imp_r2
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-359%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_chip_imputed_LD_r2_%A_%a.out
#SBATCH -e logs/08_chip_imputed_LD_r2_%A_%a.err

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

########################################
# 1. Software
########################################

module load bioinfo-tools
ml plink/1.90b4.9

########################################
# 2. Paths
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"

CHIP_ROOT="${WORK}/08_imputed_GWAS/01_chip_GRM"

IMP_ROOT="${WORK}/08_imputed_GWAS/02_imputed_inputs"

CV_ROOT="${WORK}/02_CV"

OUTROOT="${WORK}/08_imputed_GWAS/05_chip_imputed_LD_r2"

########################################
# 3. Parameters
########################################

R2_THRESHOLD=0.1
WINDOW_KB=1000
WINDOW_BP=1000000

########################################
# 4. Task mapping
########################################

TASK=${SLURM_ARRAY_TASK_ID}

FOLD_INDEX=$(( TASK / 12 ))
CHR=$(( TASK % 12 + 1 ))

OUTER=$(printf "%02d" $(( FOLD_INDEX / 3 + 1 )))
INNER=$(printf "%02d" $(( FOLD_INDEX % 3 + 1 )))
CHR_PAD=$(printf "%02d" "${CHR}")

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}/chr_${CHR_PAD}"

mkdir -p "${OUT}"

echo "========================================"
echo "Step 8.3a: Direct pairwise chip-imputed LD"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "Chromosome: ${CHR}"
echo "R2 threshold: ${R2_THRESHOLD}"
echo "Window: ${WINDOW_KB} kb"
echo "========================================"

########################################
# 5. Input/output files
########################################

TRAIN="${CV_ROOT}/outer_${OUTER}/inner_${INNER}/inner_train.samples"

CHIP_BIM="${CHIP_ROOT}/outer_${OUTER}/inner_${INNER}/plink/chip_train_qc.bim"

IMP_DIR="${IMP_ROOT}/outer_${OUTER}/inner_${INNER}"

IMP_BIM="${IMP_DIR}/imputed_train_qc.bim"

CHIP_LIST="${OUT}/chip.snps"
IMP_LIST="${OUT}/imputed.snps"
UNION_LIST="${OUT}/chip_plus_imputed.snps"

KEEP="${OUT}/train.keep"

PREFIX="${OUT}/chip_imputed_r2"
LD_FILE="${PREFIX}.ld.gz"

EXCLUDED="${OUT}/chip_ld_excluded.snps"
INDEPENDENT="${OUT}/chip_independent_imputed.snps"

SUMMARY="${OUT}/ld_summary.tsv"
FLAG="${OUT}/complete.flag"

########################################
# 6. Skip completed task
########################################

if [[ -s "${FLAG}" && -s "${SUMMARY}" \
      && -f "${EXCLUDED}" && -f "${INDEPENDENT}" ]]; then
    echo "SKIP: Previously completed"
    exit 0
fi

rm -f "${FLAG}"

########################################
# 7. Validate inputs
########################################

for EXT in bed bim fam
do
    if [[ ! -s "${MASTER}.${EXT}" ]]; then
        echo "ERROR: Missing master ${EXT}" >&2
        exit 1
    fi
done

for F in "${TRAIN}" "${CHIP_BIM}" "${IMP_BIM}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: Missing input ${F}" >&2
        exit 1
    fi
done

if [[ ! -s "${IMP_DIR}/complete.flag" ]]; then
    echo "ERROR: Step 8.2b incomplete" >&2
    exit 1
fi

########################################
# 8. Training individuals
########################################

awk '
NF != 1 {
    print "ERROR: Invalid training sample list" > "/dev/stderr"
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
    echo "ERROR: Unexpected sample count" >&2
    exit 1
fi

echo "Training individuals: ${N_TRAIN}"

########################################
# 9. Fold-specific QC SNP lists
########################################

awk -v chr="${CHR}" '
$1 == chr { print $2 }
' "${CHIP_BIM}" > "${CHIP_LIST}"

awk -v chr="${CHR}" '
$1 == chr { print $2 }
' "${IMP_BIM}" > "${IMP_LIST}"

N_CHIP=$(wc -l < "${CHIP_LIST}")
N_IMP=$(wc -l < "${IMP_LIST}")

echo "Chip SNPs: ${N_CHIP}"
echo "Imputed SNPs: ${N_IMP}"

if (( N_CHIP < 100 || N_IMP < 1000 )); then
    echo "ERROR: Unexpected SNP counts" >&2
    exit 1
fi

########################################
# 10. Include both SNP sets
########################################

cat "${CHIP_LIST}" "${IMP_LIST}" \
    | LC_ALL=C sort -u \
    > "${UNION_LIST}"

N_UNION=$(wc -l < "${UNION_LIST}")

if (( N_UNION != N_CHIP + N_IMP )); then
    echo "ERROR: Duplicate or overlapping SNP IDs" >&2
    exit 1
fi

echo "Total variants: ${N_UNION}"

########################################
# 11. Actual pairwise r2 calculation
########################################

echo "========================================"
echo "Running PLINK --r2 gz"
echo "========================================"

# --extract retains BOTH chip and imputed.
#
# --ld-snp-list restricts one side of the
# reported pair to chip reference SNPs.
#
# --ld-window 999999 avoids default
# 10-variant window truncation.
#
# --ld-window-kb 1000 imposes 1 Mb.
#
# --ld-window-r2 0.1 outputs pairs with
# r2 >= 0.1.
#
# --r2 gz compresses the pair report.

plink \
    --bfile "${MASTER}" \
    --keep "${KEEP}" \
    --chr "${CHR}" \
    --extract "${UNION_LIST}" \
    --ld-snp-list "${CHIP_LIST}" \
    --r2 gz \
    --ld-window 999999 \
    --ld-window-kb "${WINDOW_KB}" \
    --ld-window-r2 "${R2_THRESHOLD}" \
    --threads 4 \
    --memory 24000 \
    --out "${PREFIX}"

########################################
# 12. Check LD output
########################################

if [[ ! -s "${LD_FILE}" ]]; then
    echo "ERROR: Missing LD report" >&2
    exit 1
fi

gzip -t "${LD_FILE}"

echo "PASS: LD file exists and gzip integrity OK"

########################################
# 13. Parse actual chip-imputed pairs
########################################

# Check the LD table header, validate
# reference/candidate SNP membership,
# and verify the reported physical
# distance and r2.
#
# Only imputed SNPs observed in a
# chip-imputed pair are excluded.
#
# Chip-chip pairs are ignored.

gzip -dc "${LD_FILE}" | \
awk \
-v chipfile="${CHIP_LIST}" \
-v impfile="${IMP_LIST}" \
-v maxbp="${WINDOW_BP}" \
-v threshold="${R2_THRESHOLD}" \
'
BEGIN {
    while ((getline line < chipfile) > 0) {
        split(line, f)
        chip[f[1]] = 1
    }
    close(chipfile)

    while ((getline line < impfile) > 0) {
        split(line, f)
        imp[f[1]] = 1
    }
    close(impfile)
}

NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }

    if (!("SNP_A" in col) ||
        !("SNP_B" in col) ||
        !("BP_A" in col) ||
        !("BP_B" in col) ||
        !("R2" in col)) {
        print "ERROR: Unexpected LD header" > "/dev/stderr"
        bad = 1
        exit 1
    }

    next
}

{
    a = $(col["SNP_A"])
    b = $(col["SNP_B"])

    bp_a = $(col["BP_A"]) + 0
    bp_b = $(col["BP_B"]) + 0
    r2 = $(col["R2"]) + 0

    if (r2 < 0 || r2 > 1) {
        print "ERROR: Invalid r2" > "/dev/stderr"
        bad = 1
        exit 1
    }

    distance = bp_a - bp_b
    if (distance < 0) distance = -distance

    if (distance > maxbp) {
        print "ERROR: LD pair outside window" > "/dev/stderr"
        bad = 1
        exit 1
    }

    if (r2 + 0.000001 < threshold) {
        print "ERROR: LD pair below threshold" > "/dev/stderr"
        bad = 1
        exit 1
    }

    if ((a in chip) && (b in imp)) {
        excluded[b] = 1
        cross_pairs++
    }
    else if ((b in chip) && (a in imp)) {
        excluded[a] = 1
        cross_pairs++
    }
    else if ((a in chip) && (b in chip)) {
        chip_pairs++
    }
    else {
        print "ERROR: Unexpected LD pair:", a, b > "/dev/stderr"
        bad = 1
        exit 1
    }
}

END {
    if (NR == 0 || bad) {
        exit 1
    }

    for (snp in excluded) {
        print snp
    }

    print "Chip-imputed LD pairs:", cross_pairs+0 > "/dev/stderr"
    print "Chip-chip LD pairs:", chip_pairs+0 > "/dev/stderr"
}
' > "${OUT}/excluded.unsorted"

LC_ALL=C sort -u "${OUT}/excluded.unsorted" \
    > "${EXCLUDED}"

########################################
# 14. Keep imputed SNPs not excluded
########################################

awk '
NR == FNR {
    excluded[$1] = 1
    next
}
!($1 in excluded) {
    print $1
}
' "${EXCLUDED}" "${IMP_LIST}" \
    > "${INDEPENDENT}"

########################################
# 15. Counts
########################################

N_EXCLUDED=$(wc -l < "${EXCLUDED}")
N_INDEPENDENT=$(wc -l < "${INDEPENDENT}")

echo "========================================"
echo "LD filtering summary"
echo "========================================"

echo "Training samples: ${N_TRAIN}"
echo "Chip SNPs: ${N_CHIP}"
echo "Imputed SNPs: ${N_IMP}"
echo "Excluded by actual chip-imputed LD: ${N_EXCLUDED}"
echo "Remaining imputed SNPs: ${N_INDEPENDENT}"

if (( N_EXCLUDED + N_INDEPENDENT != N_IMP )); then
    echo "ERROR: SNP count mismatch" >&2
    exit 1
fi

########################################
# 16. Check partition
########################################

LC_ALL=C sort "${IMP_LIST}" \
    > "${OUT}/imputed.original.sorted"

cat "${EXCLUDED}" "${INDEPENDENT}" \
    | LC_ALL=C sort \
    > "${OUT}/imputed.partition.sorted"

if ! cmp -s \
    "${OUT}/imputed.original.sorted" \
    "${OUT}/imputed.partition.sorted"
then
    echo "ERROR: SNP partition mismatch" >&2
    exit 1
fi

echo "PASS: SNP partition"

########################################
# 17. Summary
########################################

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
        "${R2_THRESHOLD}" \
        "${WINDOW_KB}"

} > "${SUMMARY}"

########################################
# 18. Require review if all excluded
########################################

if (( N_INDEPENDENT == 0 )); then
    echo "WARNING: All imputed SNPs excluded by direct pairwise r2" >&2
    echo "WARNING: Result requires investigation before proceeding" >&2
    exit 1
fi

########################################
# 19. Complete
########################################

printf "DONE Step8.3a outer_%s inner_%s chr_%s\n" \
    "${OUTER}" "${INNER}" "${CHR_PAD}" \
    > "${FLAG}"

echo "========================================"
echo "DONE Step8.3a"
echo "========================================"

