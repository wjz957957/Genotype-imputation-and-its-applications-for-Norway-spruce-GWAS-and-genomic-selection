#!/bin/bash
#SBATCH -J cf_chip_prep
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-14%3
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/01_chip_inner_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/01_chip_inner_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"

SAMPLE_ROOT="${WORK}/Cross_family/cv_samples_v1"
MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
CHIP_MARKERS="${WORK}/01_reference/chip.marker"
OUTROOT="${WORK}/Cross_family/results/01_chip_inner_genotypes"

N_EXPECTED_SNP=43267
SCRIPT_VERSION=1

########################################
# 1. Software
########################################

module load bioinfo-tools
ml plink/1.90b4.9

########################################
# 2. Map Slurm task to outer/inner fold
########################################

TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

if [[ ! "${TASK}" =~ ^[0-9]+$ ]] || (( TASK < 0 || TASK > 14 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 0 and 14" >&2
    exit 1
fi

OUTER=$(printf '%02d' "$(( TASK / 3 + 1 ))")
INNER=$(printf '%02d' "$(( TASK % 3 + 1 ))")

TRAIN="${SAMPLE_ROOT}/outer_${OUTER}/inner_${INNER}/inner_train.samples"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}"
PREFIX="${OUT}/chip_inner_train"
KEEP="${OUT}/inner_train.keep"
FLAG="${OUT}/complete.flag"

mkdir -p "${OUT}"

echo "========================================"
echo "Cross-family chip genotype preparation"
echo "Task:  ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

########################################
# 3. Check input files
########################################

for EXT in bed bim fam; do
    [[ -s "${MASTER}.${EXT}" ]] || {
        echo "ERROR: Missing ${MASTER}.${EXT}" >&2
        exit 1
    }
done

[[ -s "${TRAIN}" ]] || {
    echo "ERROR: Missing training samples: ${TRAIN}" >&2
    exit 1
}

[[ -s "${CHIP_MARKERS}" ]] || {
    echo "ERROR: Missing chip marker list: ${CHIP_MARKERS}" >&2
    exit 1
}

########################################
# 4. Validate master samples and markers
########################################

N_MASTER=$(wc -l < "${MASTER}.fam")

if (( N_MASTER != 904 )); then
    echo "ERROR: Expected 904 master samples; found ${N_MASTER}" >&2
    exit 1
fi

awk '
{
    if (NF < 2 || $1 != $2 || seen[$2]++) {
        print "ERROR: Invalid or duplicate master sample ID at line", NR > "/dev/stderr"
        bad=1
    }
}
END { if (bad) exit 1 }
' "${MASTER}.fam"

N_MARKER_LINES=$(wc -l < "${CHIP_MARKERS}")
N_MARKER_UNIQUE=$(sort -u "${CHIP_MARKERS}" | wc -l)

if (( N_MARKER_LINES != N_EXPECTED_SNP ||
      N_MARKER_UNIQUE != N_EXPECTED_SNP )); then
    echo "ERROR: Chip marker count/uniqueness check failed" >&2
    echo "Lines=${N_MARKER_LINES}; unique=${N_MARKER_UNIQUE}" >&2
    exit 1
fi

########################################
# 5. Validate inner-training IDs
########################################

awk '
{
    if (NF != 1 || seen[$1]++) {
        print "ERROR: Invalid or duplicate training ID at line", NR > "/dev/stderr"
        bad=1
    }
}
END {
    if (NR == 0 || bad) exit 1
}
' "${TRAIN}"

N_TRAIN=$(wc -l < "${TRAIN}")

# Every training ID must occur in the master PLINK FAM.
awk '
NR == FNR {
    master[$2]=1
    next
}
!($1 in master) {
    print "ERROR: Training ID absent from master:", $1 > "/dev/stderr"
    bad=1
}
END { if (bad) exit 1 }
' "${MASTER}.fam" "${TRAIN}"

echo "Master samples: ${N_MASTER}"
echo "Inner training samples: ${N_TRAIN}"
echo "Expected chip SNPs: ${N_EXPECTED_SNP}"

########################################
# 6. Skip a previously verified result
########################################

if [[ -s "${FLAG}" &&
      -s "${PREFIX}.bed" &&
      -s "${PREFIX}.bim" &&
      -s "${PREFIX}.fam" ]] &&
   grep -qx "version=${SCRIPT_VERSION}" "${FLAG}"; then

    echo "Previously completed output detected."
    echo "Skipping extraction; existing output was not modified."
    exit 0
fi

########################################
# 7. Prepare PLINK keep file
########################################

# The master FID and IID are identical.
awk '{print $1, $1}' "${TRAIN}" > "${KEEP}"

if [[ $(wc -l < "${KEEP}") -ne ${N_TRAIN} ]]; then
    echo "ERROR: Keep-file row count mismatch" >&2
    exit 1
fi

########################################
# 8. Extract training individuals and chip SNPs
########################################

# Remove stale completion status before starting.
rm -f "${FLAG}"

plink \
    --bfile "${MASTER}" \
    --keep "${KEEP}" \
    --extract "${CHIP_MARKERS}" \
    --make-bed \
    --allow-extra-chr \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --memory 24000 \
    --out "${PREFIX}"

########################################
# 9. Validate extracted dimensions
########################################

for EXT in bed bim fam; do
    [[ -s "${PREFIX}.${EXT}" ]] || {
        echo "ERROR: Missing extracted ${EXT} file" >&2
        exit 1
    }
done

N_OUT_TRAIN=$(wc -l < "${PREFIX}.fam")
N_OUT_SNP=$(wc -l < "${PREFIX}.bim")

if (( N_OUT_TRAIN != N_TRAIN )); then
    echo "ERROR: Training sample count mismatch" >&2
    echo "Expected=${N_TRAIN}; observed=${N_OUT_TRAIN}" >&2
    exit 1
fi

if (( N_OUT_SNP != N_EXPECTED_SNP )); then
    echo "ERROR: SNP count mismatch" >&2
    echo "Expected=${N_EXPECTED_SNP}; observed=${N_OUT_SNP}" >&2
    exit 1
fi

########################################
# 10. Exact sample-ID comparison
########################################

EXPECTED_IDS="${OUT}/expected_sample_ids.sorted"
OBSERVED_IDS="${OUT}/observed_sample_ids.sorted"

sort "${TRAIN}" > "${EXPECTED_IDS}"
awk '{print $2}' "${PREFIX}.fam" | sort > "${OBSERVED_IDS}"

awk '
$1 != $2 {bad=1}
END {if (bad) exit 1}
' "${PREFIX}.fam"

if ! cmp -s "${EXPECTED_IDS}" "${OBSERVED_IDS}"; then
    echo "ERROR: Extracted sample IDs differ from inner training IDs" >&2
    diff -u "${EXPECTED_IDS}" "${OBSERVED_IDS}" | head -n 30 || true
    exit 1
fi

########################################
# 11. Exact SNP-ID comparison
########################################

EXPECTED_SNPS="${OUT}/expected_chip_snps.sorted"
OBSERVED_SNPS="${OUT}/observed_chip_snps.sorted"

sort -u "${CHIP_MARKERS}" > "${EXPECTED_SNPS}"
awk '{print $2}' "${PREFIX}.bim" | sort > "${OBSERVED_SNPS}"

if ! cmp -s "${EXPECTED_SNPS}" "${OBSERVED_SNPS}"; then
    echo "ERROR: Extracted SNP IDs differ from chip.marker" >&2
    diff -u "${EXPECTED_SNPS}" "${OBSERVED_SNPS}" | head -n 30 || true
    exit 1
fi

########################################
# 12. Publish completion flag
########################################

{
    echo "version=${SCRIPT_VERSION}"
    echo "outer=${OUTER}"
    echo "inner=${INNER}"
    echo "n_train=${N_TRAIN}"
    echo "n_chip_snps=${N_OUT_SNP}"
    echo "sample_ID_check=PASS"
    echo "SNP_ID_check=PASS"
} > "${FLAG}.tmp"

mv "${FLAG}.tmp" "${FLAG}"

echo "========================================"
echo "PASS: outer_${OUTER}/inner_${INNER}"
echo "Training samples: ${N_OUT_TRAIN}"
echo "Chip SNPs: ${N_OUT_SNP}"
echo "Sample IDs: exact match"
echo "SNP IDs: exact match"
echo "========================================"
