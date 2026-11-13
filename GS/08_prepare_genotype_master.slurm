#!/bin/bash
#SBATCH -J 08_genotype_master
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_prepare_genotype_master_%j.out
#SBATCH -e logs/08_prepare_genotype_master_%j.err

set -euo pipefail

########################################
# 1. Modules
########################################

module load bioinfo-tools
ml plink/1.90b4.9

########################################
# 2. Paths
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

VCF="${WORK}/01_reference/Imputed_array_combine.ID.vcf.gz"

CHIP_LIST="${WORK}/01_reference/chip.marker"
IMP_LIST="${WORK}/01_reference/imputed.marker"

SAMPLES="${WORK}/02_CV/all.samples"

OUTDIR="${WORK}/08_imputed_GWAS/00_reference"

PREFIX="${OUTDIR}/all_904_genotypes"

mkdir -p "${OUTDIR}"

########################################
# 3. Input checks
########################################

echo "=========================================="
echo "Step 8.1: Prepare master genotype"
echo "=========================================="

for F in "${VCF}" "${CHIP_LIST}" \
         "${IMP_LIST}" "${SAMPLES}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: Missing or empty input: ${F}" >&2
        exit 1
    fi
done

N_CHIP=$(wc -l < "${CHIP_LIST}")
N_IMP=$(wc -l < "${IMP_LIST}")
N_SAMPLE=$(wc -l < "${SAMPLES}")

echo "Chip marker IDs: ${N_CHIP}"
echo "Imputed marker IDs: ${N_IMP}"
echo "Samples: ${N_SAMPLE}"

if [[ "${N_CHIP}" -ne 43267 ]]; then
    echo "ERROR: Unexpected chip marker count" >&2
    exit 1
fi

if [[ "${N_IMP}" -ne 3121909 ]]; then
    echo "ERROR: Unexpected imputed marker count" >&2
    exit 1
fi

if [[ "${N_SAMPLE}" -ne 904 ]]; then
    echo "ERROR: Expected 904 samples" >&2
    exit 1
fi

########################################
# 4. Check marker list independence
########################################

echo "Checking marker lists"

awk '
NR==FNR {
    if (NF != 1 || seen[$1]++) {
        print "ERROR: Invalid chip marker list" > "/dev/stderr"
        exit 1
    }
    chip[$1] = 1
    next
}
{
    if (NF != 1) {
        print "ERROR: Invalid imputed marker list" > "/dev/stderr"
        exit 1
    }
    if ($1 in chip) {
        print "ERROR: Chip/imputed overlap:", $1 > "/dev/stderr"
        bad = 1
    }
}
END {
    if (bad) exit 1
}
' "${CHIP_LIST}" "${IMP_LIST}"

########################################
# 5. Convert VCF to PLINK binary
########################################

# Important:
# No phenotype is used.
# No genome-wide MAF, missingness or LD
# filtering is applied at this stage.
#
# These operations will be conducted
# separately within each training fold.

echo "Converting full VCF to PLINK"

rm -f "${PREFIX}.complete"

plink \
    --vcf "${VCF}" \
    --double-id \
    --allow-extra-chr \
    --keep-allele-order \
    --make-bed \
    --threads 4 \
    --memory 48000 \
    --out "${PREFIX}"

########################################
# 6. Check output
########################################

for EXT in bed bim fam
do
    if [[ ! -s "${PREFIX}.${EXT}" ]]; then
        echo "ERROR: Missing output ${EXT}" >&2
        exit 1
    fi
done

echo "Checking sample identities"

awk '
NR==FNR {
    if (NF != 1 || ($1 in expected)) {
        print "ERROR: Invalid sample list" > "/dev/stderr"
        exit 1
    }
    expected[$1] = 1
    total++
    next
}
{
    observed++

    if (NF < 6 || $1 != $2) {
        bad = 1
    }

    if (!($2 in expected) || seen[$2]++) {
        bad = 1
    }
}
END {
    if (total != 904 || observed != total) {
        bad = 1
    }

    for (s in expected) {
        if (!(s in seen)) {
            bad = 1
        }
    }

    if (bad) {
        print "ERROR: Master sample IDs do not match" > "/dev/stderr"
        exit 1
    }

    print "PASS: 904 matching samples"
}
' "${SAMPLES}" "${PREFIX}.fam"

########################################
# 7. Verify chip marker retention
########################################

echo "Checking chip marker retention"

awk '
NR==FNR {
    wanted[$1] = 1
    next
}
{
    delete wanted[$2]
}
END {
    for (s in wanted) print s
}
' "${CHIP_LIST}" "${PREFIX}.bim" \
    > "${OUTDIR}/missing_chip_markers.txt"

N_MISSING_CHIP=$(
    wc -l < "${OUTDIR}/missing_chip_markers.txt"
)

echo "Missing chip markers: ${N_MISSING_CHIP}"

if [[ "${N_MISSING_CHIP}" -ne 0 ]]; then
    echo "ERROR: Some chip markers were lost" >&2
    head "${OUTDIR}/missing_chip_markers.txt"
    exit 1
fi

########################################
# 8. Verify imputed marker retention
########################################

echo "Checking imputed marker retention"

awk '
NR==FNR {
    wanted[$1] = 1
    next
}
{
    delete wanted[$2]
}
END {
    for (s in wanted) print s
}
' "${IMP_LIST}" "${PREFIX}.bim" \
    > "${OUTDIR}/missing_imputed_markers.txt"

N_MISSING_IMP=$(
    wc -l < "${OUTDIR}/missing_imputed_markers.txt"
)

echo "Missing imputed markers: ${N_MISSING_IMP}"

# A small number of variants may be
# unsupported by PLINK VCF import.
# Do not silently accept a large loss.

if [[ "${N_MISSING_IMP}" -gt 1000 ]]; then
    echo "ERROR: Excessive imputed marker loss" >&2
    exit 1
fi

########################################
# 9. Report chromosome counts
########################################

echo "Chromosome distribution"

awk '
{
    count[$1]++
}
END {
    for (chr in count) {
        print chr, count[chr]
    }
}
' "${PREFIX}.bim" \
    | sort -V \
    | tee "${OUTDIR}/chromosome_counts.txt"

########################################
# 10. Final report
########################################

N_VARIANTS=$(wc -l < "${PREFIX}.bim")
N_PEOPLE=$(wc -l < "${PREFIX}.fam")

echo "=========================================="
echo "Master genotype summary"
echo "=========================================="
echo "Samples: ${N_PEOPLE}"
echo "Total variants: ${N_VARIANTS}"
echo "Missing chip SNP: ${N_MISSING_CHIP}"
echo "Missing imputed SNP: ${N_MISSING_IMP}"

touch "${PREFIX}.complete"

echo "DONE Step8.1"
