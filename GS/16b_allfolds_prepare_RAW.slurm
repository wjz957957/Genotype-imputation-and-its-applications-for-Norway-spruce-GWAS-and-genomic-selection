#!/bin/bash
#SBATCH -J 16b_crossRAW
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=0-59%5
#SBATCH -o logs/16b_crossRAW_%A_%a.out
#SBATCH -e logs/16b_crossRAW_%A_%a.err

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
SET_ROOT="${WORK}/08_imputed_GWAS/09_crosschip_LD_r2_0p8_allfolds"
IMP_ROOT="${WORK}/08_imputed_GWAS/02_imputed_inputs"
OUTROOT="${WORK}/09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_allfolds"

TASK=${SLURM_ARRAY_TASK_ID}

(( TASK >= 0 && TASK < 60 )) || { echo "ERROR task $TASK" >&2; exit 1; }
OUTER=$(printf '%02d' "$(( TASK / 6 + 1 ))")
INNER=$(printf '%02d' "$(( TASK % 3 + 1 ))")
if (( TASK % 6 < 3 )); then TRAIT=BB; else TRAIT=H6; fi

CV="${WORK}/02_CV/outer_${OUTER}/inner_${INNER}"
SETS="${SET_ROOT}/outer_${OUTER}/inner_${INNER}/${TRAIT}/r2_0p8"
IMP_BIM="${IMP_ROOT}/outer_${OUTER}/inner_${INNER}/imputed_train_qc.bim"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}/${TRAIT}"
GENO="${OUT}/genotype"
FLAG="${OUT}/step9a.complete.flag"

mkdir -p "${GENO}"

echo "========================================"
echo "Step 16b: Prepare filtered imputed genotypes ${TRAIT} (r2=0.8)"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

if [[ -s "${FLAG}" ]]; then
    READY=1

    for TRAIT in "$TRAIT"
    do
        for TYPE in train valid
        do
            [[ -s "${GENO}/imputed_${TRAIT}_top700_${TYPE}.raw" ]] \
                || READY=0
        done
    done

    if (( READY == 1 )); then
        echo "SKIP: Already completed"
        exit 0
    fi
fi

rm -f "${FLAG}"

########################################
# Check inputs
########################################

for EXT in bed bim fam
do
    [[ -s "${MASTER}.${EXT}" ]] || {
        echo "ERROR: Missing master ${EXT}" >&2
        exit 1
    }
done

[[ -s "${SETS}/complete.flag" ]] || {
    echo "ERROR: 16a cross-CHIP LD filtering incomplete" >&2
    exit 1
}

[[ -s "${IMP_BIM}" ]] || {
    echo "ERROR: Missing imputed QC BIM" >&2
    exit 1
}

########################################
# Prepare sample keep files
########################################

for TYPE in train valid
do
    SRC="${CV}/inner_${TYPE}.samples"
    DEST="${GENO}/${TYPE}.keep"

    [[ -s "${SRC}" ]] || {
        echo "ERROR: Missing sample file ${SRC}" >&2
        exit 1
    }

    awk '
    NF != 1 {
        print "ERROR: Invalid sample row" > "/dev/stderr"
        bad=1
        exit 1
    }
    {
        if (seen[$1]++) {
            print "ERROR: Duplicate sample ID" > "/dev/stderr"
            bad=1
            exit 1
        }
        print $1, $1
    }
    END {
        if (bad || NR < 1) exit 1
    }
    ' "${SRC}" > "${DEST}"
done

N_TRAIN=$(wc -l < "${GENO}/train.keep")
N_VALID=$(wc -l < "${GENO}/valid.keep")

if (( N_TRAIN < 500 || N_TRAIN > 600 )); then
    echo "ERROR: Invalid training sample count" >&2
    exit 1
fi

if (( N_VALID < 250 || N_VALID > 300 )); then
    echo "ERROR: Invalid validation sample count" >&2
    exit 1
fi

awk '
NR == FNR {
    ids[$1]=1
    next
}
($1 in ids) {
    print "ERROR: Train/valid sample overlap" > "/dev/stderr"
    bad=1
}
END {
    if (bad) exit 1
}
' "${GENO}/train.keep" "${GENO}/valid.keep"

echo "Training samples: ${N_TRAIN}"
echo "Validation samples: ${N_VALID}"

########################################
# Process traits separately
########################################

for TRAIT in "$TRAIT"
do
    echo "========================================"
    echo "Trait: ${TRAIT}"
    echo "========================================"

    SNPS="${SETS}/marker_sets/extra_700.snps"

    [[ -s "${SNPS}" ]] || {
        echo "ERROR: Missing marker list ${SNPS}" >&2
        exit 1
    }

    [[ $(wc -l < "${SNPS}") -eq 700 ]] || {
        echo "ERROR: Expected exactly 700 SNPs" >&2
        exit 1
    }

    # The filtered, ranked top700 list must agree with the published marker sets.
    cmp -s "${SNPS}" "${SETS}/top700_extra.snps" || {
        echo "ERROR: Pilot top700 markers mismatch" >&2
        exit 1
    }

    ####################################
    # Verify marker IDs belong to
    # fold-specific imputed QC variants
    ####################################

    awk '
    NR == FNR {
        if (NF != 1 || wanted[$1]++) {
            print "ERROR: Invalid candidate SNP list" > "/dev/stderr"
            bad=1
        }
        next
    }

    ($2 in wanted) {
        found[$2]++
    }

    END {
        if (bad) exit 1

        n=0

        for (snp in wanted) {
            if (found[snp] != 1) {
                print "ERROR: Candidate SNP absent from imputed QC:", snp > "/dev/stderr"
                bad=1
            }
            n++
        }

        if (n != 700 || bad) exit 1
    }
    ' "${SNPS}" "${IMP_BIM}"

    echo "PASS: All 700 SNPs belong to inner-train imputed QC"

    ####################################
    # Extract train and valid genotype
    ####################################

    for TYPE in train valid
    do
        PREFIX="${GENO}/imputed_${TRAIT}_top700_${TYPE}"

        echo "Extracting: ${TRAIT} ${TYPE}"

        plink \
            --bfile "${MASTER}" \
            --keep "${GENO}/${TYPE}.keep" \
            --extract "${SNPS}" \
            --recode A \
            --threads 2 \
            --memory 10000 \
            --out "${PREFIX}"

        RAW="${PREFIX}.raw"

        [[ -s "${RAW}" ]] || {
            echo "ERROR: Missing PLINK RAW ${RAW}" >&2
            exit 1
        }

        if [[ "${TYPE}" == "train" ]]; then
            EXPECTED="${N_TRAIN}"
        else
            EXPECTED="${N_VALID}"
        fi

        awk -v expected="${EXPECTED}" '
        NR == 1 {
            if (NF != 706) {
                print "ERROR: Expected 700 SNP columns, found:", NF-6 > "/dev/stderr"
                bad=1
            }
            next
        }

        END {
            if (NR-1 != expected) {
                print "ERROR: Incorrect sample count:", NR-1, expected > "/dev/stderr"
                bad=1
            }
            if (bad) exit 1
        }
        ' "${RAW}"

        echo "PASS: ${TRAIT} ${TYPE} ${EXPECTED} samples, 700 SNPs"
    done
done

printf "DONE Step16b outer_%s inner_%s\n" \
    "${OUTER}" "${INNER}" > "${FLAG}"

echo "========================================"
echo "DONE Step16b"
echo "========================================"
