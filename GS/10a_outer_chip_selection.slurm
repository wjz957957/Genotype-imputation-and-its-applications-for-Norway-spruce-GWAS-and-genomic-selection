#!/bin/bash
#SBATCH -J 10a_outer_chip
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-9%3
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/10a_outer_chip_%A_%a.out
#SBATCH -e logs/10a_outer_chip_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

########################################
# 1. Software
########################################

module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1

########################################
# 2. Paths
########################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"

CHIP_MARKERS="${WORK}/01_reference/chip.marker"

PHENO_ROOT="${WORK}/00_input"

OPT_ROOT="${WORK}/07_chip_RRBLUP/optimization"

OUTER=$(printf '%02d' "$(( SLURM_ARRAY_TASK_ID + 1 ))")

CV="${WORK}/02_CV/outer_${OUTER}"

TRAIN="${CV}/outer_train.samples"
TEST="${CV}/outer_test.samples"

OUT="${WORK}/10_outer_test/01_chip_selection/outer_${OUTER}"

CHIP="${OUT}/plink/chip_train"
GRM="${OUT}/grm/chip_grm"

FLAG="${OUT}/complete.flag"

SUMMARY="${OUT}/chip_selection_summary.tsv"

QC="${OUT}/gwas_pvalue_qc.tsv"

mkdir -p \
    "${OUT}/plink" \
    "${OUT}/grm" \
    "${OUT}/pheno" \
    "${OUT}/gwas" \
    "${OUT}/clumping" \
    "${OUT}/selected"

echo "========================================"
echo "Step 10a: Outer-train chip selection"
echo "Outer: ${OUTER}"
echo "========================================"

########################################
# 3. Completion check
########################################

if [[ -s "${FLAG}" \
      && -s "${SUMMARY}" \
      && -s "${QC}" \
      && -s "${OUT}/selected/BB_best_chip.snps" \
      && -s "${OUT}/selected/H6_best_chip.snps" ]]; then

    echo "SKIP: Previously completed"
    exit 0
fi

rm -f "${FLAG}"

########################################
# 4. Check inputs
########################################

for EXT in bed bim fam
do
    [[ -s "${MASTER}.${EXT}" ]] || {
        echo "ERROR: Missing master ${EXT}" >&2
        exit 1
    }
done

for F in "${CHIP_MARKERS}" "${TRAIN}" "${TEST}"
do
    [[ -s "${F}" ]] || {
        echo "ERROR: Missing input ${F}" >&2
        exit 1
    }
done

if [[ $(wc -l < "${CHIP_MARKERS}") -ne 43267 ]]; then
    echo "ERROR: Expected 43267 chip markers" >&2
    exit 1
fi

########################################
# 5. Training samples
########################################

KEEP="${OUT}/plink/outer_train.keep"

awk '
{
    sub(/\r$/, "", $0)

    if (NF != 1 || seen[$1]++) {
        print "ERROR: Invalid or duplicate training ID" > "/dev/stderr"
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

if (( N_TRAIN != 813 && N_TRAIN != 814 )); then
    echo "ERROR: Unexpected outer training count" >&2
    exit 1
fi

if (( N_TRAIN + N_TEST != 904 )); then
    echo "ERROR: Outer train/test total mismatch" >&2
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

    if ($1 in train) {
        print "ERROR: Train/test overlap:", $1 > "/dev/stderr"
        bad=1
    }
}
END {
    if (bad) exit 1
}
' "${TRAIN}" "${TEST}"

echo "Outer training samples: ${N_TRAIN}"
echo "Outer testing samples: ${N_TEST}"

########################################
# 6. Extract chip genotypes
########################################

plink \
    --bfile "${MASTER}" \
    --keep "${KEEP}" \
    --extract "${CHIP_MARKERS}" \
    --make-bed \
    --threads 4 \
    --memory 28000 \
    --out "${CHIP}"

for EXT in bed bim fam
do
    [[ -s "${CHIP}.${EXT}" ]] || {
        echo "ERROR: Missing chip ${EXT}" >&2
        exit 1
    }
done

N_CHIP=$(wc -l < "${CHIP}.bim")
N_FAM=$(wc -l < "${CHIP}.fam")

if (( N_CHIP != 43267 || N_FAM != N_TRAIN )); then
    echo "ERROR: Chip genotype dimension mismatch" >&2
    exit 1
fi

echo "PASS: ${N_CHIP} chip SNPs"
echo "PASS: ${N_FAM} outer-training individuals"

########################################
# 7. Prepare phenotypes
########################################

for TRAIT in BB H6
do

    if [[ "${TRAIT}" == "BB" ]]; then
        SOURCE="${PHENO_ROOT}/Trait-Adj_Budburst.txt"
    else
        SOURCE="${PHENO_ROOT}/Trait-Adj_H6.txt"
    fi

    [[ -s "${SOURCE}" ]] || {
        echo "ERROR: Missing phenotype ${SOURCE}" >&2
        exit 1
    }

    TARGET="${OUT}/pheno/${TRAIT}.pheno"

    awk -v expected="${N_TRAIN}" '
    BEGIN {
        OFS="\t"
    }

    NR == FNR {
        sub(/\r$/, "", $0)

        if (FNR == 1) next

        if (NF != 3 || $1 != $2 || seen[$2]++) {
            print "ERROR: Invalid phenotype row" > "/dev/stderr"
            bad=1
            exit 1
        }

        value=$3

        gsub(/^[[:space:]]+/, "", value)
        gsub(/[[:space:]]+$/, "", value)

        pheno[$2]=value
        next
    }

    {
        sub(/\r$/, "", $0)

        id=$1

        if (!(id in pheno)) {
            print "ERROR: Missing phenotype:", id > "/dev/stderr"
            bad=1
            exit 1
        }

        value=pheno[id]

        if (value !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
            print "ERROR: Non-numeric phenotype:", id > "/dev/stderr"
            bad=1
            exit 1
        }

        print id, id, value
        n++
    }

    END {
        if (n != expected || bad) {
            print "ERROR: Phenotype count mismatch" > "/dev/stderr"
            exit 1
        }
    }
    ' "${SOURCE}" "${TRAIN}" > "${TARGET}"

    echo "PASS: ${TRAIT} phenotype $(wc -l < "${TARGET}")"
done

########################################
# 8. Build or reuse GRM
########################################

if [[ -s "${GRM}.grm.bin" \
      && -s "${GRM}.grm.N.bin" \
      && -s "${GRM}.grm.id" ]]; then

    echo "REUSE: Existing outer-train GRM"

else

    gcta64 \
        --bfile "${CHIP}" \
        --make-grm \
        --thread-num 4 \
        --out "${GRM}"

fi

for EXT in grm.bin grm.N.bin grm.id
do
    [[ -s "${GRM}.${EXT}" ]] || {
        echo "ERROR: Missing GRM ${EXT}" >&2
        exit 1
    }
done

if [[ $(wc -l < "${GRM}.grm.id") -ne ${N_TRAIN} ]]; then
    echo "ERROR: GRM sample count mismatch" >&2
    exit 1
fi

########################################
# 9. Initialize summaries
########################################

printf 'Outer\tTrait\tBest_Marker\tN_train\tN_all_chip\tN_clumped_index\tN_selected_chip\n' \
    > "${SUMMARY}"

printf 'Outer\tTrait\tGWAS_Total\tValid_P\tNonfinite_P\tZero_P_Floored\n' \
    > "${QC}"

########################################
# 10. Trait-specific GWAS and clumping
########################################

for TRAIT in BB H6
do

    echo "========================================"
    echo "Trait: ${TRAIT}"
    echo "========================================"

    PREFIX="${OUT}/gwas/chip_${TRAIT}_gwas"

    GWAS="${PREFIX}.loco.mlma"

    ####################################
    # 10.1 Run or reuse GWAS
    ####################################

    if [[ -s "${GWAS}" ]]; then

        echo "REUSE: Existing ${TRAIT} GWAS"
        echo "Validating all GWAS records"

    else

        gcta64 \
            --mlma-loco \
            --bfile "${CHIP}" \
            --grm "${GRM}" \
            --pheno "${OUT}/pheno/${TRAIT}.pheno" \
            --thread-num 4 \
            --out "${PREFIX}"

    fi

    [[ -s "${GWAS}" ]] || {
        echo "ERROR: Missing GWAS ${GWAS}" >&2
        exit 1
    }

    ####################################
    # 10.2 Prepare clumping directory
    ####################################

    CLUMP="${OUT}/clumping/${TRAIT}"

    mkdir -p \
        "${CLUMP}/clump" \
        "${CLUMP}/marker_sets"

    PFILE="${CLUMP}/gwas_p.txt"

    QC_TMP="${CLUMP}/pvalue_counts.txt"

    ####################################
    # 10.3 Validate GWAS and filter
    #      non-estimable P-values
    ####################################

    awk -v expected="${N_CHIP}" -v counts="${QC_TMP}" '
    BEGIN {
        OFS="\t"
        print "SNP", "P"
    }

    NR == 1 {
        if (tolower($1) != "chr" || tolower($2) != "snp" || tolower($9) != "p") {
            print "ERROR: Unexpected GWAS header" > "/dev/stderr"
            bad=1
            exit 1
        }
        next
    }

    {
        if (NF < 9 || seen[$2]++) {
            print "ERROR: Malformed or duplicate GWAS SNP:", $2 > "/dev/stderr"
            bad=1
            exit 1
        }

        total++

        p=$9
        low=tolower(p)

        # GCTA may output -nan for SNPs
        # whose association test is not estimable.

        if (low ~ /^[-+]?nan(\([^)]*\))?$/ || low ~ /^[-+]?inf(inity)?$/ || low == "na" || low == "n/a" || low == "." || low == "-na") {
            invalid++
            next
        }

        # All other values must be
        # valid numeric probabilities.

        if (p !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || p+0 < 0 || p+0 > 1) {
            print "ERROR: Unexpected GWAS P-value for", $2, "[" p "]" > "/dev/stderr"
            bad=1
            exit 1
        }

        # Avoid zero P-values in PLINK clumping.
        # This affects only numerical underflow.

        if (p+0 == 0) {
            p="1e-300"
            zero++
        }

        print $2, p
        valid++
    }

    END {
        if (bad || total != expected || valid < 1 || valid+invalid != total) {
            print "ERROR: GWAS QC failed. Total:", total+0, "Valid:", valid+0, "Nonfinite:", invalid+0, "Expected:", expected > "/dev/stderr"
            exit 1
        }

        printf "%d\t%d\t%d\t%d\n", total, valid, invalid, zero+0 > counts
        close(counts)
    }
    ' "${GWAS}" > "${PFILE}"

    ####################################
    # 10.4 Read GWAS QC counts
    ####################################

    read -r TOTAL VALID INVALID ZERO < "${QC_TMP}"

    if [[ $(wc -l < "${PFILE}") -ne $((VALID + 1)) ]]; then
        echo "ERROR: Clumping P-value file count mismatch" >&2
        exit 1
    fi

    echo "========================================"
    echo "GWAS P-value QC: ${TRAIT}"
    echo "Total SNPs: ${TOTAL}"
    echo "Valid P-values: ${VALID}"
    echo "Nonfinite P-values excluded: ${INVALID}"
    echo "Zero P-values floored: ${ZERO}"
    echo "========================================"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${OUTER}" "${TRAIT}" \
        "${TOTAL}" "${VALID}" "${INVALID}" "${ZERO}" \
        >> "${QC}"

    ####################################
    # 10.5 PLINK chip LD clumping
    #
    # Same parameters as original Step 6.
    ####################################

    plink \
        --bfile "${CHIP}" \
        --clump "${PFILE}" \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.2 \
        --clump-kb 50 \
        --threads 4 \
        --memory 28000 \
        --out "${CLUMP}/clump/chip"

    CLUMPED="${CLUMP}/clump/chip.clumped"

    [[ -s "${CLUMPED}" ]] || {
        echo "ERROR: Missing PLINK .clumped" >&2
        exit 1
    }

    ####################################
    # 10.6 Ranked chip SNPs
    ####################################

    RANKED="${CLUMP}/ranked_chip_snps.txt"

    awk '
    NR == 1 {
        if ($1 != "CHR" || $3 != "SNP" || $5 != "P") {
            print "ERROR: Unexpected clumping header" > "/dev/stderr"
            bad=1
            exit 1
        }
        next
    }

    NF == 0 { next }

    {
        if ($3 == "NONE" || $3 == "NA" || seen[$3]++) {
            print "ERROR: Invalid clump index SNP" > "/dev/stderr"
            bad=1
            exit 1
        }

        print $3 "\t" $5
        n++
    }

    END {
        if (bad || n < 1) exit 1
    }
    ' "${CLUMPED}" \
        | LC_ALL=C sort -k2,2g -k1,1 \
        > "${RANKED}"

    N_RANKED=$(wc -l < "${RANKED}")

    echo "Valid clumped index SNPs: ${N_RANKED}"

    ####################################
    # 10.7 Reconstruct original
    #      Step 6 marker sets
    ####################################

    for N in \
        50 100 200 300 400 500 \
        1000 2000 3000 4000 5000 \
        8000 10000 20000 30000 40000
    do

        head -n "${N}" "${RANKED}" \
            | awk '{print $1}' \
            > "${CLUMP}/marker_sets/chip_${N}.snps"

    done

    awk '{print $1}' "${RANKED}" \
        > "${CLUMP}/marker_sets/chip_LD_ALL.snps"

    # chip_ALL includes all chip SNPs,
    # regardless of GWAS P-value availability.

    awk '{print $2}' "${CHIP}.bim" \
        > "${CLUMP}/marker_sets/chip_ALL.snps"

    ####################################
    # 10.8 Load previously optimized chip
    ####################################

    BEST_FILE="${OPT_ROOT}/outer_${OUTER}/${TRAIT}_best_marker.txt"

    [[ -s "${BEST_FILE}" ]] || {
        echo "ERROR: Missing best chip setting" >&2
        exit 1
    }

    BEST=$(tr -d '[:space:]' < "${BEST_FILE}")

    case "${BEST}" in
        chip_ALL|chip_LD_ALL|chip_50|chip_100|chip_200|chip_300|chip_400|chip_500|chip_1000|chip_2000|chip_3000|chip_4000|chip_5000|chip_8000|chip_10000|chip_20000|chip_30000|chip_40000)
            ;;
        *)
            echo "ERROR: Unknown best chip setting ${BEST}" >&2
            exit 1
            ;;
    esac

    SELECTED="${CLUMP}/marker_sets/${BEST}.snps"

    [[ -s "${SELECTED}" ]] || {
        echo "ERROR: Missing selected chip SNP list" >&2
        exit 1
    }

    N_SELECTED=$(wc -l < "${SELECTED}")

    if [[ "${BEST}" == "chip_ALL" ]]; then

        EXPECTED_SELECTED=${N_CHIP}

    elif [[ "${BEST}" == "chip_LD_ALL" ]]; then

        EXPECTED_SELECTED=${N_RANKED}

    else

        EXPECTED_SELECTED=${BEST#chip_}

    fi

    if (( N_SELECTED != EXPECTED_SELECTED )); then
        echo "ERROR: Selected chip SNP count mismatch" >&2
        echo "Best: ${BEST}" >&2
        echo "Expected: ${EXPECTED_SELECTED}" >&2
        echo "Actual: ${N_SELECTED}" >&2
        exit 1
    fi

    ####################################
    # 10.9 Verify selected IDs in BIM
    ####################################

    awk '
    NR == FNR {
        if (wanted[$1]++) {
            print "ERROR: Duplicate selected SNP:", $1 > "/dev/stderr"
            bad=1
        }
        next
    }

    ($2 in wanted) {
        found[$2]++
    }

    END {
        for (snp in wanted) {
            if (found[snp] != 1) {
                print "ERROR: Selected SNP absent from BIM:", snp > "/dev/stderr"
                bad=1
            }
        }

        if (bad) exit 1
    }
    ' "${SELECTED}" "${CHIP}.bim"

    ####################################
    # 10.10 Publish final chip SNP set
    ####################################

    cp "${SELECTED}" \
        "${OUT}/selected/${TRAIT}_best_chip.snps"

    printf '%s\n' "${BEST}" \
        > "${OUT}/selected/${TRAIT}_best_marker.txt"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${OUTER}" \
        "${TRAIT}" \
        "${BEST}" \
        "${N_TRAIN}" \
        "${N_CHIP}" \
        "${N_RANKED}" \
        "${N_SELECTED}" \
        >> "${SUMMARY}"

    echo "PASS: ${TRAIT}"
    echo "Best chip configuration: ${BEST}"
    echo "Selected chip SNPs: ${N_SELECTED}"

done

########################################
# 11. Final checks
########################################

[[ $(wc -l < "${SUMMARY}") -eq 3 ]] || {
    echo "ERROR: Incorrect selection summary" >&2
    exit 1
}

[[ $(wc -l < "${QC}") -eq 3 ]] || {
    echo "ERROR: Incorrect GWAS QC summary" >&2
    exit 1
}

########################################
# 12. Completion flag
########################################

printf "DONE Step10a outer_%s\n" "${OUTER}" \
    > "${FLAG}"

echo "========================================"
echo "DONE Step10a outer_${OUTER}"
echo "========================================"

