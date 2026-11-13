#!/bin/bash
#SBATCH -J 08_imp_sets
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-29%6
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_merge_imputed_marker_sets_%A_%a.out
#SBATCH -e logs/08_merge_imputed_marker_sets_%A_%a.err

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

############################################
# 1. Configuration
############################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

CLUMP_ROOT="${WORK}/08_imputed_GWAS/06_imputed_clumping"

GENO_ROOT="${WORK}/08_imputed_GWAS/02_imputed_inputs"

OUTROOT="${WORK}/08_imputed_GWAS/07_imputed_marker_sets"

EXTRA_COUNTS=(0 1 5 10 20 50 100 200 500 1000)

############################################
# 2. Task mapping: 30 inner folds
############################################

TASK=${SLURM_ARRAY_TASK_ID}

OUTER=$(printf "%02d" $((TASK / 3 + 1)))
INNER=$(printf "%02d" $((TASK % 3 + 1)))

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}"

BFILE="${GENO_ROOT}/outer_${OUTER}/inner_${INNER}/imputed_train_qc"

FLAG="${OUT}/complete.flag"

mkdir -p "${OUT}"

echo "========================================"
echo "Step 8.4: Merge imputed clumping results"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

############################################
# 3. Skip completed task
############################################

if [[ -s "${FLAG}" \
      && -s "${OUT}/merge_summary.tsv" \
      && -s "${OUT}/BB/extra_1000.snps" \
      && -s "${OUT}/H6/extra_1000.snps" \
      && -f "${OUT}/BB/extra_0.snps" \
      && -f "${OUT}/H6/extra_0.snps" ]]; then

    echo "SKIP: Already completed"
    exit 0
fi

rm -f "${FLAG}"

############################################
# 4. Input validation
############################################

if [[ ! -s "${BFILE}.bim" || ! -s "${BFILE}.fam" ]]; then
    echo "ERROR: Missing fold-specific genotype files" >&2
    exit 1
fi

N_SAMPLES=$(wc -l < "${BFILE}.fam")

if (( N_SAMPLES < 500 || N_SAMPLES > 600 )); then
    echo "ERROR: Unexpected sample count: ${N_SAMPLES}" >&2
    exit 1
fi

############################################
# 5. Staging directory
############################################

TMP=$(mktemp -d "${OUT}/.tmp.XXXXXX")

trap 'rm -rf "${TMP}"' EXIT

SUMMARY="${TMP}/merge_summary.tsv"

printf "Outer\tInner\tTrait\tSamples\tChromosomes\tTotal_Index_SNPs\tMax_Extra\n" \
    > "${SUMMARY}"

############################################
# 6. Process BB and H6
############################################

for TRAIT in BB H6
do
    echo "========================================"
    echo "Processing trait: ${TRAIT}"
    echo "========================================"

    TRAIT_TMP="${TMP}/${TRAIT}"
    mkdir -p "${TRAIT_TMP}"

    UNSORTED="${TRAIT_TMP}/unsorted.tsv"
    RANKED="${TRAIT_TMP}/ranked_imputed_snps.tsv"

    printf "SNP\tP\tChr\n" > "${UNSORTED}"

    ########################################
    # 7. Merge 12 chromosome index lists
    ########################################

    for CHR in $(seq 1 12)
    do
        CHR_PAD=$(printf "%02d" "${CHR}")

        CHR_DIR="${CLUMP_ROOT}/outer_${OUTER}/inner_${INNER}/chr_${CHR_PAD}"

        SRC="${CHR_DIR}/${TRAIT}_chr${CHR_PAD}.index_ranked.tsv"

        CHR_FLAG="${CHR_DIR}/complete.flag"

        if [[ ! -s "${CHR_FLAG}" ]]; then
            echo "ERROR: Incomplete clumping ${CHR_DIR}" >&2
            exit 1
        fi

        if [[ ! -s "${SRC}" ]]; then
            echo "ERROR: Missing ${SRC}" >&2
            exit 1
        fi

        awk -v chr="${CHR}" '
        BEGIN {
            OFS = "\t"
        }

        FNR == 1 {
            if (NF != 2 || $1 != "SNP" || $2 != "P") {
                print "ERROR: Invalid clumping header" > "/dev/stderr"
                bad = 1
                exit 1
            }
            next
        }

        {
            if (NF != 2) {
                print "ERROR: Invalid index SNP row" > "/dev/stderr"
                bad = 1
                exit 1
            }

            snp = $1
            p = $2

            split(snp, fields, ":")

            if (fields[1]+0 != chr) {
                print "ERROR: SNP chromosome mismatch:", snp > "/dev/stderr"
                bad = 1
                exit 1
            }

            if (p !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
                print "ERROR: Invalid P-value:", p > "/dev/stderr"
                bad = 1
                exit 1
            }

            if (p+0 <= 0 || p+0 > 1) {
                print "ERROR: P-value outside (0,1]" > "/dev/stderr"
                bad = 1
                exit 1
            }

            print snp, p, chr
            n++
        }

        END {
            if (bad || n < 1) {
                print "ERROR: No valid index SNPs" > "/dev/stderr"
                exit 1
            }
        }
        ' "${SRC}" >> "${UNSORTED}"

        echo "PASS: ${TRAIT} chr${CHR_PAD}"
    done

    ########################################
    # 8. Check duplicate SNP IDs
    ########################################

    awk '
    NR == 1 { next }

    {
        if (seen[$1]++) {
            print "ERROR: Duplicate index SNP:", $1 > "/dev/stderr"
            bad = 1
        }
    }

    END {
        if (bad) exit 1
    }
    ' "${UNSORTED}"

    ########################################
    # 9. Global ranking by GWAS P-value
    ########################################

    {
        printf "SNP\tP\tChr\n"

        tail -n +2 "${UNSORTED}" \
            | LC_ALL=C sort \
                -t $'\t' \
                -k2,2g \
                -k1,1

    } > "${RANKED}"

    N_TOTAL=$(awk 'END { print NR-1 }' "${RANKED}")

    echo "Total clumped index SNPs: ${N_TOTAL}"

    if (( N_TOTAL < 1000 )); then
        echo "ERROR: Fewer than 1000 candidate SNPs" >&2
        exit 1
    fi

    ########################################
    # 10. Verify all marker IDs and
    #     chromosomes against PLINK BIM
    ########################################

    awk '
    NR == FNR {
        if (FNR == 1) next

        id = $1

        if (id in wanted) {
            print "ERROR: Duplicate ranked SNP:", id > "/dev/stderr"
            bad = 1
        }

        wanted[id] = $3
        total++
        next
    }

    ($2 in wanted) {
        if ($1+0 != wanted[$2]+0) {
            print "ERROR: SNP chromosome mismatch in BIM:", $2 > "/dev/stderr"
            bad = 1
        }

        found[$2]++
        matched++
    }

    END {
        if (matched != total) {
            print "ERROR: Ranked SNPs not fully represented in BIM:", matched, total > "/dev/stderr"
            bad = 1
        }

        if (bad) exit 1
    }
    ' "${RANKED}" "${BFILE}.bim"

    echo "PASS: All candidate SNPs found in training BIM"

    ########################################
    # 11. Verify global P-value order
    ########################################

    awk '
    NR == 1 { next }

    {
        p = $2+0

        if (NR > 2 && p < previous) {
            print "ERROR: Global P-value ranking is not ascending" > "/dev/stderr"
            bad = 1
        }

        previous = p
    }

    END {
        if (bad) exit 1
    }
    ' "${RANKED}"

    ########################################
    # 12. Generate nested extra SNP sets
    ########################################

    for N in "${EXTRA_COUNTS[@]}"
    do
        SET="${TRAIT_TMP}/extra_${N}.snps"

        if (( N == 0 )); then
            : > "${SET}"
        else
            awk -v limit="${N}" '
            NR > 1 && NR <= limit+1 {
                print $1
            }
            ' "${RANKED}" > "${SET}"
        fi

        ACTUAL=$(wc -l < "${SET}")

        if (( ACTUAL != N )); then
            echo "ERROR: extra_${N} SNP count mismatch" >&2
            exit 1
        fi

        echo "${TRAIT} extra_${N}: ${ACTUAL} SNPs"
    done

    ########################################
    # 13. Publish validated trait files
    ########################################

    TRAIT_OUT="${OUT}/${TRAIT}"

    mkdir -p "${TRAIT_OUT}"

    mv -f "${RANKED}" \
        "${TRAIT_OUT}/ranked_imputed_snps.tsv"

    for N in "${EXTRA_COUNTS[@]}"
    do
        mv -f "${TRAIT_TMP}/extra_${N}.snps" \
            "${TRAIT_OUT}/extra_${N}.snps"
    done

    printf "%s\t%s\t%s\t%s\t12\t%s\t1000\n" \
        "${OUTER}" \
        "${INNER}" \
        "${TRAIT}" \
        "${N_SAMPLES}" \
        "${N_TOTAL}" >> "${SUMMARY}"

    echo "PASS: ${TRAIT} marker sets completed"

done

############################################
# 14. Save summary
############################################

mv -f "${SUMMARY}" "${OUT}/merge_summary.tsv"

############################################
# 15. Final output checks
############################################

for TRAIT in BB H6
do
    [[ -s "${OUT}/${TRAIT}/ranked_imputed_snps.tsv" ]] || exit 1

    for N in "${EXTRA_COUNTS[@]}"
    do
        F="${OUT}/${TRAIT}/extra_${N}.snps"

        [[ -f "${F}" ]] || exit 1

        ACTUAL=$(wc -l < "${F}")

        if (( ACTUAL != N )); then
            echo "ERROR: Final marker set count mismatch" >&2
            exit 1
        fi
    done
done

############################################
# 16. Completion flag
############################################

printf "DONE Step8.4 outer_%s inner_%s\n" \
    "${OUTER}" "${INNER}" \
    > "${FLAG}"

echo "========================================"
echo "DONE Step8.4"
echo "========================================"

