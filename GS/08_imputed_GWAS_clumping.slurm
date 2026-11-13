#!/bin/bash
#SBATCH -J 08_imp_clump
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-359%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_imputed_GWAS_clumping_%A_%a.out
#SBATCH -e logs/08_imputed_GWAS_clumping_%A_%a.err

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

INPUT_ROOT="${WORK}/08_imputed_GWAS/02_imputed_inputs"
GWAS_ROOT="${WORK}/08_imputed_GWAS/03_imputed_GWAS"
OUTROOT="${WORK}/08_imputed_GWAS/06_imputed_clumping"

R2=0.2
WINDOW_KB=1000

TASK=${SLURM_ARRAY_TASK_ID}
FOLD_INDEX=$((TASK / 12))
CHR=$((TASK % 12 + 1))

OUTER=$(printf "%02d" $((FOLD_INDEX / 3 + 1)))
INNER=$(printf "%02d" $((FOLD_INDEX % 3 + 1)))
CHR_PAD=$(printf "%02d" "${CHR}")

BFILE="${INPUT_ROOT}/outer_${OUTER}/inner_${INNER}/imputed_train_qc"

GWAS_DIR="${GWAS_ROOT}/outer_${OUTER}/inner_${INNER}/chr_${CHR_PAD}"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}/chr_${CHR_PAD}"

mkdir -p "${OUT}"

FLAG="${OUT}/complete.flag"
SUMMARY="${OUT}/clumping_summary.tsv"

echo "========================================"
echo "Step 8.3: Imputed GWAS clumping"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "Chromosome: ${CHR}"
echo "R2: ${R2}"
echo "Window: ${WINDOW_KB} kb"
echo "========================================"

if [[ -s "${FLAG}" && -s "${SUMMARY}" \
      && -s "${OUT}/BB_chr${CHR_PAD}.index_ranked.tsv" \
      && -s "${OUT}/H6_chr${CHR_PAD}.index_ranked.tsv" ]]; then
    echo "SKIP: Already complete"
    exit 0
fi

rm -f "${FLAG}"

for EXT in bed bim fam
do
    [[ -s "${BFILE}.${EXT}" ]] || {
        echo "ERROR: Missing ${BFILE}.${EXT}" >&2
        exit 1
    }
done

N_SAMPLES=$(wc -l < "${BFILE}.fam")

if (( N_SAMPLES < 500 || N_SAMPLES > 600 )); then
    echo "ERROR: Unexpected sample count" >&2
    exit 1
fi

N_EXPECT=$(awk -v chr="${CHR}" '
$1 == chr { n++ }
END { print n+0 }
' "${BFILE}.bim")

if (( N_EXPECT < 1000 )); then
    echo "ERROR: Unexpected chromosome SNP count" >&2
    exit 1
fi

printf 'Outer\tInner\tChr\tTrait\tSamples\tGWAS_SNPs\tClumped_Index_SNPs\tR2\tWindow_kb\n' \
    > "${SUMMARY}"

for TRAIT in BB H6
do
    echo "========================================"
    echo "Trait: ${TRAIT}"
    echo "========================================"

    GWAS="${GWAS_DIR}/${TRAIT}_chr${CHR_PAD}.mlma"

    [[ -s "${GWAS}" ]] || {
        echo "ERROR: Missing GWAS ${GWAS}" >&2
        exit 1
    }

    PREFIX="${OUT}/${TRAIT}_chr${CHR_PAD}"
    ASSOC="${PREFIX}.for_clump.tsv"
    RANKED="${PREFIX}.index_ranked.tsv"
    INDEX_LIST="${PREFIX}.index_snps.txt"

    ######################################
    # Prepare GWAS SNP + P-value input
    ######################################

    awk -v chr="${CHR}" -v expected="${N_EXPECT}" '
    BEGIN {
        OFS = "\t"
    }

    NR == 1 {
        for (i = 1; i <= NF; i++) {
            field = tolower($i)
            if (field == "snp") snpcol = i
            if (field == "p") pcol = i
            if (field == "chr") chrcol = i
        }

        if (!snpcol || !pcol || !chrcol) {
            print "ERROR: Invalid GWAS header" > "/dev/stderr"
            bad = 1
            exit 1
        }

        print "SNP", "P"
        next
    }

    {
        snp = $(snpcol)
        p = $(pcol)

        if ($(chrcol)+0 != chr) {
            print "ERROR: Wrong chromosome" > "/dev/stderr"
            bad = 1
            exit 1
        }

        if (seen[snp]++) {
            print "ERROR: Duplicate SNP in GWAS" > "/dev/stderr"
            bad = 1
            exit 1
        }

        if (p !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
            print "ERROR: Invalid GWAS P-value:", p > "/dev/stderr"
            bad = 1
            exit 1
        }

        if (p+0 < 0 || p+0 > 1) {
            print "ERROR: GWAS P-value outside [0,1]" > "/dev/stderr"
            bad = 1
            exit 1
        }

        # Handle numerical underflow.
        if (p+0 == 0) {
            p = "1e-300"
            zero_p++
        }

        print snp, p
        n++
    }

    END {
        if (bad) exit 1

        if (n != expected) {
            print "ERROR: GWAS SNP count mismatch:", n, expected > "/dev/stderr"
            exit 1
        }

        if (zero_p > 0) {
            print "WARNING: P-values floored to 1e-300:", zero_p > "/dev/stderr"
        }
    }
    ' "${GWAS}" > "${ASSOC}"

    echo "PASS: GWAS input SNP count ${N_EXPECT}"

    ######################################
    # Imputed-only clumping
    ######################################

    plink \
        --bfile "${BFILE}" \
        --chr "${CHR}" \
        --clump "${ASSOC}" \
        --clump-snp-field SNP \
        --clump-field P \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 "${R2}" \
        --clump-kb "${WINDOW_KB}" \
        --threads 4 \
        --memory 24000 \
        --out "${PREFIX}"

    ######################################
    # Extract independent index SNPs
    ######################################

    [[ -s "${PREFIX}.clumped" ]] || {
        echo "ERROR: No clumped output" >&2
        exit 1
    }

    awk -v chr="${CHR}" '
    BEGIN { OFS = "\t" }

    NR == 1 {
        if ($1 != "CHR" || $3 != "SNP" || $5 != "P") {
            print "ERROR: Unexpected clumped header" > "/dev/stderr"
            bad = 1
            exit 1
        }

        print "SNP", "P"
        next
    }

    NF == 0 { next }

    {
        if ($1+0 != chr || $3 == "NONE" || $3 == "NA") {
            print "ERROR: Invalid clumped index SNP" > "/dev/stderr"
            bad = 1
            exit 1
        }

        if (seen[$3]++) {
            print "ERROR: Duplicate index SNP" > "/dev/stderr"
            bad = 1
            exit 1
        }

        print $3, $5
        n++
    }

    END {
        if (bad || n < 1) {
            print "ERROR: No valid index SNPs" > "/dev/stderr"
            exit 1
        }
    }
    ' "${PREFIX}.clumped" > "${RANKED}"

    # Sort index SNPs by GWAS P-value.
    {
        head -n 1 "${RANKED}"
        tail -n +2 "${RANKED}" \
            | LC_ALL=C sort -g -k2,2 -k1,1
    } > "${RANKED}.tmp"

    mv "${RANKED}.tmp" "${RANKED}"

    awk 'NR > 1 {print $1}' "${RANKED}" \
        > "${INDEX_LIST}"

    N_INDEX=$(wc -l < "${INDEX_LIST}")

    if (( N_INDEX < 1 || N_INDEX > N_EXPECT )); then
        echo "ERROR: Invalid index SNP count" >&2
        exit 1
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${OUTER}" "${INNER}" "${CHR}" \
        "${TRAIT}" "${N_SAMPLES}" "${N_EXPECT}" \
        "${N_INDEX}" "${R2}" "${WINDOW_KB}" \
        >> "${SUMMARY}"

    echo "Trait: ${TRAIT}"
    echo "Input GWAS SNPs: ${N_EXPECT}"
    echo "Clumped index SNPs: ${N_INDEX}"
    echo "PASS: ${TRAIT} clumping"

done

printf 'DONE Step8.3 outer_%s inner_%s chr_%s\n' \
    "${OUTER}" "${INNER}" "${CHR_PAD}" \
    > "${FLAG}"

echo "========================================"
echo "DONE Step8.3"
echo "========================================"
