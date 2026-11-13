#!/bin/bash
#SBATCH -J cf_chip_gwas
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-14%3
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/03_chip_gwas_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/03_chip_gwas_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"

SAMPLE_ROOT="${WORK}/Cross_family/cv_samples_v1"
GENO_ROOT="${WORK}/Cross_family/results/01_chip_inner_genotypes"
PHENO_ROOT="${WORK}/Cross_family/results/02_inner_phenotypes"
OUTROOT="${WORK}/Cross_family/results/03_inner_chip_gwas_clumping"

N_EXPECTED_SNP=43267
SCRIPT_VERSION=1
THREADS="${SLURM_CPUS_PER_TASK:-4}"

########################################
# 1. Software
########################################

module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1

########################################
# 2. Map array task to outer/inner fold
########################################

TASK="${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array job}"

if [[ ! "${TASK}" =~ ^[0-9]+$ ]] || (( TASK < 0 || TASK > 14 )); then
    echo "ERROR: Task ID must be between 0 and 14" >&2
    exit 1
fi

OUTER=$(printf '%02d' "$((TASK / 3 + 1))")
INNER=$(printf '%02d' "$((TASK % 3 + 1))")

SAMPLE_DIR="${SAMPLE_ROOT}/outer_${OUTER}/inner_${INNER}"
GENO="${GENO_ROOT}/outer_${OUTER}/inner_${INNER}/chip_inner_train"
PHENO_DIR="${PHENO_ROOT}/outer_${OUTER}/inner_${INNER}"

OUT="${OUTROOT}/outer_${OUTER}/inner_${INNER}"
GRM="${OUT}/grm/chip_train_grm"
GWAS_DIR="${OUT}/gwas"
CLUMP_DIR="${OUT}/clumping"
SET_DIR="${OUT}/marker_sets"

FLAG="${OUT}/complete.flag"
QC="${OUT}/gwas_qc.tsv"
SET_SUMMARY="${OUT}/marker_set_counts.tsv"

mkdir -p "${OUT}" "${OUT}/grm" "${GWAS_DIR}" "${CLUMP_DIR}" "${SET_DIR}"

echo "========================================"
echo "Cross-family inner-train chip GWAS"
echo "Task: ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "========================================"

########################################
# 3. Input files and completion check
########################################

TRAIN="${SAMPLE_DIR}/inner_train.samples"

for F in \
    "${TRAIN}" \
    "${GENO}.bed" "${GENO}.bim" "${GENO}.fam" \
    "${PHENO_DIR}/BB.train.pheno" \
    "${PHENO_DIR}/H6.train.pheno"
do
    [[ -s "${F}" ]] || {
        echo "ERROR: Missing or empty input: ${F}" >&2
        exit 1
    }
done

if [[ -s "${FLAG}" ]] &&
   grep -qx "version=${SCRIPT_VERSION}" "${FLAG}"; then
    echo "Previously completed output detected; skipping."
    exit 0
fi

rm -f "${FLAG}"

########################################
# 4. Validate training sample IDs
########################################

awk '
{
    if (NF != 1 || seen[$1]++) {
        print "ERROR: Invalid or duplicate training ID:", $1 > "/dev/stderr"
        bad=1
    }
}
END {
    if (NR == 0 || bad) exit 1
}
' "${TRAIN}"

N_TRAIN=$(wc -l < "${TRAIN}")
N_FAM=$(wc -l < "${GENO}.fam")
N_SNP=$(wc -l < "${GENO}.bim")

if (( N_FAM != N_TRAIN )); then
    echo "ERROR: Genotype/sample-list count mismatch" >&2
    exit 1
fi

if (( N_SNP != N_EXPECTED_SNP )); then
    echo "ERROR: Expected ${N_EXPECTED_SNP} SNPs; found ${N_SNP}" >&2
    exit 1
fi

awk '
$1 != $2 || seen[$2]++ {
    print "ERROR: Invalid or duplicate FAM ID:", $2 > "/dev/stderr"
    bad=1
}
END {if (bad) exit 1}
' "${GENO}.fam"

EXPECTED_IDS="${OUT}/expected_train_ids.sorted"
OBSERVED_IDS="${OUT}/observed_genotype_ids.sorted"

sort "${TRAIN}" > "${EXPECTED_IDS}"
awk '{print $2}' "${GENO}.fam" | sort > "${OBSERVED_IDS}"

if ! cmp -s "${EXPECTED_IDS}" "${OBSERVED_IDS}"; then
    echo "ERROR: Genotype IDs differ from inner_train.samples" >&2
    diff -u "${EXPECTED_IDS}" "${OBSERVED_IDS}" | head -n 30 || true
    exit 1
fi

echo "PASS: Training sample IDs match"
echo "Training samples: ${N_TRAIN}"
echo "Chip SNPs: ${N_SNP}"

########################################
# 5. Prepare and validate training phenotypes
########################################

for TRAIT in BB H6; do

    SOURCE="${PHENO_DIR}/${TRAIT}.train.pheno"
    TARGET="${OUT}/${TRAIT}.train.gcta.pheno"

    awk -v expected="${N_TRAIN}" '
    BEGIN {
        OFS="\t"
        numeric="^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$"
    }

    NF == 0 {next}

    {
        if (tolower($1) == "fid" && tolower($2) == "iid")
            next

        if (NF != 3 || $1 != $2 || seen[$2]++) {
            print "ERROR: Invalid phenotype row:", NR > "/dev/stderr"
            bad=1
            exit 1
        }

        if ($3 !~ numeric) {
            print "ERROR: Non-numeric phenotype:", $1, $3 > "/dev/stderr"
            bad=1
            exit 1
        }

        print $1, $2, $3
        n++
    }

    END {
        if (bad || n != expected) {
            print "ERROR: Phenotype count mismatch; expected",
                  expected, "observed", n+0 > "/dev/stderr"
            exit 1
        }
    }
    ' "${SOURCE}" > "${TARGET}"

    PHENO_IDS="${OUT}/${TRAIT}.phenotype_ids.sorted"
    awk '{print $2}' "${TARGET}" | sort > "${PHENO_IDS}"

    if ! cmp -s "${EXPECTED_IDS}" "${PHENO_IDS}"; then
        echo "ERROR: ${TRAIT} phenotype IDs differ from training IDs" >&2
        diff -u "${EXPECTED_IDS}" "${PHENO_IDS}" | head -n 30 || true
        exit 1
    fi

    echo "PASS: ${TRAIT} training phenotype IDs and values"

done

########################################
# 6. Build or reuse training GRM
########################################

if [[ ! -s "${GRM}.grm.bin" ||
      ! -s "${GRM}.grm.N.bin" ||
      ! -s "${GRM}.grm.id" ]]; then

    gcta64 \
        --bfile "${GENO}" \
        --make-grm \
        --thread-num "${THREADS}" \
        --out "${GRM}"
fi

for EXT in grm.bin grm.N.bin grm.id; do
    [[ -s "${GRM}.${EXT}" ]] || {
        echo "ERROR: Missing GRM file ${GRM}.${EXT}" >&2
        exit 1
    }
done

if [[ $(wc -l < "${GRM}.grm.id") -ne ${N_TRAIN} ]]; then
    echo "ERROR: GRM sample count mismatch" >&2
    exit 1
fi

GRM_IDS="${OUT}/grm_ids.sorted"
awk '{print $2}' "${GRM}.grm.id" | sort > "${GRM_IDS}"

if ! cmp -s "${EXPECTED_IDS}" "${GRM_IDS}"; then
    echo "ERROR: GRM IDs differ from training IDs" >&2
    exit 1
fi

echo "PASS: Training GRM"

########################################
# 7. GWAS, P-value QC and LD clumping
########################################

printf 'Trait\tTotal_SNPs\tValid_P\tNonfinite_P\tZero_P_Floored\n' > "${QC}"
printf 'Trait\tMarker_Set\tRequested_N\tActual_N\n' > "${SET_SUMMARY}"

for TRAIT in BB H6; do

    echo "========================================"
    echo "Trait: ${TRAIT}"
    echo "========================================"

    PHENO="${OUT}/${TRAIT}.train.gcta.pheno"
    PREFIX="${GWAS_DIR}/chip_${TRAIT}_gwas"
    GWAS="${PREFIX}.loco.mlma"

    ####################################
    # 7.1 GWAS
    ####################################

    if [[ ! -s "${GWAS}" ]]; then
        gcta64 \
            --mlma-loco \
            --bfile "${GENO}" \
            --grm "${GRM}" \
            --pheno "${PHENO}" \
            --thread-num "${THREADS}" \
            --out "${PREFIX}"
    else
        echo "REUSE: Existing GWAS ${GWAS}"
    fi

    [[ -s "${GWAS}" ]] || {
        echo "ERROR: Missing GWAS output: ${GWAS}" >&2
        exit 1
    }

    ####################################
    # 7.2 Validate P-values
    ####################################

    TC="${CLUMP_DIR}/${TRAIT}"
    mkdir -p "${TC}/clump"

    PFILE="${TC}/gwas_p.txt"
    COUNTS="${TC}/pvalue_counts.txt"

    awk -v expected="${N_EXPECTED_SNP}" -v counts="${COUNTS}" '
    BEGIN {
        OFS="\t"
        numeric="^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$"
        print "SNP", "P"
    }

    NR == 1 {
        if (tolower($1) != "chr" ||
            tolower($2) != "snp" ||
            tolower($9) != "p") {
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

        if (low ~ /^[-+]?nan(\([^)]*\))?$/ ||
            low ~ /^[-+]?inf(inity)?$/ ||
            low == "na" || low == "n/a" ||
            low == "." || low == "-na") {
            invalid++
            next
        }

        if (p !~ numeric || p+0 < 0 || p+0 > 1) {
            print "ERROR: Invalid P-value for", $2, p > "/dev/stderr"
            bad=1
            exit 1
        }

        # Floor numerical zero only to avoid underflow in clumping.
        if (p+0 == 0) {
            p="1e-300"
            zero++
        }

        print $2, p
        valid++
    }

    END {
        if (bad || total != expected ||
            valid < 1 || valid+invalid != total) {
            print "ERROR: GWAS P-value QC failed:",
                  "total", total+0,
                  "valid", valid+0,
                  "invalid", invalid+0,
                  "expected", expected > "/dev/stderr"
            exit 1
        }

        printf "%d\t%d\t%d\t%d\n",
               total, valid, invalid, zero+0 > counts
        close(counts)
    }
    ' "${GWAS}" > "${PFILE}"

    read -r TOTAL VALID INVALID ZERO < "${COUNTS}"

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${TRAIT}" "${TOTAL}" "${VALID}" "${INVALID}" "${ZERO}" \
        >> "${QC}"

    echo "GWAS SNPs: ${TOTAL}"
    echo "Valid P-values: ${VALID}"
    echo "Nonfinite P-values excluded: ${INVALID}"
    echo "Zero P-values floored: ${ZERO}"

    ####################################
    # 7.3 LD clumping
    ####################################

    plink \
        --bfile "${GENO}" \
        --clump "${PFILE}" \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.2 \
        --clump-kb 50 \
        --allow-extra-chr \
        --threads "${THREADS}" \
        --memory 24000 \
        --out "${TC}/clump/chip"

    CLUMPED="${TC}/clump/chip.clumped"

    [[ -s "${CLUMPED}" ]] || {
        echo "ERROR: Missing or empty clumping output" >&2
        exit 1
    }

    ####################################
    # 7.4 Rank independent index SNPs
    ####################################

    RANKED="${TC}/ranked_chip_snps.tsv"

    awk '
    NR == 1 {
        if (toupper($1) != "CHR" ||
            toupper($3) != "SNP" ||
            toupper($5) != "P") {
            print "ERROR: Unexpected .clumped header" > "/dev/stderr"
            bad=1
            exit 1
        }
        next
    }

    NF == 0 {next}

    {
        if ($3 == "NONE" || $3 == "NA" || seen[$3]++) {
            print "ERROR: Invalid or duplicate index SNP:", $3 > "/dev/stderr"
            bad=1
            exit 1
        }

        print $3 "\t" $5
        n++
    }

    END {
        if (bad || n < 1) {
            print "ERROR: No valid clumped index SNPs" > "/dev/stderr"
            exit 1
        }
    }
    ' "${CLUMPED}" | sort -k2,2g -k1,1 > "${RANKED}"

    N_RANKED=$(wc -l < "${RANKED}")
    echo "Independent index SNPs: ${N_RANKED}"

    ####################################
    # 7.5 Generate candidate chip sets
    ####################################

    for N in \
        50 100 200 300 400 500 \
        1000 2000 3000 4000 5000 \
        8000 10000 20000 30000 40000
    do
        FILE="${SET_DIR}/${TRAIT}_chip_${N}.snps"

        head -n "${N}" "${RANKED}" |
            awk '{print $1}' > "${FILE}"

        ACTUAL=$(wc -l < "${FILE}")

        printf '%s\t%s\t%s\t%s\n' \
            "${TRAIT}" "chip_${N}" "${N}" "${ACTUAL}" \
            >> "${SET_SUMMARY}"
    done

    # All LD-clumped index SNPs, ranked by P-value.
    awk '{print $1}' "${RANKED}" \
        > "${SET_DIR}/${TRAIT}_chip_LD_ALL.snps"

    printf '%s\t%s\t%s\t%s\n' \
        "${TRAIT}" "chip_LD_ALL" "ALL_CLUMPED" "${N_RANKED}" \
        >> "${SET_SUMMARY}"

    # All original chip SNPs, without GWAS ranking or clumping.
    awk '{print $2}' "${GENO}.bim" \
        > "${SET_DIR}/${TRAIT}_chip_ALL.snps"

    printf '%s\t%s\t%s\t%s\n' \
        "${TRAIT}" "chip_ALL" "${N_EXPECTED_SNP}" "${N_EXPECTED_SNP}" \
        >> "${SET_SUMMARY}"

    ####################################
    # 7.6 Validate every candidate set
    ####################################

    for LIST in "${SET_DIR}/${TRAIT}_"*.snps; do

        [[ -s "${LIST}" ]] || {
            echo "ERROR: Empty candidate set: ${LIST}" >&2
            exit 1
        }

        N_LIST=$(wc -l < "${LIST}")
        N_UNIQUE=$(sort -u "${LIST}" | wc -l)

        if (( N_LIST != N_UNIQUE )); then
            echo "ERROR: Duplicate SNP IDs in ${LIST}" >&2
            exit 1
        fi

        awk '
        NR == FNR {
            available[$2]=1
            next
        }
        !($1 in available) {
            print "ERROR: SNP absent from genotype BIM:", $1 > "/dev/stderr"
            bad=1
        }
        END {if (bad) exit 1}
        ' "${GENO}.bim" "${LIST}"

    done

    echo "PASS: ${TRAIT} candidate marker sets"
done

########################################
# 8. Publish completion flag
########################################

{
    echo "version=${SCRIPT_VERSION}"
    echo "outer=${OUTER}"
    echo "inner=${INNER}"
    echo "n_train=${N_TRAIN}"
    echo "n_chip_snps=${N_SNP}"
    echo "sample_ID_check=PASS"
    echo "phenotype_ID_check=PASS"
    echo "GRM_check=PASS"
    echo "GWAS_P_QC=PASS"
    echo "candidate_SNP_check=PASS"
} > "${FLAG}.tmp"

mv "${FLAG}.tmp" "${FLAG}"

echo "========================================"
echo "OVERALL PASS: outer_${OUTER}/inner_${INNER}"
echo "Training samples: ${N_TRAIN}"
echo "Chip SNPs: ${N_SNP}"
echo "Traits: BB, H6"
echo "========================================"
