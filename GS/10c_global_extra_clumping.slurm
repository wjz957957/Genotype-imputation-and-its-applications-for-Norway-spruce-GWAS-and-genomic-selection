#!/bin/bash
#SBATCH -J 10c_global_extra
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-19%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/10c_global_extra_%A_%a.out
#SBATCH -e logs/10c_global_extra_%A_%a.err

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
GWAS_ROOT="${WORK}/10_outer_test/03_imputed_GWAS"
INPUT_ROOT="${WORK}/10_outer_test/02_imputed_inputs"
OUTROOT="${WORK}/10_outer_test/04_imputed_clumping"
CHIP_LIST="${WORK}/01_reference/chip.marker"
CV_ROOT="${WORK}/02_CV"

########################################
# 1. Task mapping: 10 outer x 2 traits
########################################

TASK=${SLURM_ARRAY_TASK_ID:?ERROR: submit as a Slurm array}

if (( TASK < 0 || TASK > 19 )); then
    echo "ERROR: array task must be 0..19" >&2
    exit 1
fi

OUTER=$(printf '%02d' "$((TASK / 2 + 1))")

if (( TASK % 2 == 0 )); then
    TRAIT=BB
else
    TRAIT=H6
fi

BFILE="${INPUT_ROOT}/outer_${OUTER}/plink/imputed_train_qc"
CHR_COUNTS="${INPUT_ROOT}/outer_${OUTER}/imp_chr_counts.tsv"
TRAIN_IDS="${CV_ROOT}/outer_${OUTER}/outer_train.samples"
OUT="${OUTROOT}/outer_${OUTER}/${TRAIT}"

printf 'Step10c: task=%s outer=%s trait=%s\n' \
    "$TASK" "$OUTER" "$TRAIT"

printf 'GWAS source: %s\nGenotype source: %s\n' \
    "$GWAS_ROOT" "$BFILE"

########################################
# 2. Validate inputs
########################################

for f in \
    "$CHIP_LIST" \
    "$CHR_COUNTS" \
    "$TRAIN_IDS" \
    "${BFILE}.bed" \
    "${BFILE}.bim" \
    "${BFILE}.fam"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing/empty input: $f" >&2
        exit 1
    fi
done

# Remove legacy output for THIS outer/trait.
# An incomplete task will never receive a completion flag.

mkdir -p "${OUTROOT}/outer_${OUTER}"
rm -rf -- "$OUT"
mkdir -p "$OUT/clump" "$OUT/marker_sets"

GWAS_P="$OUT/all_gwas.p"
EXTRA_P="$OUT/imputed_only.p"
FILTER_QC="$OUT/filter_qc.tsv"

N_BIM=$(wc -l < "${BFILE}.bim")
N_TRAIN=$(wc -l < "${BFILE}.fam")
N_EXPECTED_TRAIN=$(wc -l < "$TRAIN_IDS")

if (( N_TRAIN != N_EXPECTED_TRAIN ||
      N_TRAIN < 700 ||
      N_TRAIN > 900 )); then
    echo "ERROR: outer-training sample count inconsistent" >&2
    exit 1
fi

awk '
NF != 6 || $1 != $2 { bad=1 }
{ print $2 }
END { if(bad) exit 1 }
' "${BFILE}.fam" \
    | LC_ALL=C sort > "$OUT/genotype.samples.sorted"

LC_ALL=C sort "$TRAIN_IDS" \
    > "$OUT/expected.samples.sorted"

if ! cmp -s \
    "$OUT/genotype.samples.sorted" \
    "$OUT/expected.samples.sorted"
then
    echo "ERROR: genotype individuals do not match outer train" >&2
    exit 1
fi

########################################
# 3. Merge 12 chromosome GWAS reports
#
# GCTA MLMA:
# Chr SNP bp A1 A2 Freq b se p
#
# SNP = column 2
# P   = column 9
########################################

printf 'SNP\tP\n' > "$GWAS_P"

for CHR in $(seq 1 12)
do
    CHR_PAD=$(printf '%02d' "$CHR")

    MLMA="${GWAS_ROOT}/outer_${OUTER}/chr_${CHR_PAD}/${TRAIT}/${TRAIT}_chr${CHR_PAD}.mlma"

    if [[ ! -s "$MLMA" ]]; then
        echo "ERROR: missing GWAS: $MLMA" >&2
        exit 1
    fi

    EXPECTED_CHR=$(awk -v chr="$CHR" '
        $1 == chr { n=$2; hits++ }

        END {
            if(hits != 1 || n !~ /^[0-9]+$/ || n < 1000) exit 1
            print n
        }
    ' "$CHR_COUNTS")

    awk -v chr="$CHR" -v expected="$EXPECTED_CHR" '
        BEGIN { OFS="\t" }

        FNR == 1 {
            if (NF != 9 ||
                $1 != "Chr" ||
                $2 != "SNP" ||
                $9 != "p") {
                print "ERROR: invalid MLMA header: " FILENAME > "/dev/stderr"
                bad=1
                exit 1
            }
            next
        }

        {
            if (NF != 9 ||
                $1 != chr ||
                index($2, chr ":") != 1) {
                print "ERROR: invalid GWAS row at " FILENAME ":" FNR > "/dev/stderr"
                bad=1
                exit 1
            }

            print $2, $9
            n++
        }

        END {
            if (n != expected) {
                print "ERROR: GWAS marker count mismatch: chr=" chr ", got=" n ", expected=" expected > "/dev/stderr"
                bad=1
            }

            if(bad) exit 1
        }
    ' "$MLMA" >> "$GWAS_P"

    echo "PASS chr${CHR_PAD}: ${EXPECTED_CHR} GWAS rows"
done

########################################
# 4. Validate globally merged SNPs
########################################

N_GWAS=$(($(wc -l < "$GWAS_P") - 1))

if (( N_GWAS != N_BIM )); then
    echo "ERROR: merged GWAS rows ${N_GWAS} != BIM SNPs ${N_BIM}" >&2
    exit 1
fi

# Check:
# - exactly two columns
# - header SNP / P
# - no duplicate SNP IDs
# - every SNP exists in imputed BIM

awk '
    FILENAME == ARGV[1] {
        bim[$2]=1
        next
    }

    FILENAME == ARGV[2] && FNR == 1 {
        if(NF != 2 || $1 != "SNP" || $2 != "P") {
            bad=1
            print "ERROR: GWAS_P header" > "/dev/stderr"
        }
        next
    }

    FILENAME == ARGV[2] {
        if(NF != 2 || !($1 in bim) || seen[$1]++) {
            if(++bad <= 5) {
                print "ERROR: invalid/duplicate GWAS SNP: " $1 > "/dev/stderr"
            }
        }
        n++
    }

    END {
        if(bad || n < 1000) exit 1
    }
' "${BFILE}.bim" "$GWAS_P"

########################################
# 5. Remove ALL original chip SNPs
#
# Not best chip!
# Do not filter imputed-chip LD.
#
# Invalid/nonfinite P values are skipped.
########################################

awk -v qc="$FILTER_QC" '
    BEGIN {
        OFS="\t"
        print "SNP", "P"
    }

    FILENAME == ARGV[1] {
        if(NF != 1 || chip[$1]++) {
            bad=1
            print "ERROR: invalid/duplicate chip SNP " $1 > "/dev/stderr"
        }
        next
    }

    FILENAME == ARGV[2] && FNR == 1 {
        if($1 != "SNP" || $2 != "P" || NF != 2) {
            bad=1
        }
        next
    }

    FILENAME == ARGV[2] {
        total++
        p=$2

        if (p !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ ||
            p+0 < 0 ||
            p+0 > 1) {
            invalid++
            next
        }

        valid++

        if ($1 in chip) {
            removed++
            next
        }

        if (p+0 == 0) {
            p="1e-300"
            floored++
        }

        print $1, p
        retained++
    }

    END {
        print "N_total_GWAS\tN_valid_P\tN_invalid_P\tN_chip_removed\tN_extra_retained\tN_zero_P_floored" > qc

        printf "%d\t%d\t%d\t%d\t%d\t%d\n", total, valid, invalid, removed, retained, floored >> qc

        if(bad ||
           total != valid + invalid ||
           valid != retained + removed ||
           retained < 1000) {
            exit 1
        }
    }
' "$CHIP_LIST" "$GWAS_P" > "$EXTRA_P"

read -r \
    N_FILTERED_TOTAL \
    N_VALID_P \
    N_INVALID_P \
    N_REMOVED \
    N_EXTRA \
    N_FLOORED \
    < <(tail -n 1 "$FILTER_QC")

if (( N_FILTERED_TOTAL != N_GWAS ||
      N_VALID_P != N_EXTRA + N_REMOVED )); then
    echo "ERROR: SNP accounting failed" >&2
    exit 1
fi

########################################
# 6. Independent chip-overlap check
########################################

awk '
    FILENAME == ARGV[1] {
        chip[$1]=1
        next
    }

    FILENAME == ARGV[2] && FNR > 1 {
        if($1 in chip) bad++
    }

    END {
        if(bad) {
            print "ERROR: chip SNPs remain in extra pool: " bad > "/dev/stderr"
            exit 1
        }
    }
' "$CHIP_LIST" "$EXTRA_P"

printf 'GWAS=%d valid_P=%d invalid_P=%d matched_chip_removed=%d extra=%d\n' \
    "$N_GWAS" \
    "$N_VALID_P" \
    "$N_INVALID_P" \
    "$N_REMOVED" \
    "$N_EXTRA"

head -n 4 "$EXTRA_P"

########################################
# 7. Global LD clumping
#
# Match Step8:
# r2 = 0.2
# window = 1000 kb
#
# Outer-training genotypes only.
########################################

plink \
    --bfile "$BFILE" \
    --clump "$EXTRA_P" \
    --clump-p1 1 \
    --clump-p2 1 \
    --clump-r2 0.2 \
    --clump-kb 1000 \
    --threads 4 \
    --memory 32000 \
    --out "$OUT/clump/extra"

CLUMP="$OUT/clump/extra.clumped"

if [[ ! -s "$CLUMP" ]]; then
    echo "ERROR: missing clumped output" >&2
    exit 1
fi

########################################
# 8. Extract clumping INDEX SNPs
#
# PLINK .clumped column 3 = SNP
#
# Do not select top SNPs directly
# from all GWAS variants.
########################################

INDEX_IDS="$OUT/clump_index.snps"

awk '
    NF >= 5 &&
    $1 ~ /^[0-9]+$/ &&
    $3 != "NONE" &&
    $3 != "." {
        if(seen[$3]++) {
            bad=1
            print "ERROR: duplicate clump index " $3 > "/dev/stderr"
        }

        print $3
        n++
    }

    END {
        if(bad || n < 1000) {
            print "ERROR: too few/invalid clump indices: " n > "/dev/stderr"
            exit 1
        }
    }
' "$CLUMP" > "$INDEX_IDS"

N_INDEX=$(wc -l < "$INDEX_IDS")

########################################
# 9. Retrieve original P values
#    for all clump index SNPs
########################################

UNSORTED="$OUT/clump_index_with_p.unsorted.tsv"

awk -v expected="$N_INDEX" '
    BEGIN { OFS="\t" }

    FILENAME == ARGV[1] {
        index_snp[$1]=1
        next
    }

    FILENAME == ARGV[2] && FNR == 1 {
        next
    }

    FILENAME == ARGV[2] && ($1 in index_snp) {
        print $1, $2
        n++
    }

    END {
        if(n != expected) {
            print "ERROR: clump index association P join mismatch: " n "/" expected > "/dev/stderr"
            exit 1
        }
    }
' "$INDEX_IDS" "$EXTRA_P" > "$UNSORTED"

########################################
# 10. GLOBAL P-value ranking
########################################

RANKED="$OUT/ranked_extra_snps.tsv"

{
    printf 'SNP\tP\n'

    LC_ALL=C sort \
        -t $'\t' \
        -k2,2g \
        -k1,1 \
        "$UNSORTED"

} > "$RANKED"

########################################
# 11. Generate QTN sets
#
# Counts refer to whole genome,
# NOT individual chromosomes.
########################################

: > "$OUT/marker_sets/extra_0.snps"

for N in 1 5 10 20 50 100 200 500 1000
do
    awk -v n="$N" '
        NR > 1 && NR <= n+1 { print $1 }
    ' "$RANKED" > "$OUT/marker_sets/extra_${N}.snps"

    GOT=$(wc -l < "$OUT/marker_sets/extra_${N}.snps")

    if (( GOT != N )); then
        echo "ERROR: extra_${N} contains ${GOT}, expected ${N}" >&2
        exit 1
    fi
done

awk '
    NR > 1 { print $1 }
' "$RANKED" > "$OUT/marker_sets/extra_LD_ALL.snps"

if (( $(wc -l < "$OUT/marker_sets/extra_LD_ALL.snps") != N_INDEX )); then
    echo "ERROR: LD_ALL count mismatch" >&2
    exit 1
fi

########################################
# 12. Summary
########################################

{
    printf 'Outer\tTrait\tN_train\tN_GWAS\tN_valid_P\tN_invalid_P\tN_chip_removed\tN_extra_pool\tN_LD_index\tClump_r2\tClump_kb\n'

    printf '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t0.2\t1000\n' \
        "$OUTER" \
        "$TRAIT" \
        "$N_TRAIN" \
        "$N_GWAS" \
        "$N_VALID_P" \
        "$N_INVALID_P" \
        "$N_REMOVED" \
        "$N_EXTRA" \
        "$N_INDEX"

} > "$OUT/summary.tsv"

########################################
# 13. Completion
########################################

printf 'DONE Step10c outer_%s trait_%s\n' \
    "$OUTER" "$TRAIT" > "$OUT/complete.flag"

echo "PASS: outer=${OUTER} trait=${TRAIT} extra_indices=${N_INDEX}"
