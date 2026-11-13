#!/bin/bash
#SBATCH -J cf09_imp_gwas
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-59%6
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09d_imputed_gwas_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09d_imputed_gwas_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 09d at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"

module load bioinfo-tools
ml gcta/1.94.1
GCTA_BIN="$(command -v gcta64 || command -v gcta || true)"
[[ -n "$GCTA_BIN" ]] || { echo "ERROR: GCTA executable not found" >&2; exit 127; }

THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1

INPUT_ROOT="$WORK/Cross_family/results/09_outer_imputed_gwas_inputs"
PHENO_ROOT="$WORK/Cross_family/results/09_outer_inputs"
OUT_ROOT="$WORK/Cross_family/results/09_outer_imputed_gwas"
TASK="${SLURM_ARRAY_TASK_ID:?Run through sbatch array tasks 0-59}"
(( TASK >= 0 && TASK < 60 )) || { echo "Invalid task id: $TASK" >&2; exit 2; }

OUTER="$(printf '%02d' "$(( TASK / 12 + 1 ))")"
CHR="$(( TASK % 12 + 1 ))"
PAD="$(printf '%02d' "$CHR")"
INPUT="$INPUT_ROOT/outer_${OUTER}"
PHENO_DIR="$PHENO_ROOT/outer_${OUTER}"
OUT="$OUT_ROOT/outer_${OUTER}/chr_${PAD}"
BFILE="$INPUT/imputed_train_qc"
GRM_ALL="$INPUT/chip_grm/chip_all"
GRM_CHR="$INPUT/chip_grm/chip_chr${PAD}"
SCRIPT_VERSION=1

mkdir -p "$OUT"

# 09c must be complete, and the outer-test phenotype must never be read.
[[ -s "$INPUT/complete.flag" ]] || { echo "ERROR: missing 09c complete.flag: $INPUT/complete.flag" >&2; exit 1; }
for ENTRY in \
    "outer=${OUTER}" \
    "status=PASS" \
    "training_only=PASS" \
    "outer_test_phenotype_read=NO" \
    "full_chip_grm=PASS" \
    "chromosome_chip_grms=12/12"
do
    grep -qx "$ENTRY" "$INPUT/complete.flag" || {
        echo "ERROR: 09c flag missing or mismatched entry: $ENTRY" >&2
        exit 1
    }
done

for F in "${BFILE}.bed" "${BFILE}.bim" "${BFILE}.fam" \
         "${GRM_ALL}.grm.bin" "${GRM_ALL}.grm.N.bin" "${GRM_ALL}.grm.id" \
         "${GRM_CHR}.grm.bin" "${GRM_CHR}.grm.N.bin" "${GRM_CHR}.grm.id" \
         "$PHENO_DIR/BB.train.pheno" "$PHENO_DIR/H6.train.pheno"; do
    [[ -s "$F" ]] || { echo "ERROR: missing/empty prerequisite: $F" >&2; exit 1; }
done

cmp -s "${GRM_ALL}.grm.id" "${GRM_CHR}.grm.id" || {
    echo "ERROR: full and chromosome GRM sample order mismatch: outer_${OUTER}, chr${PAD}" >&2
    exit 1
}

N_TRAIN="$(wc -l < "${BFILE}.fam")"
N_GRM="$(wc -l < "${GRM_ALL}.grm.id")"
N_EXPECT="$(awk -v chr="$CHR" '$1==chr{n++} END{print n+0}' "${BFILE}.bim")"
N_TRAIN_FLAG="$(awk -F= '$1=="n_training"{print $2}' "$INPUT/complete.flag")"
(( N_TRAIN >= 1 && N_TRAIN == N_GRM && N_TRAIN == N_TRAIN_FLAG )) || {
    echo "ERROR: training sample count mismatch: fam=$N_TRAIN GRM=$N_GRM flag=$N_TRAIN_FLAG" >&2
    exit 1
}
(( N_EXPECT >= 1000 )) || {
    echo "ERROR: too few imputed SNPs on chr${PAD}: ${N_EXPECT}" >&2
    exit 1
}

# Guard against accidentally reusing completed results with changed inputs.
# Large binary files are checked by size/mtime rather than re-hashed for every array task.
SIG="$OUT/input_signature.txt"
SIG_NOW="$OUT/input_signature.now"
{
    sha256sum "$INPUT/complete.flag" "${BFILE}.fam" \
        "${GRM_ALL}.grm.id" "${GRM_CHR}.grm.id" \
        "$PHENO_DIR/BB.train.pheno" "$PHENO_DIR/H6.train.pheno"
    stat -c 'file=%n size=%s mtime=%Y' \
        "${BFILE}.bed" "${BFILE}.bim" \
        "${GRM_ALL}.grm.bin" "${GRM_ALL}.grm.N.bin" \
        "${GRM_CHR}.grm.bin" "${GRM_CHR}.grm.N.bin"
} > "$SIG_NOW"

if [[ -s "$SIG" ]]; then
    cmp -s "$SIG" "$SIG_NOW" || {
        rm -f "$SIG_NOW"
        echo "ERROR: 09d inputs changed for outer_${OUTER}/chr_${PAD}; do not reuse old results." >&2
        exit 1
    }
    rm -f "$SIG_NOW"
else
    if [[ -e "$OUT/complete.flag" || -e "$OUT/BB_chr${PAD}.mlma" || -e "$OUT/H6_chr${PAD}.mlma" ]]; then
        rm -f "$SIG_NOW"
        echo "ERROR: existing outputs lack input_signature.txt; inspect $OUT before rerunning." >&2
        exit 1
    fi
    mv "$SIG_NOW" "$SIG"
fi

# A completed task is reusable only if its flag and both GWAS files are present.
if [[ -s "$OUT/complete.flag" ]] &&
   grep -qx "status=PASS" "$OUT/complete.flag" &&
   grep -qx "version=${SCRIPT_VERSION}" "$OUT/complete.flag" &&
   [[ -s "$OUT/gwas_qc.tsv" && -s "$OUT/BB_chr${PAD}.mlma" && -s "$OUT/H6_chr${PAD}.mlma" ]]; then
    echo "[SKIP] outer_${OUTER}/chr_${PAD} already complete with matching input signature"
    exit 0
fi

rm -f "$OUT/complete.flag" "$OUT/gwas_qc.tsv"
printf 'Outer\tTrait\tChromosome\tExpected_SNPs\tTested_SNPs\tValid_P\tInvalid_P\n' > "$OUT/gwas_qc.tsv"

for TRAIT in BB H6; do
    PHENO="$PHENO_DIR/${TRAIT}.train.pheno"
    PREFIX="$OUT/${TRAIT}_chr${PAD}"
    RESULT="${PREFIX}.mlma"

    # Remove partial/stale files for this chromosome before rerunning GCTA.
    rm -f "$RESULT" "${PREFIX}.log"
    echo "===== GCTA MLMA-LOCO: outer_${OUTER}, ${TRAIT}, chr${CHR} ====="
    "$GCTA_BIN" --mlma --bfile "$BFILE" --autosome-num 12 --chr "$CHR" \
        --grm "$GRM_ALL" --mlma-subtract-grm "$GRM_CHR" \
        --pheno "$PHENO" --thread-num "$THREADS" --out "$PREFIX"

    [[ -s "$RESULT" ]] || { echo "ERROR: GCTA result missing: $RESULT" >&2; exit 1; }

    awk -v outer="$OUTER" -v trait="$TRAIT" -v chr="$CHR" -v expected="$N_EXPECT" '
        function is_probability(x) {
            if (x !~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ &&
                x !~ /^[.][0-9]+([eE][-+]?[0-9]+)?$/) return 0
            return (x+0 >= 0 && x+0 <= 1)
        }
        NR==1 {
            if (NF<9 || $1!="Chr" || $2!="SNP" || $9!="p") {
                print "ERROR: unexpected GCTA MLMA header" > "/dev/stderr"
                bad=1
                exit 1
            }
            next
        }
        {
            if (NF<9 || $1!=chr) {
                print "ERROR: invalid chromosome/columns in GCTA result" > "/dev/stderr"
                bad=1
                exit 1
            }
            tested++
            if (is_probability($9)) valid++
            else invalid++
        }
        END {
            if (bad) exit 1
            if (tested<1000 || valid<1000 || tested>expected || valid+invalid!=tested) {
                printf "ERROR: GWAS result QC failed: tested=%d valid=%d invalid=%d expected=%d\n", \
                    tested+0,valid+0,invalid+0,expected > "/dev/stderr"
                exit 1
            }
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\n", \
                outer,trait,chr,expected,tested,valid,invalid
        }
    ' "$RESULT" >> "$OUT/gwas_qc.tsv"
done

cat > "$OUT/complete.flag.tmp" <<EOF_FLAG
status=PASS
version=${SCRIPT_VERSION}
outer=${OUTER}
chromosome=${CHR}
chromosome_padded=${PAD}
n_training=${N_TRAIN}
n_imputed_snps_on_chromosome=${N_EXPECT}
traits=BB,H6
model=GCTA_MLMA_with_mlma_subtract_grm
outer_test_phenotype_read=NO
EOF_FLAG
mv "$OUT/complete.flag.tmp" "$OUT/complete.flag"

echo "PASS 09d outer_${OUTER}/chr_${PAD}: training=${N_TRAIN}, expected_imputed_SNPs=${N_EXPECT}, traits=BB,H6"
