#!/bin/bash
#SBATCH -J cf_qtn_mlma
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-179%6
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06b_qtn_mlma_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06b_qtn_mlma_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"
module load bioinfo-tools
ml gcta/1.94.1
GCTA_BIN="$(command -v gcta64 || command -v gcta || true)"
[[ -n "$GCTA_BIN" ]] || { echo "ERROR: GCTA executable not found" >&2; exit 127; }
CPUS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="$CPUS"
export OPENBLAS_NUM_THREADS=1

INPUT_ROOT="$WORK/Cross_family/results/06_inner_imputed_gwas_inputs"
PHENO_ROOT="$WORK/Cross_family/results/02_inner_phenotypes"
OUT_ROOT="$WORK/Cross_family/results/06_inner_imputed_gwas"
TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
(( TASK >= 0 && TASK < 180 )) || { echo "Invalid task id: $TASK" >&2; exit 2; }
FOLD_IDX=$((TASK / 12))
CHR=$((TASK % 12 + 1))
OUTER=$(printf '%02d' "$((FOLD_IDX / 3 + 1))")
INNER=$(printf '%02d' "$((FOLD_IDX % 3 + 1))")
PAD=$(printf '%02d' "$CHR")
FOLD="outer_${OUTER}/inner_${INNER}"
INPUT="$INPUT_ROOT/$FOLD"
OUT="$OUT_ROOT/$FOLD/chr_${PAD}"
mkdir -p "$OUT"

[[ -s "$INPUT/complete.flag" ]] && grep -qx 'status=PASS' "$INPUT/complete.flag" || { echo "ERROR: Stage 06a incomplete: $INPUT/complete.flag" >&2; exit 1; }
BFILE="$INPUT/imputed_train_qc"
GRM_ALL="$INPUT/chip_grm/chip_all"
GRM_CHR="$INPUT/chip_grm/chip_chr${PAD}"
for f in "${BFILE}.bed" "${BFILE}.bim" "${BFILE}.fam"; do [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done
for p in "$GRM_ALL" "$GRM_CHR"; do for ext in grm.bin grm.N.bin grm.id; do [[ -s "${p}.${ext}" ]] || { echo "ERROR: missing ${p}.${ext}" >&2; exit 1; }; done; done
cmp -s "${GRM_ALL}.grm.id" "${GRM_CHR}.grm.id" || { echo "ERROR: full and chr GRM sample order mismatch" >&2; exit 1; }
N_TRAIN=$(wc -l < "${BFILE}.fam")
[[ $(wc -l < "${GRM_ALL}.grm.id") -eq "$N_TRAIN" ]] || { echo "ERROR: GRM sample count mismatch" >&2; exit 1; }
N_EXPECT=$(awk -v chr="$CHR" '$1==chr{n++} END{print n+0}' "${BFILE}.bim")
(( N_EXPECT >= 1000 )) || { echo "ERROR: too few imputed SNPs on chr $CHR: $N_EXPECT" >&2; exit 1; }

rm -f "$OUT/complete.flag" "$OUT/gwas_qc.tsv"
printf 'Outer\tInner\tTrait\tChromosome\tExpected_SNPs\tTested_SNPs\tValid_P\tInvalid_P\n' > "$OUT/gwas_qc.tsv"
for TRAIT in BB H6; do
  PHENO="$PHENO_ROOT/$FOLD/${TRAIT}.train.pheno"
  [[ -s "$PHENO" ]] || { echo "ERROR: missing phenotype $PHENO" >&2; exit 1; }
  PREFIX="$OUT/${TRAIT}_chr${PAD}"
  echo "===== GCTA MLMA LOCO by GRM subtraction: outer_${OUTER}/inner_${INNER} ${TRAIT} chr${CHR} ====="
  "$GCTA_BIN" --mlma --bfile "$BFILE" --autosome-num 12 --chr "$CHR" \
    --grm "$GRM_ALL" --mlma-subtract-grm "$GRM_CHR" \
    --pheno "$PHENO" --thread-num "$CPUS" --out "$PREFIX"
  RESULT="${PREFIX}.mlma"
  [[ -s "$RESULT" ]] || { echo "ERROR: missing GCTA result $RESULT" >&2; exit 1; }
  awk -v outer="$OUTER" -v inner="$INNER" -v trait="$TRAIT" -v chr="$CHR" -v expected="$N_EXPECT" '
    function is_probability(x) {
      if (x !~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ && x !~ /^[.][0-9]+([eE][-+]?[0-9]+)?$/) return 0
      return (x+0 >= 0 && x+0 <= 1)
    }
    NR==1 { if (NF<9 || $1!="Chr" || $2!="SNP" || $9!="p") {print "ERROR: unexpected GCTA MLMA header" > "/dev/stderr"; bad=1; exit 1}; next }
    { if (NF<9 || $1!=chr) {print "ERROR: invalid chromosome/columns" > "/dev/stderr"; bad=1; exit 1}; tested++; if(is_probability($9)) valid++; else invalid++ }
    END { if(bad)exit 1; if(tested<1000 || valid<1000 || tested>expected) {print "ERROR: unexpected GWAS result count" > "/dev/stderr"; exit 1}; printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",outer,inner,trait,chr,expected,tested,valid,invalid }
  ' "$RESULT" >> "$OUT/gwas_qc.tsv"
done
printf 'status=PASS\nouter=%s\ninner=%s\nchr=%s\ntraits=BB,H6\nmodel=GCTA_MLMA_with_mlma_subtract_grm\n' "$OUTER" "$INNER" "$CHR" > "$OUT/complete.flag"
echo "PASS GWAS outer_${OUTER}/inner_${INNER}/chr_${PAD}"
