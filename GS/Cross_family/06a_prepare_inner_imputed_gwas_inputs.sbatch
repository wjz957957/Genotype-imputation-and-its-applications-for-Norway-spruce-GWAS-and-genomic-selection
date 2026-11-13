#!/bin/bash
#SBATCH -J cf_qtn_prep
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-14%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06a_qtn_prep_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06a_qtn_prep_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1

PLINK_BIN="$(command -v plink)"
if command -v gcta64 >/dev/null 2>&1; then GCTA_BIN="$(command -v gcta64)"
elif command -v gcta >/dev/null 2>&1; then GCTA_BIN="$(command -v gcta)"
else echo "ERROR: GCTA executable not found" >&2; exit 127; fi

MASTER="$WORK/08_imputed_GWAS/00_reference/all_904_genotypes"
PHENO_ROOT="$WORK/Cross_family/results/02_inner_phenotypes"
CHIP_ROOT="$WORK/Cross_family/results/03_inner_chip_gwas_clumping"
OUT_ROOT="$WORK/Cross_family/results/06_inner_imputed_gwas_inputs"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="$CPUS"
export OPENBLAS_NUM_THREADS=1

for ext in bed bim fam; do [[ -s "${MASTER}.${ext}" ]] || { echo "ERROR: missing current master ${MASTER}.${ext}" >&2; exit 1; }; done

TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
(( TASK >= 0 && TASK < 15 )) || { echo "Invalid task id: $TASK" >&2; exit 2; }
OUTER=$(printf '%02d' "$((TASK / 3 + 1))")
INNER=$(printf '%02d' "$((TASK % 3 + 1))")
FOLD="outer_${OUTER}/inner_${INNER}"
FOLD_PHENO="$PHENO_ROOT/$FOLD"
FOLD_CHIP="$CHIP_ROOT/$FOLD/marker_sets"
OUT="$OUT_ROOT/$FOLD"
mkdir -p "$OUT"
rm -f "$OUT/complete.flag"

BB_PHENO="$FOLD_PHENO/BB.train.pheno"
H6_PHENO="$FOLD_PHENO/H6.train.pheno"
CHIP_BB="$FOLD_CHIP/BB_chip_ALL.snps"
CHIP_H6="$FOLD_CHIP/H6_chip_ALL.snps"
for f in "$BB_PHENO" "$H6_PHENO" "$CHIP_BB" "$CHIP_H6"; do [[ -s "$f" ]] || { echo "ERROR: missing current Cross_family input: $f" >&2; exit 1; }; done
cmp -s "$CHIP_BB" "$CHIP_H6" || { echo "ERROR: BB/H6 chip_ALL SNP lists differ in $FOLD" >&2; exit 1; }
CHIP_LIST="$OUT/chip_ALL.snps"
cp "$CHIP_BB" "$CHIP_LIST"
N_CHIP_LIST=$(awk 'NF!=1{bad=1} {if(seen[$1]++)bad=1; n++} END{if(bad)exit 1; print n+0}' "$CHIP_LIST")
[[ "$N_CHIP_LIST" -eq 43267 ]] || { echo "ERROR: expected 43267 chip SNP IDs, got $N_CHIP_LIST" >&2; exit 1; }

for pheno in "$BB_PHENO" "$H6_PHENO"; do
  awk 'NF!=3 || $1!=$2 || seen[$2]++ {exit 1}' "$pheno" || { echo "ERROR: expected unique headerless FID IID Y with FID=IID: $pheno" >&2; exit 1; }
done
awk 'BEGIN{OFS="\t"}{print $1,$2}' "$BB_PHENO" > "$OUT/inner_train.keep"
awk 'BEGIN{OFS="\t"}{print $1,$2}' "$H6_PHENO" > "$OUT/h6_samples.check"
cmp -s "$OUT/inner_train.keep" "$OUT/h6_samples.check" || { echo "ERROR: BB/H6 training sample IDs/order differ" >&2; exit 1; }
rm -f "$OUT/h6_samples.check"
N_TRAIN=$(wc -l < "$OUT/inner_train.keep")
(( N_TRAIN >= 100 )) || { echo "ERROR: unexpected training sample count $N_TRAIN" >&2; exit 1; }

# Remove incomplete outputs so a rerun is deterministic.
IMP_BFILE="$OUT/imputed_train_qc"
CHIP_BFILE="$OUT/chip_train_qc"
rm -f "${IMP_BFILE}.bed" "${IMP_BFILE}.bim" "${IMP_BFILE}.fam" "${IMP_BFILE}.log" \
      "${CHIP_BFILE}.bed" "${CHIP_BFILE}.bim" "${CHIP_BFILE}.fam" "${CHIP_BFILE}.log"
rm -f "$OUT/chip_grm"/*.grm.* 2>/dev/null || true
mkdir -p "$OUT/chip_grm"

# Original analysis QC: MAF >= 0.01 and SNP missingness <= 0.05.
# Imputed SNPs are defined from the current master by excluding the current
# 43,267-SNP chip_ALL list, rather than reusing old 10-fold marker-list paths.
echo "[PLINK] Current inner-training imputed SNP set: MAF 0.01, GENO 0.05"
"$PLINK_BIN" --bfile "$MASTER" --keep "$OUT/inner_train.keep" --exclude "$CHIP_LIST" \
  --maf 0.01 --geno 0.05 --keep-allele-order --allow-extra-chr --make-bed \
  --threads "$CPUS" --memory 36000 --out "$IMP_BFILE"

echo "[PLINK] Current inner-training full chip SNP set for GRM: MAF 0.01, GENO 0.05"
"$PLINK_BIN" --bfile "$MASTER" --keep "$OUT/inner_train.keep" --extract "$CHIP_LIST" \
  --maf 0.01 --geno 0.05 --keep-allele-order --allow-extra-chr --make-bed \
  --threads "$CPUS" --memory 16000 --out "$CHIP_BFILE"

for prefix in "$IMP_BFILE" "$CHIP_BFILE"; do
  for ext in bed bim fam; do [[ -s "${prefix}.${ext}" ]] || { echo "ERROR: missing ${prefix}.${ext}" >&2; exit 1; }; done
done
N_IMP=$(wc -l < "${IMP_BFILE}.bim")
N_CHIP_QC=$(wc -l < "${CHIP_BFILE}.bim")
[[ $(wc -l < "${IMP_BFILE}.fam") -eq "$N_TRAIN" ]] || { echo "ERROR: imputed genotype sample count mismatch" >&2; exit 1; }
[[ $(wc -l < "${CHIP_BFILE}.fam") -eq "$N_TRAIN" ]] || { echo "ERROR: chip genotype sample count mismatch" >&2; exit 1; }
(( N_IMP >= 100000 )) || { echo "ERROR: unexpectedly few imputed SNPs after QC: $N_IMP" >&2; exit 1; }
(( N_CHIP_QC >= 10000 )) || { echo "ERROR: unexpectedly few chip SNPs after QC: $N_CHIP_QC" >&2; exit 1; }

# Check the current training sample IDs, independent of PLINK's .fam ordering.
for prefix in "$IMP_BFILE" "$CHIP_BFILE"; do
  awk 'BEGIN{OFS="\t"}{if($1!=$2)exit 1; print $1,$2}' "${prefix}.fam" | sort -k1,1 -k2,2 > "$OUT/observed.samples.sorted"
  sort -k1,1 -k2,2 "$OUT/inner_train.keep" > "$OUT/expected.samples.sorted"
  cmp -s "$OUT/observed.samples.sorted" "$OUT/expected.samples.sorted" || { echo "ERROR: sample IDs mismatch for $prefix" >&2; exit 1; }
done

# Ensure the imputed test set excludes all chip IDs and is limited to chr 1-12.
awk 'NR==FNR{chip[$1]=1;next} ($2 in chip){bad=1} END{exit bad}' "$CHIP_LIST" "${IMP_BFILE}.bim" || { echo "ERROR: chip SNP overlap remains in imputed GWAS set" >&2; exit 1; }
awk '$1 !~ /^(1|2|3|4|5|6|7|8|9|10|11|12)$/ {bad=1} END{exit bad}' "${IMP_BFILE}.bim" || { echo "ERROR: non-autosomal chromosome labels in imputed input" >&2; exit 1; }
awk 'BEGIN{OFS="\t"}{n[$1]++} END{for(c=1;c<=12;c++)print c,n[c]+0}' "${IMP_BFILE}.bim" > "$OUT/imputed_chr_counts.tsv"
if awk '$2<1000{bad=1} END{exit !bad}' "$OUT/imputed_chr_counts.tsv"; then echo "ERROR: at least one chromosome has fewer than 1000 imputed SNPs" >&2; exit 1; fi

# Build the full-chip GRM and chromosome-specific chip GRMs for the original
# LOCO implementation: --mlma plus --mlma-subtract-grm.
GRM_DIR="$OUT/chip_grm"
GRM_ALL="$GRM_DIR/chip_all"
"$GCTA_BIN" --bfile "$CHIP_BFILE" --autosome-num 12 --make-grm --thread-num "$CPUS" --out "$GRM_ALL"
for ext in grm.bin grm.N.bin grm.id; do [[ -s "${GRM_ALL}.${ext}" ]] || { echo "ERROR: missing full-chip GRM ${GRM_ALL}.${ext}" >&2; exit 1; }; done
for CHR in $(seq 1 12); do
  PAD=$(printf '%02d' "$CHR")
  PREFIX="$GRM_DIR/chip_chr${PAD}"
  N_THIS=$(awk -v c="$CHR" '$1==c{n++} END{print n+0}' "${CHIP_BFILE}.bim")
  (( N_THIS > 0 )) || { echo "ERROR: no QC chip SNPs on chromosome $CHR" >&2; exit 1; }
  "$GCTA_BIN" --bfile "$CHIP_BFILE" --autosome-num 12 --chr "$CHR" --make-grm --thread-num "$CPUS" --out "$PREFIX"
  for ext in grm.bin grm.N.bin grm.id; do [[ -s "${PREFIX}.${ext}" ]] || { echo "ERROR: missing chromosome GRM ${PREFIX}.${ext}" >&2; exit 1; }; done
  cmp -s "${GRM_ALL}.grm.id" "${PREFIX}.grm.id" || { echo "ERROR: full/chr GRM IDs or order differ on chr $CHR" >&2; exit 1; }
done
[[ $(wc -l < "${GRM_ALL}.grm.id") -eq "$N_TRAIN" ]] || { echo "ERROR: GRM sample count mismatch" >&2; exit 1; }

cat > "$OUT/input_summary.tsv" <<EOF
Item	Value
outer	$OUTER
inner	$INNER
n_training	$N_TRAIN
chip_markers_requested	$N_CHIP_LIST
chip_markers_after_qc	$N_CHIP_QC
imputed_markers_after_qc	$N_IMP
maf_filter	0.01
geno_missingness_max	0.05
hwe_filter	NONE
chip_grm	full_chip_plus_chromosome_specific_grms
loco_implementation	GCTA_MLMA_with_mlma_subtract_grm
EOF
cat > "$OUT/complete.flag" <<EOF
version=2
outer=$OUTER
inner=$INNER
n_training=$N_TRAIN
n_chip_markers_requested=$N_CHIP_LIST
n_chip_markers_after_qc=$N_CHIP_QC
n_imputed_markers_after_qc=$N_IMP
maf_filter=0.01
geno_filter=0.05
hwe_filter=NONE
full_chip_grm=PASS
chromosome_chip_grms=12/12
training_only=PASS
status=PASS
EOF

echo "PASS: prepared outer_${OUTER}/inner_${INNER}; n_train=$N_TRAIN; chip_QC=$N_CHIP_QC; imputed_QC=$N_IMP"
