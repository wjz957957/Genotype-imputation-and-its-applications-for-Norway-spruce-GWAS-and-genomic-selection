#!/bin/bash
#SBATCH -J cf09_imp_prep
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-4%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09c_imp_prep_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09c_imp_prep_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 09c line=${LINENO} exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1

PLINK_BIN="$(command -v plink)"
if command -v gcta64 >/dev/null 2>&1; then
    GCTA_BIN="$(command -v gcta64)"
elif command -v gcta >/dev/null 2>&1; then
    GCTA_BIN="$(command -v gcta)"
else
    echo "ERROR: GCTA executable not found" >&2
    exit 127
fi

MASTER="$WORK/08_imputed_GWAS/00_reference/all_904_genotypes"
INPUT_ROOT="$WORK/Cross_family/results/09_outer_inputs"
CHIP_ROOT="$WORK/Cross_family/results/09_outer_chip_genotypes"
OUT_ROOT="$WORK/Cross_family/results/09_outer_imputed_gwas_inputs"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="$CPUS"
export OPENBLAS_NUM_THREADS=1

TASK="${SLURM_ARRAY_TASK_ID:?Run as Slurm array task 0-4}"
[[ "$TASK" =~ ^[0-4]$ ]] || { echo "ERROR: invalid array task id: $TASK" >&2; exit 2; }
OUTER="$(printf '%02d' "$((TASK + 1))")"
INPUT="$INPUT_ROOT/outer_${OUTER}"
CHIP_PREP="$CHIP_ROOT/outer_${OUTER}"
CHIP_ALL_BFILE="$CHIP_PREP/chip_ALL_outer_train"
OUT="$OUT_ROOT/outer_${OUTER}"

for ext in bed bim fam; do
    [[ -s "${MASTER}.${ext}" ]] || { echo "ERROR: missing master ${MASTER}.${ext}" >&2; exit 1; }
    [[ -s "${CHIP_ALL_BFILE}.${ext}" ]] || { echo "ERROR: missing prepared outer-training chip genotype ${CHIP_ALL_BFILE}.${ext}" >&2; exit 1; }
done
for f in "$INPUT/prepare.complete" "$INPUT/outer_train.keep" "$INPUT/outer_test.ids" \
         "$INPUT/chip_ALL.snps" "$INPUT/BB.train.pheno" "$INPUT/H6.train.pheno" \
         "$INPUT/frozen_inner_selection.tsv" "$CHIP_PREP/complete.flag"; do
    [[ -s "$f" ]] || { echo "ERROR: missing/empty prerequisite: $f" >&2; exit 1; }
done
grep -qx 'status=PASS' "$INPUT/prepare.complete"
grep -qx 'outer_test_phenotype_read=NO' "$INPUT/prepare.complete"
grep -qx 'status=PASS' "$CHIP_PREP/complete.flag"
grep -qx "outer=$OUTER" "$CHIP_PREP/complete.flag"
grep -qx 'training_only=PASS' "$CHIP_PREP/complete.flag"
grep -qx 'outer_test_phenotype_read=NO' "$CHIP_PREP/complete.flag"
grep -qx 'chip_ALL=43267' "$CHIP_PREP/complete.flag"

mkdir -p "$OUT"

# Validate the outer-training / outer-test partition against the 904-sample master.
awk 'NF<2 || $1!=$2 {bad=1} END{exit bad}' "${MASTER}.fam" || {
    echo "ERROR: master FAM must have FID=IID for this workflow" >&2; exit 1;
}
[[ "$(wc -l < "${MASTER}.fam")" -eq 904 ]] || {
    echo "ERROR: expected 904 master individuals" >&2; exit 1;
}
awk 'NF!=2 || $1!=$2 {bad=1} {if(seen[$2]++)bad=1; n++} END{if(bad || n<100)exit 1}' \
    "$INPUT/outer_train.keep" || { echo "ERROR: invalid outer_train.keep" >&2; exit 1; }
awk 'NF!=1 {bad=1} {if(seen[$1]++)bad=1; n++} END{if(bad || n<1)exit 1}' \
    "$INPUT/outer_test.ids" || { echo "ERROR: invalid outer_test.ids" >&2; exit 1; }
for trait in BB H6; do
    pheno="$INPUT/${trait}.train.pheno"
    awk '
      BEGIN{num="^[+-]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eE][+-]?[0-9]+)?$"}
      NF!=3 || $1!=$2 || $3 !~ num || seen[$2]++ {bad=1}
      END{exit bad || NR==0}
    ' "$pheno" || { echo "ERROR: invalid phenotype file $pheno" >&2; exit 1; }
    awk 'BEGIN{OFS="\t"}{print $1,$2}' "$pheno" > "$OUT/.${trait}.samples.tmp"
    cmp -s "$OUT/.${trait}.samples.tmp" "$INPUT/outer_train.keep" || {
        echo "ERROR: ${trait} phenotype sample IDs/order differ from outer_train.keep" >&2; exit 1;
    }
done
N_TRAIN=$(wc -l < "$INPUT/outer_train.keep")
N_TEST=$(wc -l < "$INPUT/outer_test.ids")
[[ "$(wc -l < "$INPUT/BB.train.pheno")" -eq "$N_TRAIN" && \
   "$(wc -l < "$INPUT/H6.train.pheno")" -eq "$N_TRAIN" ]] || {
    echo "ERROR: phenotype/training sample count mismatch" >&2; exit 1;
}

awk 'NR==FNR{train[$2]=1; next} ($1 in train){overlap=1} END{exit overlap}' \
    "$INPUT/outer_train.keep" "$INPUT/outer_test.ids" || {
    echo "ERROR: outer training/test overlap" >&2; exit 1;
}
awk 'NR==FNR{master[$2]=1; next} !($2 in master){bad=1} END{exit bad}' \
    "${MASTER}.fam" "$INPUT/outer_train.keep" || {
    echo "ERROR: outer training ID absent from master FAM" >&2; exit 1;
}
awk 'NR==FNR{master[$2]=1; next} !($1 in master){bad=1} END{exit bad}' \
    "${MASTER}.fam" "$INPUT/outer_test.ids" || {
    echo "ERROR: outer test ID absent from master FAM" >&2; exit 1;
}
awk '{print $2}' "${MASTER}.fam" | sort -u > "$OUT/.master_ids.sorted"
awk '{print $2}' "$INPUT/outer_train.keep" | sort -u > "$OUT/.train_ids.sorted"
sort -u "$INPUT/outer_test.ids" > "$OUT/.test_ids.sorted"
cat "$OUT/.train_ids.sorted" "$OUT/.test_ids.sorted" | sort -u > "$OUT/.partition_ids.sorted"
cmp -s "$OUT/.master_ids.sorted" "$OUT/.partition_ids.sorted" || {
    echo "ERROR: outer train/test do not exactly partition the 904 master individuals" >&2; exit 1;
}

# Validate the frozen full-chip SNP roster and the previously extracted outer-training chip set.
N_CHIP_LIST=$(awk 'NF!=1 || seen[$1]++ {bad=1} {n++} END{if(bad)exit 1; print n+0}' "$INPUT/chip_ALL.snps") || {
    echo "ERROR: malformed or duplicated chip_ALL.snps" >&2; exit 1;
}
[[ "$N_CHIP_LIST" -eq 43267 ]] || {
    echo "ERROR: expected 43267 chip SNP IDs, found $N_CHIP_LIST" >&2; exit 1;
}
awk 'NR==FNR{known[$2]=1; next} !($1 in known){bad=1} END{exit bad}' \
    "${MASTER}.bim" "$INPUT/chip_ALL.snps" || {
    echo "ERROR: chip_ALL contains SNP IDs absent from master BIM" >&2; exit 1;
}
[[ "$(wc -l < "${CHIP_ALL_BFILE}.fam")" -eq "$N_TRAIN" ]] || {
    echo "ERROR: prepared chip genotype sample count mismatch" >&2; exit 1;
}
[[ "$(wc -l < "${CHIP_ALL_BFILE}.bim")" -eq 43267 ]] || {
    echo "ERROR: prepared chip genotype does not contain exactly 43267 SNPs" >&2; exit 1;
}
awk '{print $2}' "${CHIP_ALL_BFILE}.bim" | sort -u > "$OUT/.chip_bim_ids.sorted"
sort -u "$INPUT/chip_ALL.snps" > "$OUT/.chip_list_ids.sorted"
cmp -s "$OUT/.chip_bim_ids.sorted" "$OUT/.chip_list_ids.sorted" || {
    echo "ERROR: prepared chip BIM differs from frozen chip_ALL list" >&2; exit 1;
}
awk 'NF<2 || $1!=$2 {bad=1} {print $1 "\t" $2}' "${CHIP_ALL_BFILE}.fam" | \
    sort -k1,1 -k2,2 > "$OUT/.chip_fam_ids.sorted"
awk '{print $1 "\t" $2}' "$INPUT/outer_train.keep" | \
    sort -k1,1 -k2,2 > "$OUT/.expected_chip_fam_ids.sorted"
cmp -s "$OUT/.chip_fam_ids.sorted" "$OUT/.expected_chip_fam_ids.sorted" || {
    echo "ERROR: prepared chip FAM IDs differ from outer-training keep list" >&2; exit 1;
}

# A successful rerun is a no-op. Otherwise, rebuild incomplete products deterministically.
if [[ -s "$OUT/complete.flag" ]] && grep -qx 'status=PASS' "$OUT/complete.flag" && \
   grep -qx "outer=$OUTER" "$OUT/complete.flag"; then
    for f in "$OUT/imputed_train_qc.bed" "$OUT/imputed_train_qc.bim" "$OUT/imputed_train_qc.fam" \
             "$OUT/chip_train_qc.bed" "$OUT/chip_train_qc.bim" "$OUT/chip_train_qc.fam" \
             "$OUT/chip_grm/chip_all.grm.bin" "$OUT/chip_grm/chip_all.grm.N.bin" "$OUT/chip_grm/chip_all.grm.id"; do
        [[ -s "$f" ]] || { echo "ERROR: completion flag exists but output missing: $f" >&2; exit 1; }
    done
    for CHR in $(seq 1 12); do
        PAD=$(printf '%02d' "$CHR")
        for ext in grm.bin grm.N.bin grm.id; do
            [[ -s "$OUT/chip_grm/chip_chr${PAD}.${ext}" ]] || {
                echo "ERROR: completion flag exists but chromosome GRM missing: chip_chr${PAD}.${ext}" >&2; exit 1;
            }
        done
    done
    rm -f "$OUT"/.*.tmp "$OUT/.master_ids.sorted" "$OUT/.train_ids.sorted" \
          "$OUT/.test_ids.sorted" "$OUT/.partition_ids.sorted" "$OUT/.chip_bim_ids.sorted" \
          "$OUT/.chip_list_ids.sorted" "$OUT/.chip_fam_ids.sorted" "$OUT/.expected_chip_fam_ids.sorted" \
          "$OUT/.observed_samples.sorted" "$OUT/.expected_samples.sorted"
    echo "[SKIP] 09c outer_${OUTER} already complete"
    exit 0
fi
rm -f "$OUT/complete.flag"
IMP_BFILE="$OUT/imputed_train_qc"
CHIP_BFILE="$OUT/chip_train_qc"
GRM_DIR="$OUT/chip_grm"
mkdir -p "$GRM_DIR"
rm -f "${IMP_BFILE}.bed" "${IMP_BFILE}.bim" "${IMP_BFILE}.fam" "${IMP_BFILE}.log" "${IMP_BFILE}.nosex"
rm -f "${CHIP_BFILE}.bed" "${CHIP_BFILE}.bim" "${CHIP_BFILE}.fam" "${CHIP_BFILE}.log" "${CHIP_BFILE}.nosex"
rm -f "$GRM_DIR"/chip_all.grm.* "$GRM_DIR"/chip_chr*.grm.*

# Match inner 06a QC exactly: MAF >= 0.01, SNP missingness <= 0.05, no HWE filter.
echo "[PLINK] outer_${OUTER}: preparing imputed SNPs (exclude all 43267 chip SNPs)"
"$PLINK_BIN" --bfile "$MASTER" --keep "$INPUT/outer_train.keep" \
    --exclude "$INPUT/chip_ALL.snps" --maf 0.01 --geno 0.05 \
    --keep-allele-order --allow-extra-chr --make-bed \
    --threads "$CPUS" --memory 36000 --out "$IMP_BFILE"

echo "[PLINK] outer_${OUTER}: QC chip SNPs for GRM"
"$PLINK_BIN" --bfile "$CHIP_ALL_BFILE" --maf 0.01 --geno 0.05 \
    --keep-allele-order --allow-extra-chr --make-bed \
    --threads "$CPUS" --memory 16000 --out "$CHIP_BFILE"

for prefix in "$IMP_BFILE" "$CHIP_BFILE"; do
    for ext in bed bim fam; do
        [[ -s "${prefix}.${ext}" ]] || { echo "ERROR: missing ${prefix}.${ext}" >&2; exit 1; }
    done
    [[ "$(wc -l < "${prefix}.fam")" -eq "$N_TRAIN" ]] || {
        echo "ERROR: sample count mismatch in $prefix" >&2; exit 1;
    }
    awk '{if(NF<2 || $1!=$2)bad=1; print $1 "\t" $2}' "${prefix}.fam" | sort -k1,1 -k2,2 > "$OUT/.observed_samples.sorted"
    sort -k1,1 -k2,2 "$INPUT/outer_train.keep" > "$OUT/.expected_samples.sorted"
    cmp -s "$OUT/.observed_samples.sorted" "$OUT/.expected_samples.sorted" || {
        echo "ERROR: sample IDs differ from outer-training keep list for $prefix" >&2; exit 1;
    }
done
N_IMP=$(wc -l < "${IMP_BFILE}.bim")
N_CHIP_QC=$(wc -l < "${CHIP_BFILE}.bim")
(( N_IMP >= 100000 )) || { echo "ERROR: unexpectedly few imputed SNPs after QC: $N_IMP" >&2; exit 1; }
(( N_CHIP_QC >= 10000 )) || { echo "ERROR: unexpectedly few chip SNPs after QC: $N_CHIP_QC" >&2; exit 1; }
awk 'NR==FNR{chip[$1]=1; next} ($2 in chip){bad=1} END{exit bad}' \
    "$INPUT/chip_ALL.snps" "${IMP_BFILE}.bim" || {
    echo "ERROR: chip SNP overlap remains in imputed QC set" >&2; exit 1;
}
awk '$1 !~ /^(1|2|3|4|5|6|7|8|9|10|11|12)$/ {bad=1} END{exit bad}' "${IMP_BFILE}.bim" || {
    echo "ERROR: non-autosomal chromosome label in imputed QC set" >&2; exit 1;
}
awk 'BEGIN{OFS="\t"}{n[$1]++} END{for(c=1;c<=12;c++)print c,n[c]+0}' \
    "${IMP_BFILE}.bim" > "$OUT/imputed_chr_counts.tsv"
if awk '$2<1000{bad=1} END{exit !bad}' "$OUT/imputed_chr_counts.tsv"; then
    echo "ERROR: at least one chromosome has fewer than 1000 imputed SNPs" >&2; exit 1
fi

# Full-chip and chromosome-specific GRMs for GCTA MLMA + --mlma-subtract-grm LOCO.
echo "[GCTA] outer_${OUTER}: full-chip GRM"
"$GCTA_BIN" --bfile "$CHIP_BFILE" --autosome-num 12 --make-grm \
    --thread-num "$CPUS" --out "$GRM_DIR/chip_all"
for ext in grm.bin grm.N.bin grm.id; do
    [[ -s "$GRM_DIR/chip_all.${ext}" ]] || { echo "ERROR: missing full-chip GRM: $GRM_DIR/chip_all.${ext}" >&2; exit 1; }
done
for CHR in $(seq 1 12); do
    PAD=$(printf '%02d' "$CHR")
    PREFIX="$GRM_DIR/chip_chr${PAD}"
    N_THIS=$(awk -v c="$CHR" '$1==c{n++} END{print n+0}' "${CHIP_BFILE}.bim")
    (( N_THIS > 0 )) || { echo "ERROR: no QC chip SNPs on chromosome $CHR" >&2; exit 1; }
    echo "[GCTA] outer_${OUTER}: chromosome $CHR GRM ($N_THIS SNPs)"
    "$GCTA_BIN" --bfile "$CHIP_BFILE" --autosome-num 12 --chr "$CHR" \
        --make-grm --thread-num "$CPUS" --out "$PREFIX"
    for ext in grm.bin grm.N.bin grm.id; do
        [[ -s "${PREFIX}.${ext}" ]] || { echo "ERROR: missing chromosome GRM ${PREFIX}.${ext}" >&2; exit 1; }
    done
    cmp -s "$GRM_DIR/chip_all.grm.id" "${PREFIX}.grm.id" || {
        echo "ERROR: full/chr GRM sample IDs/order differ on chromosome $CHR" >&2; exit 1;
    }
done
[[ "$(wc -l < "$GRM_DIR/chip_all.grm.id")" -eq "$N_TRAIN" ]] || {
    echo "ERROR: GRM sample count mismatch" >&2; exit 1;
}

printf 'Item\tValue\n' > "$OUT/input_summary.tsv"
printf 'outer\t%s\n' "$OUTER" >> "$OUT/input_summary.tsv"
printf 'n_training\t%s\n' "$N_TRAIN" >> "$OUT/input_summary.tsv"
printf 'n_outer_test_excluded\t%s\n' "$N_TEST" >> "$OUT/input_summary.tsv"
printf 'chip_markers_requested\t%s\n' "$N_CHIP_LIST" >> "$OUT/input_summary.tsv"
printf 'chip_markers_after_qc\t%s\n' "$N_CHIP_QC" >> "$OUT/input_summary.tsv"
printf 'imputed_markers_after_qc\t%s\n' "$N_IMP" >> "$OUT/input_summary.tsv"
printf 'maf_filter\t0.01\n' >> "$OUT/input_summary.tsv"
printf 'geno_missingness_max\t0.05\n' >> "$OUT/input_summary.tsv"
printf 'hwe_filter\tNONE\n' >> "$OUT/input_summary.tsv"
printf 'grm\tfull_chip_plus_12_chromosome_specific_grms\n' >> "$OUT/input_summary.tsv"
printf 'loco_implementation\tGCTA_MLMA_with_mlma_subtract_grm\n' >> "$OUT/input_summary.tsv"
printf 'training_only\tPASS\n' >> "$OUT/input_summary.tsv"

cat > "$OUT/complete.flag.tmp" <<EOF_FLAG
version=1
outer=$OUTER
n_training=$N_TRAIN
n_outer_test_excluded=$N_TEST
n_chip_markers_requested=$N_CHIP_LIST
n_chip_markers_after_qc=$N_CHIP_QC
n_imputed_markers_after_qc=$N_IMP
maf_filter=0.01
geno_filter=0.05
hwe_filter=NONE
full_chip_grm=PASS
chromosome_chip_grms=12/12
training_only=PASS
outer_test_phenotype_read=NO
status=PASS
EOF_FLAG
mv "$OUT/complete.flag.tmp" "$OUT/complete.flag"
rm -f "$OUT"/.*.tmp "$OUT/.master_ids.sorted" "$OUT/.train_ids.sorted" \
      "$OUT/.test_ids.sorted" "$OUT/.partition_ids.sorted" "$OUT/.chip_bim_ids.sorted" \
      "$OUT/.chip_list_ids.sorted" "$OUT/.chip_fam_ids.sorted" "$OUT/.expected_chip_fam_ids.sorted" \
      "$OUT/.observed_samples.sorted" "$OUT/.expected_samples.sorted"

echo "PASS 09c outer_${OUTER}: train=${N_TRAIN}, imputed_SNPs=${N_IMP}, chip_SNPs_after_QC=${N_CHIP_QC}, GRM=13/13"
