#!/bin/bash
#SBATCH -J cf09g_xchipLD
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-9%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09g_xchipLD_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09g_xchipLD_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 09g line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9

PLINK_BIN="$(command -v plink)"
PYTHON_BIN="$(command -v python3)"
FILTER="$WORK/Cross_family/scripts/09g_filter_outer_crosschip.py"
MASTER="$WORK/08_imputed_GWAS/00_reference/all_904_genotypes"
OUTER_INPUT_ROOT="$WORK/Cross_family/results/09_outer_inputs"
IMPUTED_INPUT_ROOT="$WORK/Cross_family/results/09_outer_imputed_gwas_inputs"
CHIP_ROOT="$WORK/Cross_family/results/09_outer_chip_gwas_clumping"
MERGED_ROOT="$WORK/Cross_family/results/09_outer_imputed_clumping_merged"
MERGE_SUMMARY="$MERGED_ROOT/merge_summary.tsv"
MERGE_FLAG="$MERGED_ROOT/complete.flag"
OUT_ROOT="$WORK/Cross_family/results/09_outer_crosschip_ld_filtering"
R2=0.8
CPUS="${SLURM_CPUS_PER_TASK:-4}"

for F in "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam" "$FILTER" \
         "$MERGE_SUMMARY" "$MERGE_FLAG"; do
    [[ -s "$F" ]] || { echo "ERROR: required input missing/empty: $F" >&2; exit 1; }
done
grep -qx 'status=PASS' "$MERGE_FLAG"
grep -qx 'outer_folds=5/5' "$MERGE_FLAG"
grep -qx 'trait_combinations=10/10' "$MERGE_FLAG"
grep -qx 'training_only=PASS' "$MERGE_FLAG"
grep -qx 'outer_test_phenotype_read=NO' "$MERGE_FLAG"

TASK="${SLURM_ARRAY_TASK_ID:?Run as Slurm array 0-9}"
(( TASK >= 0 && TASK < 10 )) || { echo "Invalid task id: $TASK" >&2; exit 2; }
OUTER="$(printf '%02d' "$((TASK / 2 + 1))")"
if (( TASK % 2 == 0 )); then TRAIT=BB; else TRAIT=H6; fi

FOLD_INPUT="$OUTER_INPUT_ROOT/outer_${OUTER}"
IMPUTED_INPUT="$IMPUTED_INPUT_ROOT/outer_${OUTER}"
FROZEN="$FOLD_INPUT/frozen_inner_selection.tsv"
PREP_FLAG="$FOLD_INPUT/prepare.complete"
IMPUTED_FLAG="$IMPUTED_INPUT/complete.flag"
TRAIN_KEEP="$FOLD_INPUT/outer_train.keep"
TEST_IDS="$FOLD_INPUT/outer_test.ids"
RANKED="$MERGED_ROOT/outer_${OUTER}/${TRAIT}_index_ranked.tsv"
FOLD_MERGE_FLAG="$MERGED_ROOT/outer_${OUTER}/complete.flag"
CHIP_FOLD="$CHIP_ROOT/outer_${OUTER}"
OUT="$OUT_ROOT/outer_${OUTER}/$TRAIT"

for F in "$FROZEN" "$PREP_FLAG" "$IMPUTED_FLAG" "$TRAIN_KEEP" "$TEST_IDS" \
         "$RANKED" "$FOLD_MERGE_FLAG" "$CHIP_FOLD/complete.flag"; do
    [[ -s "$F" ]] || { echo "ERROR: required input missing/empty: $F" >&2; exit 1; }
done
for ENTRY in 'status=PASS' "outer=$OUTER" 'outer_test_phenotype_read=NO' 'full_chip_consistency=PASS'; do
    grep -qx "$ENTRY" "$PREP_FLAG" || { echo "ERROR: 09a flag mismatch: $ENTRY" >&2; exit 1; }
done
for ENTRY in "outer=$OUTER" 'status=PASS' 'training_only=PASS' 'outer_test_phenotype_read=NO'; do
    grep -qx "$ENTRY" "$IMPUTED_FLAG" || { echo "ERROR: 09c flag mismatch: $ENTRY" >&2; exit 1; }
done
for ENTRY in "status=PASS" "outer=$OUTER" 'chromosomes=12/12' 'training_only=PASS' 'outer_test_phenotype_read=NO'; do
    grep -qx "$ENTRY" "$FOLD_MERGE_FLAG" || { echo "ERROR: 09f flag mismatch: $ENTRY" >&2; exit 1; }
done

read -r BEST_SET TARGET < <(awk -F '\t' -v o="$OUTER" -v t="$TRAIT" '
    NR==1 {if ($1!="Outer" || $2!="Trait" || $3!="Best_Chip" || $4!="Best_QTN") exit 2; next}
    $1==o && $2==t {if (++n>1) exit 3; chip=$3; qtn=$4}
    END {if (n!=1) exit 4; print chip, qtn}
' "$FROZEN")
[[ "$BEST_SET" =~ ^chip_(50|100|200|300|400|500|1000|2000|3000|4000|5000|8000|10000|20000|30000|40000)$ ]] || {
    echo "ERROR: invalid frozen chip strategy $BEST_SET for outer_${OUTER}/${TRAIT}" >&2; exit 1;
}
[[ "$TARGET" =~ ^(0|1|2|3|4|5|10|20|30|40|50|75|100|150|200|300|400|500|600|700)$ ]] || {
    echo "ERROR: invalid frozen QTN count $TARGET for outer_${OUTER}/${TRAIT}" >&2; exit 1;
}
CHIP="$CHIP_FOLD/marker_sets/${TRAIT}_${BEST_SET}.snps"
[[ -s "$CHIP" ]] || { echo "ERROR: selected chip list missing: $CHIP" >&2; exit 1; }
N_CHIP="$(wc -l < "$CHIP")"
EXPECTED_CHIP="${BEST_SET#chip_}"
[[ "$N_CHIP" -eq "$EXPECTED_CHIP" ]] || {
    echo "ERROR: selected chip list has $N_CHIP rows, expected $EXPECTED_CHIP" >&2; exit 1;
}

N_SUMMARY="$(awk -F '\t' -v o="$OUTER" -v t="$TRAIT" '
    NR==1 {if ($1!="outer" || $2!="trait" || $3!="n_candidate_snps") exit 2; next}
    $1==o && $2==t {n++; value=$3}
    END {if (n!=1 || value !~ /^[0-9]+$/) exit 3; print value}
' "$MERGE_SUMMARY")" || { echo "ERROR: invalid 09f merge_summary row for outer_${OUTER}/${TRAIT}" >&2; exit 1; }
[[ "$(head -n 1 "$RANKED")" == $'SNP\tP' ]] || { echo "ERROR: invalid ranked header: $RANKED" >&2; exit 1; }
N_RANKED_FILE=$(( $(wc -l < "$RANKED") - 1 ))
[[ "$N_RANKED_FILE" -eq "$N_SUMMARY" ]] || {
    echo "ERROR: 09f summary says $N_SUMMARY SNPs but $RANKED contains $N_RANKED_FILE rows" >&2; exit 1;
}

mkdir -p "$OUT"
SIG="$OUT/input_signature.txt"
SIG_NOW="$OUT/input_signature.now"
{
    sha256sum "$RANKED" "$CHIP" "$FROZEN" "$TRAIN_KEEP" "$TEST_IDS" \
        "$MERGE_SUMMARY" "$MERGE_FLAG" "$PREP_FLAG" "$IMPUTED_FLAG" "$FOLD_MERGE_FLAG" "$CHIP_FOLD/complete.flag"
    stat -c 'file=%n size=%s mtime=%Y' "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam"
} > "$SIG_NOW"
if [[ -s "$SIG" ]]; then
    cmp -s "$SIG" "$SIG_NOW" || {
        rm -f "$SIG_NOW"
        echo "ERROR: inputs changed since 09g started for outer_${OUTER}/${TRAIT}; inspect outputs before rerun." >&2
        exit 1
    }
    rm -f "$SIG_NOW"
else
    if find "$OUT" -mindepth 1 -maxdepth 1 -type f ! -name 'input_signature.now' | grep -q . || \
       find "$OUT" -mindepth 1 -maxdepth 1 -type d | grep -q .; then
        rm -f "$SIG_NOW"
        echo "ERROR: output directory has files but no input_signature.txt: $OUT" >&2
        exit 1
    fi
    mv "$SIG_NOW" "$SIG"
fi

TRAIN_IDS="$OUT/outer_train.ids"
awk 'NF!=2 || $1!=$2 {exit 1} {print $2}' "$TRAIN_KEEP" > "$TRAIN_IDS"
[[ -s "$TRAIN_IDS" ]] || { echo "ERROR: could not derive outer_train IDs" >&2; exit 1; }
N_TRAIN="$(wc -l < "$TRAIN_IDS")"
N_TEST="$(wc -l < "$TEST_IDS")"
(( N_TRAIN >= 500 && N_TEST >= 50 && N_TRAIN + N_TEST == 904 )) || {
    echo "ERROR: unexpected outer partition sizes train=$N_TRAIN test=$N_TEST" >&2; exit 1;
}

if [[ -s "$OUT/complete.flag" ]]; then
    if grep -qx 'status=PASS' "$OUT/complete.flag" && \
       grep -qx "outer=$OUTER" "$OUT/complete.flag" && \
       grep -qx "trait=$TRAIT" "$OUT/complete.flag" && \
       grep -qx "chip_strategy=$BEST_SET" "$OUT/complete.flag" && \
       grep -qx "selected_qtn=$TARGET" "$OUT/complete.flag" && \
       [[ -f "$OUT/marker_sets/extra_${TARGET}.snps" ]] && \
       [[ "$(wc -l < "$OUT/marker_sets/extra_${TARGET}.snps")" -eq "$TARGET" ]]; then
        echo "[SKIP] outer_${OUTER}/${TRAIT} already complete; chip=$BEST_SET, QTN=$TARGET"
        exit 0
    fi
    echo "ERROR: existing completion flag or selected marker set conflicts with current inputs: $OUT" >&2
    exit 1
fi

# Reuse the selected chip and candidate rank, but never use outer-test phenotypes.
"$PYTHON_BIN" "$FILTER" rank --ranked "$RANKED" --chip "$CHIP" \
    --out "$OUT/ranked_nonchip.tsv" --summary "$OUT/rank_nonchip_summary.tsv"
N_NONCHIP=$(( $(wc -l < "$OUT/ranked_nonchip.tsv") - 1 ))
(( N_NONCHIP > 0 )) || { echo "ERROR: no non-chip candidates remain" >&2; exit 1; }
echo "===== outer_${OUTER}/${TRAIT}: chip=$BEST_SET ($N_CHIP), frozen_QTN=$TARGET, index_SNPs=$N_RANKED_FILE, nonchip=$N_NONCHIP, r2=$R2 ====="

if (( TARGET == 0 )); then
    mkdir -p "$OUT/marker_sets"
    : > "$OUT/marker_sets/extra_0.snps"
    : > "$OUT/selected_top_QTN.snps"
    printf '0\n' > "$OUT/best_extra.txt"
    printf 'SNP\tP\n' > "$OUT/filtered_ranked_extra_snps.tsv"
    printf 'Original_Rank\tSNP\tP\tChr\tLD_status\tMax_reported_high_r2\tTop_linked_CHIP\tKeep\tFiltered_Rank\n' > "$OUT/candidate_ld_annotated.tsv"
    printf 'N_QTN\tAvailable\tReason\n0\t1\tFROZEN_TARGET_ZERO_NO_LD_REQUIRED\n' > "$OUT/candidate_count_availability.tsv"
    printf 'Item\tValue\nFrozen_Best_QTN\t0\nCHIP_count\t%s\nOuter_train_size\t%s\nLD_threshold\t%s\nLD_filter_performed\tNO_QTN_REQUIRED\nRetained\t0\nBest_Extra\t0\n' \
        "$N_CHIP" "$N_TRAIN" "$R2" > "$OUT/filter_summary.tsv"
    cat > "$OUT/complete.flag.tmp" <<EOF_FLAG
status=PASS
version=1
outer=$OUTER
trait=$TRAIT
chip_strategy=$BEST_SET
chip_markers=$N_CHIP
selected_qtn=0
ld_threshold=$R2
ld_filter_performed=NO_QTN_REQUIRED
n_training=$N_TRAIN
n_outer_test=$N_TEST
n_clump_indices=$N_RANKED_FILE
n_nonchip_candidates=$N_NONCHIP
n_retained_after_ld=0
training_only=PASS
outer_test_phenotype_read=NO
EOF_FLAG
    mv "$OUT/complete.flag.tmp" "$OUT/complete.flag"
    echo "PASS outer_${OUTER}/${TRAIT}: frozen QTN count is zero; wrote empty extra_0.snps without LD computation"
    exit 0
fi

POOL="$N_NONCHIP"
STAGE="$OUT/stages/top${POOL}"
mkdir -p "$STAGE"
"$PYTHON_BIN" "$FILTER" prepare --ranked "$OUT/ranked_nonchip.tsv" --chip "$CHIP" \
    --train "$TRAIN_IDS" --test "$TEST_IDS" --pool "$POOL" --stage "$STAGE"

LD_BFILE="$STAGE/train_chip_plus_candidates"
if [[ ! -s "$STAGE/genotype_extract.complete" || ! -s "${LD_BFILE}.bed" || \
      ! -s "${LD_BFILE}.bim" || ! -s "${LD_BFILE}.fam" ]]; then
    rm -f "$STAGE/genotype_extract.complete"
    "$PLINK_BIN" --bfile "$MASTER" --keep "$STAGE/train.keep" \
        --extract "$STAGE/union.snps" --keep-allele-order --make-bed \
        --threads "$CPUS" --memory 36000 --out "$LD_BFILE"
    for EXT in bed bim fam; do
        [[ -s "${LD_BFILE}.${EXT}" ]] || { echo "ERROR: LD genotype output missing ${LD_BFILE}.${EXT}" >&2; exit 1; }
    done
    printf 'DONE outer-training genotype extraction\n' > "$STAGE/genotype_extract.complete"
fi

LD_ARGS=()
for CHR in $(seq 1 12); do
    PAD="$(printf '%02d' "$CHR")"
    CHR_DIR="$STAGE/chr_${PAD}"
    mkdir -p "$CHR_DIR"
    CHR_CHIP="$CHR_DIR/selected_chip.snps"
    awk -v chr="$CHR" 'NR==FNR {selected[$1]=1; next} $1==chr && ($2 in selected) {print $2}' \
        "$STAGE/selected_chip.snps" "${LD_BFILE}.bim" > "$CHR_CHIP"
    N_THIS_CHR="$(wc -l < "$CHR_CHIP")"
    if (( N_THIS_CHR == 0 )); then
        continue
    fi
    PREFIX="$CHR_DIR/chip_qtn_sparse_LD"
    LD_GZ="${PREFIX}.ld.gz"
    if [[ -s "${PREFIX}.complete" && -s "$LD_GZ" ]]; then
        gzip -t "$LD_GZ"
    else
        echo "[LD] chr${CHR}: selected_chip_markers=$N_THIS_CHR, r2 >= $R2"
        "$PLINK_BIN" --bfile "$LD_BFILE" --chr "$CHR" --ld-snp-list "$CHR_CHIP" \
            --r2 gz --ld-window 999999 --ld-window-kb 10000000 \
            --ld-window-r2 "$R2" --threads "$CPUS" --memory 36000 --out "$PREFIX"
        [[ -s "$LD_GZ" ]] || { echo "ERROR: sparse LD file missing: $LD_GZ" >&2; exit 1; }
        gzip -t "$LD_GZ"
        printf 'DONE sparse LD for chromosome %s\n' "$CHR" > "${PREFIX}.complete"
    fi
    LD_ARGS+=(--ld "$LD_GZ")
done
(( ${#LD_ARGS[@]} > 0 )) || { echo "ERROR: no selected chip SNP found on chr1-12" >&2; exit 1; }

"$PYTHON_BIN" "$FILTER" report --ranked "$OUT/ranked_nonchip.tsv" --stage "$STAGE" \
    --bim "${LD_BFILE}.bim" --fam "${LD_BFILE}.fam" --train "$TRAIN_IDS" \
    --pool "$POOL" --threshold "$R2" --target "$TARGET" "${LD_ARGS[@]}"
N_RETAIN="$(awk -F '\t' '$1=="Retained" {print $2}' "$STAGE/filter_summary.tsv")"
[[ "$N_RETAIN" =~ ^[0-9]+$ ]] || { echo "ERROR: unable to parse retained count" >&2; exit 1; }
if (( N_RETAIN < TARGET )); then
    printf 'Pool_tested\t%s\nRetained\t%s\nFrozen_Best_QTN\t%s\nR2_threshold\t%s\n' \
        "$POOL" "$N_RETAIN" "$TARGET" "$R2" > "$OUT/insufficient_candidates.tsv"
    echo "ERROR: insufficient candidates after cross-chip LD: retained=$N_RETAIN < frozen target=$TARGET" >&2
    exit 3
fi

"$PYTHON_BIN" "$FILTER" publish --stage "$STAGE" --out "$OUT" --target "$TARGET"
N_RETAINED_FILE=$(( $(wc -l < "$OUT/filtered_ranked_extra_snps.tsv") - 1 ))
N_SELECTED="$(wc -l < "$OUT/marker_sets/extra_${TARGET}.snps")"
[[ "$N_SELECTED" -eq "$TARGET" && "$N_RETAINED_FILE" -eq "$N_RETAIN" ]] || {
    echo "ERROR: published marker/ranking count mismatch" >&2; exit 1;
}
cat > "$OUT/complete.flag.tmp" <<EOF_FLAG
status=PASS
version=1
outer=$OUTER
trait=$TRAIT
chip_strategy=$BEST_SET
chip_markers=$N_CHIP
selected_qtn=$TARGET
ld_threshold=$R2
ld_filter_performed=YES
n_training=$N_TRAIN
n_outer_test=$N_TEST
n_clump_indices=$N_RANKED_FILE
n_nonchip_candidates=$N_NONCHIP
n_retained_after_ld=$N_RETAIN
training_only=PASS
outer_test_phenotype_read=NO
EOF_FLAG
mv "$OUT/complete.flag.tmp" "$OUT/complete.flag"
echo "PASS outer_${OUTER}/${TRAIT}: chip=$BEST_SET, target_QTN=$TARGET, retained=$N_RETAIN, selected=$N_SELECTED"
