#!/bin/bash
#SBATCH -J 17a_outerLD
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-19%2
#SBATCH -o logs/17a_outerLD_%A_%a.out
#SBATCH -e logs/17a_outerLD_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
TASK=${SLURM_ARRAY_TASK_ID:?Run as SLURM array 0-19}
(( TASK >= 0 && TASK < 20 )) || { echo 'Invalid task ID' >&2; exit 1; }
OUTER=$(printf '%02d' "$((TASK / 2 + 1))")
if (( TASK % 2 == 0 )); then TRAIT=BB; else TRAIT=H6; fi

R2=0.8
PY="${WORK}/scripts/17_filter_outer_crosschip.py"
MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
CV="${WORK}/02_CV/outer_${OUTER}"
TRAIN="${CV}/outer_train.samples"
TEST="${CV}/outer_test.samples"
RANKED="${WORK}/10_outer_test/04_imputed_clumping/outer_${OUTER}/${TRAIT}/ranked_extra_snps.tsv"
CHIP_DIR="${WORK}/10_outer_test/01_chip_selection/outer_${OUTER}/selected"
CHIP="${CHIP_DIR}/${TRAIT}_best_chip.snps"
CHIP_SETTING="${CHIP_DIR}/${TRAIT}_best_marker.txt"
OPT="${WORK}/09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_allfolds/optimization/best_extra_settings.tsv"
OUT="${WORK}/10_outer_test/06_crosschip_LD_r2_0p8_allfolds/outer_${OUTER}/${TRAIT}"

for F in "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam" "$PY" "$TRAIN" "$TEST" "$RANKED" "$CHIP" "$CHIP_SETTING" "$OPT"   "${WORK}/10_outer_test/04_imputed_clumping/outer_${OUTER}/${TRAIT}/complete.flag"; do
    [[ -s "$F" ]] || { echo "ERROR: required file missing/empty: $F" >&2; exit 1; }
done

read -r OPT_CHIP TARGET < <(awk -F '\t' -v o="$OUTER" -v t="$TRAIT" '
  NR==1 { if ($1!="Outer" || $2!="Trait" || $3!="Best_Chip" || $4!="Best_Extra") exit 1; next }
  $1==o && $2==t { if (++n>1) exit 2; print $3, $4 }
  END { if (n!=1) exit 3 }
' "$OPT")
[[ "$TARGET" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid Best_Extra=$TARGET" >&2; exit 1; }
case "$TARGET" in 0|1|2|3|4|5|10|20|30|40|50|75|100|150|200|300|400|500|600|700) ;; *) echo "ERROR: unsupported count=$TARGET" >&2; exit 1;; esac

ACTUAL_CHIP_SETTING=$(tr -d '[:space:]' < "$CHIP_SETTING")
[[ "$OPT_CHIP" == "$ACTUAL_CHIP_SETTING" ]] || { echo "ERROR: CHIP optimization mismatch: $OPT_CHIP != $ACTUAL_CHIP_SETTING" >&2; exit 1; }
N_CHIP=$(wc -l < "$CHIP")
(( N_CHIP > 0 )) || exit 1
N_RANKED=$(( $(wc -l < "$RANKED") - 1 ))
(( N_RANKED >= 1 )) || { echo 'ERROR: no clumped candidates' >&2; exit 1; }

mkdir -p "$OUT"
if [[ -s "$OUT/complete.flag" ]]; then
    [[ -f "$OUT/marker_sets/extra_${TARGET}.snps" && -s "$OUT/best_extra.txt" ]] || { echo 'ERROR: incomplete prior output' >&2; exit 1; }
    PRIOR=$(tr -d '[:space:]' < "$OUT/best_extra.txt")
    COUNT=$(wc -l < "$OUT/marker_sets/extra_${TARGET}.snps")
    [[ "$PRIOR" == "$TARGET" && "$COUNT" -eq "$TARGET" ]] || { echo 'ERROR: prior result conflicts with current best settings' >&2; exit 1; }
    echo "[SKIP] Valid completed outer_${OUTER} $TRAIT selected=$TARGET"
    exit 0
fi

echo "===== Step 17a: outer_${OUTER} ${TRAIT} | CHIP=${OPT_CHIP} N=${N_CHIP} | Best_Extra=${TARGET} ====="

# If inner CV chooses zero additional QTN, build the empty marker file without
# using outer_test, and keep a clear separate result manifest.
if (( TARGET == 0 )); then
    mkdir -p "$OUT/marker_sets"
    : > "$OUT/marker_sets/extra_0.snps"
    printf '0\n' > "$OUT/best_extra.txt"
    printf 'Item\tValue\nBest_Extra\t0\nLD_threshold\t0.8\n' > "$OUT/filter_summary.tsv"
    printf 'DONE zero QTN\n' > "$OUT/complete.flag"
    exit 0
fi

# Candidate pool expands deterministically until enough filtered markers remain.
POOL_CHOICES=(10000 20000 40000 80000 100000)
LAST_POOL=0
for REQUESTED in "${POOL_CHOICES[@]}"; do
    POOL=$REQUESTED
    (( POOL > N_RANKED )) && POOL=$N_RANKED
    (( POOL > LAST_POOL )) || continue
    LAST_POOL=$POOL
    STAGE="$OUT/stages/top${POOL}"
    mkdir -p "$STAGE"
    echo "[POOL] ${POOL} / ${N_RANKED} available clump indices"

    python3 "$PY" prepare --ranked "$RANKED" --chip "$CHIP"       --train "$TRAIN" --test "$TEST" --pool "$POOL" --stage "$STAGE"

    LD_BFILE="$STAGE/train_chip_plus_candidates"
    # Explicitly use the outer_train samples only. No second --maf/--geno filters.
    if [[ ! -s "$STAGE/genotype_extract.complete" ]]; then
        plink --bfile "$MASTER" --keep "$STAGE/train.keep"           --extract "$STAGE/union.snps" --keep-allele-order --make-bed           --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 36000 --out "$LD_BFILE"
        for EXT in bed bim fam; do [[ -s "${LD_BFILE}.${EXT}" ]] || exit 1; done
        printf 'DONE outer_train genotype extraction\n' > "$STAGE/genotype_extract.complete"
    else
        echo "[REUSE] $LD_BFILE"
    fi

    LD_ARGS=()
    for CHR in $(seq 1 12); do
        PAD=$(printf '%02d' "$CHR")
        CHR_DIR="$STAGE/chr_${PAD}"
        mkdir -p "$CHR_DIR"
        CHR_CHIP="$CHR_DIR/selected_chip.snps"
        awk -v chr="$CHR" 'NR==FNR {selected[$1]=1; next} $1==chr && ($2 in selected) {print $2}'           "$STAGE/selected_chip.snps" "${LD_BFILE}.bim" > "$CHR_CHIP"
        N_THIS_CHR=$(wc -l < "$CHR_CHIP")
        if (( N_THIS_CHR == 0 )); then continue; fi
        PREFIX="$CHR_DIR/chip_qtn_sparse_LD"
        LD_GZ="${PREFIX}.ld.gz"
        if [[ -s "${PREFIX}.complete" && -s "$LD_GZ" ]]; then
            gzip -t "$LD_GZ"
        else
            echo "[LD] chromosome $CHR: $N_THIS_CHR CHIP markers"
            plink --bfile "$LD_BFILE" --chr "$CHR" --ld-snp-list "$CHR_CHIP"               --r2 gz --ld-window 999999 --ld-window-kb 10000000               --ld-window-r2 "$R2" --threads "${SLURM_CPUS_PER_TASK:-4}"               --memory 36000 --out "$PREFIX"
            [[ -s "$LD_GZ" ]] || { echo "ERROR: missing sparse LD $LD_GZ" >&2; exit 1; }
            gzip -t "$LD_GZ"
            printf 'DONE sparse LD for chromosome %s\n' "$CHR" > "${PREFIX}.complete"
        fi
        LD_ARGS+=(--ld "$LD_GZ")
    done

    (( ${#LD_ARGS[@]} > 0 )) || { echo 'ERROR: no selected CHIP SNP on chromosomes 1-12' >&2; exit 1; }

    python3 "$PY" report --ranked "$RANKED" --stage "$STAGE"       --bim "${LD_BFILE}.bim" --fam "${LD_BFILE}.fam" --train "$TRAIN"       --pool "$POOL" --threshold "$R2" --target "$TARGET" "${LD_ARGS[@]}"

    N_RETAIN=$(awk -F '\t' '$1=="Retained" {print $2}' "$STAGE/filter_summary.tsv")
    [[ "$N_RETAIN" =~ ^[0-9]+$ ]] || { echo 'ERROR: unable to parse retained count' >&2; exit 1; }
    echo "[RESULT] pool=$POOL retained=$N_RETAIN best_extra=$TARGET"
    if (( N_RETAIN >= TARGET )); then
        python3 "$PY" publish --stage "$STAGE" --out "$OUT" --target "$TARGET"
        echo "[DONE] Step17a outer_${OUTER} ${TRAIT}"
        exit 0
    fi
    (( POOL < N_RANKED )) || break
done
printf 'Pool_tested\t%s\nRetained\t%s\nBest_Extra\t%s\nR2_threshold\t%s\n'   "$LAST_POOL" "$N_RETAIN" "$TARGET" "$R2" > "$OUT/insufficient_candidates.tsv"
echo "ERROR: insufficient low-LD QTN for outer_${OUTER} ${TRAIT}: $N_RETAIN < $TARGET" >&2
exit 3
