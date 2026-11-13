#!/bin/bash
#SBATCH -J cf_qtn_ld
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-29%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/07a_qtn_ld_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/07a_qtn_ld_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9

PLINK_BIN="$(command -v plink)"
PYTHON_BIN="$(command -v python3)"
FILTER="$WORK/Cross_family/scripts/07a_filter_inner_crosschip.py"
CONFIG="$WORK/Cross_family/config/best_chip_by_outer_trait.tsv"
MASTER="$WORK/08_imputed_GWAS/00_reference/all_904_genotypes"
PHENO_ROOT="$WORK/Cross_family/results/02_inner_phenotypes"
CHIP_ROOT="$WORK/Cross_family/results/03_inner_chip_gwas_clumping"
MERGED_ROOT="$WORK/Cross_family/results/06_inner_imputed_clumping_merged"
MERGE_SUMMARY="$MERGED_ROOT/merge_summary.tsv"
OUT_ROOT="$WORK/Cross_family/results/07_inner_crosschip_ld_filtering"
R2=0.8
CPUS="${SLURM_CPUS_PER_TASK:-4}"
POOL_CHOICES=(10000 20000 40000 80000 100000)

for f in "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam" "$FILTER" "$CONFIG"; do
    [[ -s "$f" ]] || { echo "ERROR: missing/empty input: $f" >&2; exit 1; }
done

TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
(( TASK >= 0 && TASK < 30 )) || { echo "Invalid task id: $TASK" >&2; exit 2; }
OUTER=$(printf '%02d' "$((TASK / 6 + 1))")
REMAINDER=$((TASK % 6))
INNER=$(printf '%02d' "$((REMAINDER / 2 + 1))")
if (( REMAINDER % 2 == 0 )); then TRAIT=BB; else TRAIT=H6; fi

BEST_SET=$(awk -F '\t' -v o="$OUTER" -v t="$TRAIT" '
  NR==1 {if ($1!="Outer" || $2!="Trait" || $3!="Best_Chip") exit 2; next}
  $1==o && $2==t {if (++n>1) exit 3; value=$3}
  END {if (n!=1) exit 4; print value}
' "$CONFIG")
[[ "$BEST_SET" =~ ^chip_(50|100|200|300|400|500|1000|2000|3000|4000|5000|8000|10000|20000|30000|40000)$ ]] || {
    echo "ERROR: invalid selected chip strategy '$BEST_SET' for outer_${OUTER}/${TRAIT}" >&2; exit 1;
}

MERGED_FOLD="$MERGED_ROOT/outer_${OUTER}/inner_${INNER}"
RANKED="$MERGED_FOLD/${TRAIT}_index_ranked.tsv"
PHENO="$PHENO_ROOT/outer_${OUTER}/inner_${INNER}/${TRAIT}.train.pheno"
CHIP="$CHIP_ROOT/outer_${OUTER}/inner_${INNER}/marker_sets/${TRAIT}_${BEST_SET}.snps"
OUT="$OUT_ROOT/outer_${OUTER}/inner_${INNER}/$TRAIT"

for f in "$MERGE_SUMMARY" "$RANKED" "$PHENO" "$CHIP"; do
    [[ -s "$f" ]] || { echo "ERROR: required input missing/empty: $f" >&2; exit 1; }
done
# The actual 06d merge layout has one summary table and per-trait *_index_ranked.tsv files;
# it does not create per-trait complete.flag files.
N_SUMMARY=$(awk -F '\t' -v o="$OUTER" -v i="$INNER" -v t="$TRAIT" '
    NR==1 {if (tolower($1)!="outer" || tolower($2)!="inner" || tolower($3)!="trait" || tolower($4)!="n_candidate_snps") exit 2; next}
    $1==o && $2==i && $3==t {n++; count=$4}
    END {if (n!=1 || count !~ /^[0-9]+$/) exit 3; print count}
' "$MERGE_SUMMARY") || { echo "ERROR: expected exactly one valid merge_summary.tsv row for outer_${OUTER}/inner_${INNER}/${TRAIT}" >&2; exit 1; }
RANK_HEADER=$(head -n 1 "$RANKED")
[[ "$RANK_HEADER" == $'SNP\tP' ]] || { echo "ERROR: unexpected ranked file header in $RANKED: $RANK_HEADER" >&2; exit 1; }
N_RANKED_FILE=$(( $(wc -l < "$RANKED") - 1 ))
[[ "$N_RANKED_FILE" -eq "$N_SUMMARY" ]] || {
    echo "ERROR: merge summary says $N_SUMMARY candidates but $RANKED contains $N_RANKED_FILE rows" >&2; exit 1;
}

if [[ -s "$OUT/complete.flag" ]]; then
    grep -qx "outer=$OUTER" "$OUT/complete.flag" &&
    grep -qx "inner=$INNER" "$OUT/complete.flag" &&
    grep -qx "trait=$TRAIT" "$OUT/complete.flag" &&
    grep -qx "chip_strategy=$BEST_SET" "$OUT/complete.flag" || {
        echo "ERROR: existing output flag does not match this task: $OUT/complete.flag" >&2; exit 1;
    }
    echo "[SKIP] Already completed outer_${OUTER}/inner_${INNER}/${TRAIT} with $BEST_SET"
    exit 0
fi
mkdir -p "$OUT"
rm -f "$OUT/complete.flag"

N_CHIP=$(wc -l < "$CHIP")
EXPECTED_CHIP=${BEST_SET#chip_}
[[ "$N_CHIP" -eq "$EXPECTED_CHIP" ]] || {
    echo "ERROR: selected chip list count $N_CHIP != $BEST_SET expected $EXPECTED_CHIP" >&2; exit 1;
}

RANKED_NONCHIP="$OUT/ranked_nonchip.tsv"
"$PYTHON_BIN" "$FILTER" rank --ranked "$RANKED" --chip "$CHIP" \
    --out "$RANKED_NONCHIP" --summary "$OUT/rank_nonchip_summary.tsv"
N_RANKED=$(( $(wc -l < "$RANKED_NONCHIP") - 1 ))
(( N_RANKED > 0 )) || { echo "ERROR: no non-chip candidates after exact-overlap removal" >&2; exit 1; }
echo "===== outer_${OUTER}/inner_${INNER}/${TRAIT} | chip=$BEST_SET ($N_CHIP) | non-chip candidates=$N_RANKED ====="

LAST_POOL=0
LAST_STAGE=""
N_RETAIN=0
for REQUESTED in "${POOL_CHOICES[@]}"; do
    POOL=$REQUESTED
    (( POOL > N_RANKED )) && POOL=$N_RANKED
    (( POOL > LAST_POOL )) || continue
    LAST_POOL=$POOL
    STAGE="$OUT/stages/top${POOL}"
    LAST_STAGE="$STAGE"
    mkdir -p "$STAGE"
    echo "[POOL] top $POOL of $N_RANKED non-chip clump index SNPs"

    "$PYTHON_BIN" "$FILTER" prepare --ranked "$RANKED_NONCHIP" --chip "$CHIP" \
        --train-pheno "$PHENO" --pool "$POOL" --stage "$STAGE"

    LD_BFILE="$STAGE/train_chip_plus_candidates"
    if [[ ! -s "$STAGE/genotype_extract.complete" ]]; then
        "$PLINK_BIN" --bfile "$MASTER" --keep "$STAGE/train.keep" \
            --extract "$STAGE/union.snps" --keep-allele-order --make-bed \
            --threads "$CPUS" --memory 36000 --out "$LD_BFILE"
        for ext in bed bim fam; do [[ -s "${LD_BFILE}.${ext}" ]] || { echo "ERROR: LD genotype output missing ${LD_BFILE}.${ext}" >&2; exit 1; }; done
        printf 'DONE inner-training genotype extraction\n' > "$STAGE/genotype_extract.complete"
    fi

    LD_ARGS=()
    for CHR in $(seq 1 12); do
        PAD=$(printf '%02d' "$CHR")
        CHR_DIR="$STAGE/chr_${PAD}"
        mkdir -p "$CHR_DIR"
        CHR_CHIP="$CHR_DIR/selected_chip.snps"
        awk -v chr="$CHR" 'NR==FNR {selected[$1]=1; next} $1==chr && ($2 in selected) {print $2}' \
            "$STAGE/selected_chip.snps" "${LD_BFILE}.bim" > "$CHR_CHIP"
        N_THIS_CHR=$(wc -l < "$CHR_CHIP")
        if (( N_THIS_CHR == 0 )); then continue; fi
        PREFIX="$CHR_DIR/chip_qtn_sparse_LD"
        LD_GZ="${PREFIX}.ld.gz"
        if [[ -s "${PREFIX}.complete" && -s "$LD_GZ" ]]; then
            gzip -t "$LD_GZ"
        else
            echo "[LD] chromosome $CHR: $N_THIS_CHR chip SNPs, r2 >= $R2"
            "$PLINK_BIN" --bfile "$LD_BFILE" --chr "$CHR" \
                --ld-snp-list "$CHR_CHIP" --r2 gz \
                --ld-window 999999 --ld-window-kb 10000000 \
                --ld-window-r2 "$R2" --threads "$CPUS" \
                --memory 36000 --out "$PREFIX"
            [[ -s "$LD_GZ" ]] || { echo "ERROR: sparse LD file missing: $LD_GZ" >&2; exit 1; }
            gzip -t "$LD_GZ"
            printf 'DONE sparse LD for chromosome %s\n' "$CHR" > "${PREFIX}.complete"
        fi
        LD_ARGS+=(--ld "$LD_GZ")
    done
    (( ${#LD_ARGS[@]} > 0 )) || { echo "ERROR: no selected chip SNP found on chromosomes 1-12" >&2; exit 1; }

    "$PYTHON_BIN" "$FILTER" report --ranked "$RANKED_NONCHIP" --stage "$STAGE" \
        --bim "${LD_BFILE}.bim" --fam "${LD_BFILE}.fam" --train-pheno "$PHENO" \
        --pool "$POOL" --threshold "$R2" "${LD_ARGS[@]}"
    N_RETAIN=$(awk -F '\t' '$1=="Retained" {print $2}' "$STAGE/filter_summary.tsv")
    [[ "$N_RETAIN" =~ ^[0-9]+$ ]] || { echo "ERROR: cannot parse retained count" >&2; exit 1; }
    echo "[RESULT] pool=$POOL retained_after_crosschip_LD=$N_RETAIN"
    if (( N_RETAIN >= 700 )); then
        break
    fi
    (( POOL < N_RANKED )) || break
done

[[ -n "$LAST_STAGE" ]] || { echo "ERROR: no candidate pool was processed" >&2; exit 1; }
"$PYTHON_BIN" "$FILTER" publish --stage "$LAST_STAGE" --out "$OUT" \
    --outer "$OUTER" --inner "$INNER" --trait "$TRAIT" \
    --chip-set "$BEST_SET" --threshold "$R2"

if (( N_RETAIN < 700 )); then
    echo "WARNING: only $N_RETAIN QTN candidates remain after cross-chip LD filtering; some preset counts will be unavailable." >&2
else
    echo "PASS: at least 700 QTN candidates available after cross-chip LD filtering."
fi
