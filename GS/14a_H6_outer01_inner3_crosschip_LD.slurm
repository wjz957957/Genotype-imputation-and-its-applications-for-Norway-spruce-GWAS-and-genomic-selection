#!/bin/bash
#SBATCH -J 14a_H601_innerLD
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-2%1
#SBATCH -o logs/14a_H601_innerLD_%A_%a.out
#SBATCH -e logs/14a_H601_innerLD_%A_%a.err

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR line ${LINENO}, status=${rc}: ${BASH_COMMAND}" >&2' ERR
export LC_ALL=C

module load bioinfo-tools
ml plink/1.90b4.9

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
PY="${WORK}/scripts/14a_filter_H6_inner_qtn_ld.py"
MASTER="${WORK}/08_imputed_GWAS/00_reference/all_904_genotypes"
OUTER=01
INNER=$(printf '%02d' "$(( SLURM_ARRAY_TASK_ID + 1 ))")
TRAIT=H6
R2=${R2:-0.8}

# One threshold gets its own directory; does not affect any previous analysis.
case "$R2" in
  0.8|0.5|0.2) ;;
  *) echo "ERROR: R2 must equal 0.8, 0.5, or 0.2" >&2; exit 2;;
esac
R2_TAG=${R2//./p}
FINAL="${WORK}/08_imputed_GWAS/08_crosschip_LD_pilot/outer_${OUTER}/inner_${INNER}/${TRAIT}/r2_${R2_TAG}"
RANKED="${WORK}/08_imputed_GWAS/07_imputed_marker_sets/outer_${OUTER}/inner_${INNER}/${TRAIT}/ranked_imputed_snps.tsv"
BEST_FILE="${WORK}/07_chip_RRBLUP/optimization/outer_${OUTER}/${TRAIT}_best_marker.txt"
CV="${WORK}/02_CV/outer_${OUTER}/inner_${INNER}"
TRAIN="${CV}/inner_train.samples"
VALID="${CV}/inner_valid.samples"
OUTER_TEST="${WORK}/02_CV/outer_${OUTER}/outer_test.samples"

for F in "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam" "$PY" \
  "$RANKED" "$BEST_FILE" "$TRAIN" "$VALID" "$OUTER_TEST"; do
    [[ -s "$F" ]] || { echo "ERROR: missing input $F" >&2; exit 1; }
done

BEST=$(tr -d '\r\n[:space:]' < "$BEST_FILE")
[[ "$BEST" =~ ^chip_(50|100|200|300|400|500|1000|2000|3000|4000|5000|8000|10000|20000|30000|40000|LD_ALL|ALL)$ ]] || {
  echo "ERROR: unexpected H6 best marker setting: '$BEST'" >&2; exit 1;
}
CHIP="${WORK}/06_chip_clumping/outer_${OUTER}/inner_${INNER}/${TRAIT}/marker_sets/${BEST}.snps"
[[ -s "$CHIP" ]] || { echo "ERROR: missing inner_train selected CHIP list $CHIP" >&2; exit 1; }

N_CHIP=$(wc -l < "$CHIP")
(( N_CHIP >= 30000 && N_CHIP <= 50000 )) || {
  echo "ERROR: unexpected selected H6 CHIP count $N_CHIP; expected ~40k" >&2; exit 1;
}
N_RANKED=$(( $(wc -l < "$RANKED") - 1 ))
(( N_RANKED >= 1000 )) || { echo "ERROR: only $N_RANKED ranked imputed SNPs" >&2; exit 1; }

mkdir -p "$FINAL"
if [[ -s "$FINAL/complete.flag" && -s "$FINAL/filter_summary.tsv" \
      && -s "$FINAL/top1000_extra.snps" ]]; then
  N_FINAL=$(wc -l < "$FINAL/top1000_extra.snps")
  [[ "$N_FINAL" -eq 1000 ]] || { echo 'ERROR: invalid prior top1000' >&2; exit 1; }
  echo "[SKIP] Completed: $FINAL"; exit 0
fi

printf '\n===== H6 outer_01 inner_%s / cross CHIP-imputed LD >= %s =====\n' "$INNER" "$R2"
echo "Best CHIP setting: $BEST ($N_CHIP SNPs)"
echo "Full GWAS-clumped candidate list: $N_RANKED SNPs"
echo 'Training samples only; no GWAS/clumping rerun; outputs in separate H6 pilot directory.'

# As in the BB pilot, expand only until 1000 QTN remain. There may be fewer
# than 10000 original index SNPs in the complete ranked file.
POOL_CHOICES=(10000 20000 40000 80000 100000)
LAST_POOL=0
for REQUESTED in "${POOL_CHOICES[@]}"; do
  POOL=$REQUESTED
  (( POOL > N_RANKED )) && POOL=$N_RANKED
  (( POOL > LAST_POOL )) || continue
  LAST_POOL=$POOL
  STAGE="${FINAL}/stages/top${POOL}"
  mkdir -p "$STAGE"
  echo "===== Candidate pool $POOL of $N_RANKED ====="

  python3 "$PY" prepare \
    --ranked "$RANKED" --chip "$CHIP" \
    --train "$TRAIN" --valid "$VALID" --outer-test "$OUTER_TEST" \
    --pool "$POOL" --stage "$STAGE"

  LD_BFILE="${STAGE}/train_chip_plus_candidates"
  # IMPORTANT: do not apply --maf/--geno to this extraction. H6 uses nearly
  # all chip SNPs, including some with low MAF in this training fold; silently
  # dropping markers would make the reference list inconsistent.
  if [[ ! -s "${STAGE}/genotype_extract.complete" ]]; then
    plink --bfile "$MASTER" \
      --keep "${STAGE}/train.keep" \
      --extract "${STAGE}/union.snps" \
      --keep-allele-order --make-bed \
      --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 36000 \
      --out "$LD_BFILE"
    for E in bed bim fam; do
      [[ -s "${LD_BFILE}.${E}" ]] || { echo "ERROR: genotype .${E} missing" >&2; exit 1; }
    done
    echo 'DONE genotype extraction on inner_train' > "${STAGE}/genotype_extract.complete"
  else
    echo "[REUSE] $LD_BFILE"
  fi

  # Per-chromosome calculation prevents treating the ~40k H6 CHIP markers
  # as one huge whole-genome LD task. Only sparse pairs with r2>=R2 are stored.
  # The 10,000,000 kb window covers each whole spruce chromosome; no 1Mb cap.
  LD_ARGS=()
  for CHR in $(seq 1 12); do
    CHR_PAD=$(printf '%02d' "$CHR")
    CHR_DIR="${STAGE}/chr_${CHR_PAD}"
    mkdir -p "$CHR_DIR"
    CHR_CHIP="${CHR_DIR}/selected_chip.snps"
    awk -v chr="$CHR" '
      NR==FNR { chip[$1]=1; next }
      $1==chr && ($2 in chip) { print $2 }
    ' "${STAGE}/selected_chip.snps" "${LD_BFILE}.bim" > "$CHR_CHIP"

    N_CHR_CHIP=$(wc -l < "$CHR_CHIP")
    # Without a CHIP on the candidate chromosome, there is nothing to test.
    if (( N_CHR_CHIP == 0 )); then
      echo "[CHR $CHR_PAD] No selected CHIP; skip LD"
      continue
    fi

    LD_PREFIX="${CHR_DIR}/chip_qtn_sparse_LD"
    LD_GZ="${LD_PREFIX}.ld.gz"
    LD_DONE="${LD_PREFIX}.complete"
    if [[ -s "$LD_DONE" && -s "$LD_GZ" ]]; then
      gzip -t "$LD_GZ"
      echo "[REUSE CHR $CHR_PAD] Sparse LD already completed"
    else
      echo "[LD CHR $CHR_PAD] selected CHIP=$N_CHR_CHIP started $(date)"
      plink --bfile "$LD_BFILE" --chr "$CHR" \
        --ld-snp-list "$CHR_CHIP" \
        --r2 gz --ld-window 999999 --ld-window-kb 10000000 \
        --ld-window-r2 "$R2" \
        --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 36000 \
        --out "$LD_PREFIX"
      [[ -s "$LD_GZ" ]] || { echo "ERROR: LD output missing for chr $CHR" >&2; exit 1; }
      gzip -t "$LD_GZ"
      echo "DONE chr $CHR at $(date)" > "$LD_DONE"
    fi
    LD_ARGS+=(--ld "$LD_GZ")
  done

  (( ${#LD_ARGS[@]} > 0 )) || {
    echo 'ERROR: H6 selected CHIP list contains no chromosome 1-12 markers' >&2; exit 1;
  }

  python3 "$PY" report \
    --ranked "$RANKED" --stage "$STAGE" \
    --bim "${LD_BFILE}.bim" --fam "${LD_BFILE}.fam" \
    "${LD_ARGS[@]}" \
    --train "$TRAIN" --pool "$POOL" --threshold "$R2"

  N_RETAIN=$(awk -F '\t' '$1=="Retained" {print $2}' "$STAGE/filter_summary.tsv")
  [[ "$N_RETAIN" =~ ^[0-9]+$ ]] || { echo 'ERROR: invalid retained count' >&2; exit 1; }
  echo "[PROGRESS] H6 inner_${INNER}: pool=${POOL}, retained=${N_RETAIN}"
  if (( N_RETAIN >= 1000 )); then
    python3 "$PY" publish --stage "$STAGE" --out "$FINAL"
    echo "[DONE] H6 inner_${INNER} $FINAL"
    exit 0
  fi
  if (( POOL == N_RANKED )); then break; fi
done

printf 'pool_tested\t%s\nretained\t%s\nrequired\t1000\nthreshold\t%s\n' \
    "$LAST_POOL" "$N_RETAIN" "$R2" > "${FINAL}/insufficient_candidates.tsv"
echo "[INSUFFICIENT] Only $N_RETAIN survive LD screening (pool=$LAST_POOL)." >&2
echo 'No relaxed threshold or synthetic QTN were introduced.' >&2
exit 3
