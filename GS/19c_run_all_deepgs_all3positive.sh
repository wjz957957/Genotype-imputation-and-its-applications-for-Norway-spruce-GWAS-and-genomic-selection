#!/usr/bin/env bash
set -Eeuo pipefail
cd "$HOME/GS-904-final"
mkdir -p logs

DATA=08_final_genotypes_crosschip_r2_0p8_all3positive
RES=deepgs_results_crosschip_r2_0p8_all3positive
for O in {01..10}; do
  [[ -s "deepgs_results/validation/outer_${O}_valid.samples" ]] || { echo "[ERROR] Missing existing validation IDs outer_$O" >&2; exit 1; }
  for TRAIT in BB H6; do
    D="$DATA/outer_${O}/${TRAIT}"
    for f in chip_train.raw chip_test.raw combined_train.raw combined_test.raw final_genotype_manifest.tsv; do
      [[ -s "$D/$f" ]] || { echo "[ERROR] Missing input $D/$f" >&2; exit 1; }
    done
    for f in metrics.tsv predictions.tsv; do
      [[ -s "deepgs_results/${TRAIT}/chip/outer_${O}_${f}" ]] || { echo "[ERROR] Missing CHIP baseline ${TRAIT} outer_${O} (${f})" >&2; exit 1; }
    done
    [[ -s "deepgs_results/${TRAIT}/chip/outer_${O}.complete" ]] || { echo "[ERROR] Missing CHIP baseline complete ${TRAIT} outer_${O}" >&2; exit 1; }
  done
done

N_RUN=0; N_SKIP=0
for TRAIT in BB H6; do
  for O in {01..10}; do
    D="${RES}/${TRAIT}/combined"
    if [[ -s "$D/outer_${O}.complete" && -s "$D/outer_${O}_metrics.tsv" && -s "$D/outer_${O}_predictions.tsv" ]]; then
      echo "[SKIP] $TRAIT outer_$O already completed"
      N_SKIP=$((N_SKIP+1))
      continue
    fi
    echo "[RUN] $TRAIT outer_$O at $(date)"
    docker run --rm \
        --user "$(id -u):$(id -g)" \
        -w /home/data \
        -v "$PWD":/home/data \
        malab/deepgs_cpu \
        Rscript /home/data/19c_run_deepgs_all3positive.R "$TRAIT" combined "$O"
    for suffix in _metrics.tsv _predictions.tsv .complete; do
      [[ -s "$D/outer_${O}${suffix}" ]] || { echo "[ERROR] Missing output $D/outer_${O}${suffix}" >&2; exit 1; }
    done
    N_RUN=$((N_RUN+1))
  done
done
echo "COMPLETE: ran=$N_RUN skipped=$N_SKIP"
