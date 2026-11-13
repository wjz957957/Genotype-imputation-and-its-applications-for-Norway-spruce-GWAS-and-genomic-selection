#!/usr/bin/env bash
set -u

cd "$HOME/GS-904-final" || exit 1
mkdir -p logs

echo "=================================================="
echo " RR-BLUP: original CHIP baseline vs expanded CHIP+QTN"
echo " Start: $(date)"
echo "=================================================="

# Do not start unless the baseline and input files are present.
for O in {01..10}; do
    for TRAIT in BB H6; do
        [[ -s "results/rrBLUP/${TRAIT}/chip/outer_${O}_metrics.tsv" ]] || {
            echo "[ERROR] Missing old CHIP metrics for ${TRAIT} outer_${O}"; exit 1; }
        [[ -s "results/rrBLUP/${TRAIT}/chip/outer_${O}_predictions.tsv" ]] || {
            echo "[ERROR] Missing old CHIP predictions for ${TRAIT} outer_${O}"; exit 1; }
        [[ -f "07_final_genotypes_crosschip_r2_0p8/outer_${O}/${TRAIT}/combined_train.raw" ]] || {
            echo "[ERROR] Missing combined_train.raw for ${TRAIT} outer_${O}. Run make_combined_crosschip_r08.py first."; exit 1; }
        [[ -f "07_final_genotypes_crosschip_r2_0p8/outer_${O}/${TRAIT}/combined_test.raw" ]] || {
            echo "[ERROR] Missing combined_test.raw for ${TRAIT} outer_${O}. Run make_combined_crosschip_r08.py first."; exit 1; }
    done
done

N_RUN=0
N_SKIP=0

for TRAIT in BB H6; do
    for O in {01..10}; do
        D="results_crosschip_r2_0p8/rrBLUP/${TRAIT}/combined"
        COMPLETE="${D}/outer_${O}.complete"
        METRICS="${D}/outer_${O}_metrics.tsv"
        PREDS="${D}/outer_${O}_predictions.tsv"

        if [[ -s "$COMPLETE" && -s "$METRICS" && -s "$PREDS" ]]; then
            echo "[SKIP] ${TRAIT} combined outer_${O}"
            N_SKIP=$((N_SKIP + 1))
            continue
        fi

        echo
        echo "[RUN] RR-BLUP ${TRAIT} combined outer_${O} at $(date)"
        N_RUN=$((N_RUN + 1))

        docker run --rm \
            --user "$(id -u):$(id -g)" \
            -w /home/data \
            -v "$PWD":/home/data \
            malab/deepgs_cpu \
            Rscript /home/data/run_rrblup_crosschip_r08.R \
            "$TRAIT" combined "$O"

        STATUS=$?
        if [[ "$STATUS" -ne 0 ]]; then
            echo "[ERROR] RR-BLUP failed: ${TRAIT} combined outer_${O}, status=${STATUS}"
            exit "$STATUS"
        fi

        if [[ ! -s "$COMPLETE" || ! -s "$METRICS" || ! -s "$PREDS" ]]; then
            echo "[ERROR] Output missing after run: ${TRAIT} combined outer_${O}"
            exit 1
        fi
    done
done

echo
echo "===== RR-BLUP CROSS-CHIP LD r2=0.8 RUNS FINISHED ====="
echo "New jobs run : $N_RUN"
echo "Jobs skipped : $N_SKIP"
echo "End time     : $(date)"
