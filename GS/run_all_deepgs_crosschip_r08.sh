#!/usr/bin/env bash
set -u

cd "$HOME/GS-904-final" || exit 1
mkdir -p logs

echo "=================================================="
echo " DeepGS: original CHIP baseline vs expanded CHIP+QTN"
echo " Start: $(date)"
echo "=================================================="

# The new combined-only workflow relies on existing CHIP results and validation splits.
for O in {01..10}; do
    for TRAIT in BB H6; do
        [[ -s "deepgs_results/${TRAIT}/chip/outer_${O}_metrics.tsv" ]] || {
            echo "[ERROR] Missing old CHIP metrics for ${TRAIT} outer_${O}"; exit 1; }
        [[ -s "deepgs_results/${TRAIT}/chip/outer_${O}_predictions.tsv" ]] || {
            echo "[ERROR] Missing old CHIP predictions for ${TRAIT} outer_${O}"; exit 1; }
        [[ -s "deepgs_results/${TRAIT}/chip/outer_${O}.complete" ]] || {
            echo "[ERROR] Missing old CHIP completion flag for ${TRAIT} outer_${O}"; exit 1; }
        [[ -s "deepgs_results/validation/outer_${O}_valid.samples" ]] || {
            echo "[ERROR] Missing original validation IDs for outer_${O}"; exit 1; }
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
        D="deepgs_results_crosschip_r2_0p8/${TRAIT}/combined"
        COMPLETE="${D}/outer_${O}.complete"
        METRICS="${D}/outer_${O}_metrics.tsv"
        PREDS="${D}/outer_${O}_predictions.tsv"

        if [[ -s "$COMPLETE" && -s "$METRICS" && -s "$PREDS" ]]; then
            echo "[SKIP] ${TRAIT} combined outer_${O}"
            N_SKIP=$((N_SKIP + 1))
            continue
        fi

        echo
        echo "[RUN] DeepGS ${TRAIT} combined outer_${O} at $(date)"
        N_RUN=$((N_RUN + 1))

        docker run --rm \
            --user "$(id -u):$(id -g)" \
            -w /home/data \
            -v "$PWD":/home/data \
            malab/deepgs_cpu \
            Rscript /home/data/run_deepgs_crosschip_r08.R \
            "$TRAIT" combined "$O"

        STATUS=$?
        if [[ "$STATUS" -ne 0 ]]; then
            echo "[ERROR] DeepGS failed: ${TRAIT} combined outer_${O}, status=${STATUS}"
            exit "$STATUS"
        fi

        if [[ ! -s "$COMPLETE" || ! -s "$METRICS" || ! -s "$PREDS" ]]; then
            echo "[ERROR] Output missing after run: ${TRAIT} combined outer_${O}"
            exit 1
        fi
    done
done

echo
echo "===== DEEPGS CROSS-CHIP LD r2=0.8 RUNS FINISHED ====="
echo "New jobs run : $N_RUN"
echo "Jobs skipped : $N_SKIP"
echo "End time     : $(date)"
