#!/usr/bin/env bash
set -u

cd "$HOME/GS-904-final" || exit 1

OLD="final_genotypes"
NEW="07_final_genotypes_crosschip_r2_0p8"
FAIL=0

echo "===== CROSS-CHIP LD r2=0.8 PREFLIGHT ====="
echo "Old CHIP root: $OLD"
echo "New data root: $NEW"
echo

if [[ ! -d "$OLD" || ! -d "$NEW" ]]; then
    echo "[FAIL] Required genotype directory is missing."
    exit 1
fi

for O in {01..10}; do
    for TRAIT in BB H6; do
        OLD_D="$OLD/outer_${O}/${TRAIT}"
        NEW_D="$NEW/outer_${O}/${TRAIT}"

        for F in chip_train.raw chip_test.raw; do
            if [[ ! -f "$OLD_D/$F" || ! -f "$NEW_D/$F" ]]; then
                echo "[FAIL] Missing CHIP file: $OLD_D/$F or $NEW_D/$F"
                FAIL=1
            elif cmp -s "$OLD_D/$F" "$NEW_D/$F"; then
                echo "[PASS] CHIP identical: ${TRAIT} outer_${O} $F"
            else
                echo "[FAIL] CHIP differs: ${TRAIT} outer_${O} $F"
                FAIL=1
            fi
        done

        for FF in chip_train.raw chip_test.raw extra_train.raw extra_test.raw final_genotype_manifest.tsv; do
            if [[ ! -s "$NEW_D/$FF" ]]; then
                echo "[FAIL] Missing/empty new input: $NEW_D/$FF"
                FAIL=1
            fi
        done

        ET="$NEW_D/extra_train.raw"
        ES="$NEW_D/extra_test.raw"
        if [[ -f "$ET" && ! -f "$ES" ]] || [[ ! -f "$ET" && -f "$ES" ]]; then
            echo "[FAIL] Only one EXTRA train/test file exists: $NEW_D"
            FAIL=1
        fi

        if [[ ! -f "$NEW_D/complete.flag" ]]; then
            echo "[FAIL] Missing 10d completion flag: $NEW_D/complete.flag"
            FAIL=1
        fi

        if [[ ! -f "CV/outer_${O}/outer_train.samples" ||
              ! -f "CV/outer_${O}/outer_test.samples" ]]; then
            echo "[FAIL] Missing CV sample lists for outer_${O}"
            FAIL=1
        fi
    done
done

echo
echo "===== EXISTING CHIP RESULTS ====="
for O in {01..10}; do
    for TRAIT in BB H6; do
        R1="results/rrBLUP/${TRAIT}/chip/outer_${O}_metrics.tsv"
        R2="results/rrBLUP/${TRAIT}/chip/outer_${O}_predictions.tsv"
        D1="deepgs_results/${TRAIT}/chip/outer_${O}_metrics.tsv"
        D2="deepgs_results/${TRAIT}/chip/outer_${O}_predictions.tsv"
        D3="deepgs_results/${TRAIT}/chip/outer_${O}.complete"

        for F in "$R1" "$R2" "$D1" "$D2" "$D3"; do
            if [[ ! -s "$F" ]]; then
                echo "[FAIL] Missing/empty baseline result: $F"
                FAIL=1
            fi
        done
    done
done

echo
if [[ "$FAIL" -ne 0 ]]; then
    echo "[STOP] Preflight failed."
    echo "Do not reuse the old CHIP baseline or start the combined-only workflow."
    echo "Resolve the mismatches/missing files first."
    exit 1
fi

echo "[PASS] All new CHIP train/test files are byte-identical to the original CHIP files."
echo "[PASS] Old CHIP results and required CV files are present."
echo "Safe to reuse the old CHIP baseline and run only the new CHIP+QTN models."
