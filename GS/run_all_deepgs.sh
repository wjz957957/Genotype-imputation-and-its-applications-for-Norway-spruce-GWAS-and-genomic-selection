#!/bin/bash

cd ~/GS-904-final

mkdir -p logs

echo "=================================================="
echo "        GS-904 DeepGS"
echo "        CHIP vs CHIP+QTN"
echo "=================================================="
echo "Start: $(date)"
echo

N_RUN=0
N_SKIP=0

for TRAIT in BB H6; do
    for GENO in chip combined; do
        for O in {01..10}; do

            COMPLETE="deepgs_results/${TRAIT}/${GENO}/outer_${O}.complete"

            if [ -f "$COMPLETE" ]; then
                echo "[SKIP] $TRAIT / $GENO / outer_${O}"
                N_SKIP=$((N_SKIP + 1))
                continue
            fi

            echo
            echo "=================================================="
            echo "[RUN] $TRAIT / $GENO / outer_${O}"
            echo "Time: $(date)"
            echo "=================================================="

            N_RUN=$((N_RUN + 1))

            docker run --rm \
                --user "$(id -u):$(id -g)" \
                -w /home/data \
                -v "$PWD":/home/data \
                malab/deepgs_cpu \
                Rscript /home/data/run_deepgs_one.R \
                "$TRAIT" "$GENO" "$O"

            STATUS=$?

            if [ $STATUS -ne 0 ]; then
                echo
                echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                echo "[ERROR] $TRAIT / $GENO / outer_${O}"
                echo "[ERROR] Exit status = $STATUS"
                echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                exit $STATUS
            fi

            if [ ! -f "$COMPLETE" ]; then
                echo
                echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                echo "[ERROR] .complete file missing"
                echo "$TRAIT / $GENO / outer_${O}"
                echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                exit 1
            fi

            echo "[DONE] $TRAIT / $GENO / outer_${O}"
            echo "Time: $(date)"

        done
    done
done

echo
echo "=================================================="
echo "                 ALL DONE"
echo "=================================================="
echo "New jobs run : $N_RUN"
echo "Jobs skipped  : $N_SKIP"
echo "End time     : $(date)"
echo "=================================================="
