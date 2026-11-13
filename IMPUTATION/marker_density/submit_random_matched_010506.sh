#!/bin/bash
set -euo pipefail

mkdir -p logs

echo -e "Chr\tRep\tPrepare\tImpute\tEval" > random_matched_chr010506_jobs.tsv

for c in 01 05 06
do
    case "${c}" in
        01) ARRAY_RANGE="0-16" ;;
        05) ARRAY_RANGE="0-14" ;;
        06) ARRAY_RANGE="0-10" ;;
    esac

    for r in 1 2 3
    do
        p=$(sbatch --parsable \
            38_prepare_random_matched_chr.sh \
            "${c}" "${r}")

        i=$(sbatch --parsable \
            --dependency=afterok:${p} \
            --array=${ARRAY_RANGE} \
            39_impute_random_matched_chr.sh \
            "${c}" "${r}")

        e=$(sbatch --parsable \
            --dependency=afterok:${i} \
            40_eval_random_matched_chr.sh \
            "${c}" "${r}")

        echo -e "${c}\t${r}\t${p}\t${i}\t${e}" \
            | tee -a random_matched_chr010506_jobs.tsv
    done
done
