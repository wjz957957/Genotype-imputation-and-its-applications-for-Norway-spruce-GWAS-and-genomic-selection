#!/bin/bash

set -euo pipefail

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

cd "${WORKDIR}"


DENSITIES=(
1
5
10
15
20
25
30
40
50
75
100
)

REPS=(1 2 3)


# Limit simultaneous chunk tasks within each array
MAX_CHUNK_CONCURRENT=4


# ============================================================
# Chromosomes supplied on command line
#
# Examples:
#
# bash 17_submit_marker_density_batch.sh 01
#
# bash 17_submit_marker_density_batch.sh 01 05 06 11
#
# bash 17_submit_marker_density_batch.sh all
# ============================================================

if [[ "$#" -eq 0 ]]; then
    echo "Usage:"
    echo "bash $0 01"
    echo "bash $0 01 05 06 11"
    echo "bash $0 all"
    exit 1
fi


if [[ "$1" == "all" ]]; then
    CHRS=(01 05 06 11)
else
    CHRS=("$@")
fi


mkdir -p logs

SUMMARY=${WORKDIR}/submitted_jobs_$(date +%Y%m%d_%H%M%S).tsv

echo -e \
"Chr\tDensity\tRep\tPrepare_job\tImpute_job\tFinal_job\tChunks" \
> "${SUMMARY}"


for CHR2 in "${CHRS[@]}"
do

    case "${CHR2}" in
        01|05|06|11)
            ;;
        *)
            echo "ERROR: unsupported chromosome ${CHR2}"
            exit 1
            ;;
    esac


    REGION_FILE=${WORKDIR}/Chr${CHR2}/Chr${CHR2}.impute_regions.tsv
    CANDIDATES=${WORKDIR}/Chr${CHR2}/base/Chr${CHR2}.candidate_positions.tsv


    if [[ ! -s "${REGION_FILE}" ]]; then
        echo "ERROR: missing region file ${REGION_FILE}"
        exit 1
    fi


    if [[ ! -s "${CANDIDATES}" ]]; then
        echo "ERROR: missing candidate file ${CANDIDATES}"
        exit 1
    fi


    NCHUNKS=$(wc -l < "${REGION_FILE}")
    LAST_TASK=$((NCHUNKS - 1))


    echo
    echo "===================================================="
    echo "Chr${CHR2}: ${NCHUNKS} chunks"
    echo "===================================================="


    for DENSITY in "${DENSITIES[@]}"
    do

        for REP in "${REPS[@]}"
        do

            echo
            echo "Submitting Chr${CHR2} D${DENSITY} rep${REP}"


            # ------------------------------------------------
            # Prepare + phase
            # ------------------------------------------------

            PREP_JOB=$(sbatch \
                --parsable \
                --job-name="P${CHR2}_D${DENSITY}_R${REP}" \
                14_prepare_phase_condition.sh \
                "${CHR2}" \
                "${DENSITY}" \
                "${REP}")

            PREP_JOB=${PREP_JOB%%;*}


            # ------------------------------------------------
            # IMPUTE5 chromosome chunks
            # ------------------------------------------------

            IMP_JOB=$(sbatch \
                --parsable \
                --dependency=afterok:${PREP_JOB} \
                --array=0-${LAST_TASK}%${MAX_CHUNK_CONCURRENT} \
                --job-name="I${CHR2}_D${DENSITY}_R${REP}" \
                15_impute_condition_array.sh \
                "${CHR2}" \
                "${DENSITY}" \
                "${REP}")

            IMP_JOB=${IMP_JOB%%;*}


            # ------------------------------------------------
            # CONCAT + accuracy
            # ------------------------------------------------

            FINAL_JOB=$(sbatch \
                --parsable \
                --dependency=afterok:${IMP_JOB} \
                --job-name="E${CHR2}_D${DENSITY}_R${REP}" \
                16_concat_evaluate_condition.sh \
                "${CHR2}" \
                "${DENSITY}" \
                "${REP}")

            FINAL_JOB=${FINAL_JOB%%;*}


            echo -e \
"${CHR2}\t${DENSITY}\t${REP}\t${PREP_JOB}\t${IMP_JOB}\t${FINAL_JOB}\t${NCHUNKS}" \
                >> "${SUMMARY}"


            echo \
"prepare=${PREP_JOB} impute=${IMP_JOB} final=${FINAL_JOB}"

        done

    done

done


echo
echo "===================================================="
echo "Submission completed."
echo "Job table:"
echo "${SUMMARY}"
echo "===================================================="
