#!/bin/bash

set -euo pipefail

OUT=marker_density_accuracy_summary.tsv

echo -e \
"Chr\tDensity_SNP_Mb\tRep\tInput_SNPs\tPooled_R2\tMean_variant_dosage_R2\tGP_concordance\tHard_concordance\tNRC" \
> "${OUT}"


for CHR in 01 05 06 11
do

    for D in 1 5 10 15 20 25 30 40 50 75 100
    do

        for REP in 1 2 3
        do

            DIR=Chr${CHR}/D${D}/rep${REP}

            FILE=${DIR}/Chr${CHR}.D${D}.rep${REP}.imputed_only.accuracy.tsv

            SITES=${DIR}/Chr${CHR}.D${D}.rep${REP}.sites.tsv


            if [[ ! -f "${FILE}" ]]; then
                continue
            fi


            N=$(wc -l < "${SITES}")


            get_value () {
                awk -v key="$1" \
                    '$1==key {print $2}' \
                    "${FILE}"
            }


            POOLED=$(get_value "Pooled_empirical_R2")

            DOSAGE=$(get_value "Mean_per_variant_dosage_R2")

            GPCONC=$(get_value "Mean_GP_weighted_concordance")

            HARDCONC=$(get_value "Hard_genotype_concordance")

            NRC=$(get_value "Non_reference_concordance")


            echo -e \
"${CHR}\t${D}\t${REP}\t${N}\t${POOLED}\t${DOSAGE}\t${GPCONC}\t${HARDCONC}\t${NRC}" \
            >> "${OUT}"

        done

    done

done


echo "Results written to ${OUT}"

echo -n "Conditions found: "
tail -n +2 "${OUT}" | wc -l
