#!/bin/bash
set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 09d verification at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
ROOT="$WORK/Cross_family/results/09_outer_imputed_gwas"
INPUT_ROOT="$WORK/Cross_family/results/09_outer_imputed_gwas_inputs"

PASS=0
FAIL=0

for O in 01 02 03 04 05; do
    OUTER_OK=1
    EXPECTED_SUM=0
    TESTED_SUM=0
    VALID_SUM=0
    INVALID_SUM=0
    echo "===== outer_${O} ====="

    INPUT="$INPUT_ROOT/outer_${O}"
    INPUT_FLAG="$INPUT/complete.flag"
    BFILE="$INPUT/imputed_train_qc"

    if [[ ! -s "$INPUT_FLAG" || ! -s "${BFILE}.bim" ]]; then
        echo "FAIL outer_${O}: missing 09c input flag or imputed BIM"
        OUTER_OK=0
    fi

    for CHR in $(seq 1 12); do
        PAD="$(printf '%02d' "$CHR")"
        D="$ROOT/outer_${O}/chr_${PAD}"
        FLAG="$D/complete.flag"
        QC="$D/gwas_qc.tsv"

        if [[ ! -s "$FLAG" || ! -s "$QC" ]]; then
            echo "FAIL outer_${O}/chr_${PAD}: missing complete.flag or gwas_qc.tsv"
            OUTER_OK=0
            continue
        fi

        for ENTRY in \
            "status=PASS" \
            "version=1" \
            "outer=${O}" \
            "chromosome=${CHR}" \
            "chromosome_padded=${PAD}" \
            "traits=BB,H6" \
            "model=GCTA_MLMA_with_mlma_subtract_grm" \
            "outer_test_phenotype_read=NO"
        do
            if ! grep -qx "$ENTRY" "$FLAG"; then
                echo "FAIL outer_${O}/chr_${PAD}: flag mismatch: $ENTRY"
                OUTER_OK=0
            fi
        done

        if ! awk -F '\t' '
            NR==1 {
                if (NF!=7 || $1!="Outer" || $2!="Trait" || $3!="Chromosome" ||
                    $4!="Expected_SNPs" || $5!="Tested_SNPs" ||
                    $6!="Valid_P" || $7!="Invalid_P") bad=1
                next
            }
            {
                if (NF!=7 || ($2!="BB" && $2!="H6") ||
                    $5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/ ||
                    $7 !~ /^[0-9]+$/ || ($6+$7)!=$5 || $5>$4) bad=1
                if ($2=="BB") bb++
                if ($2=="H6") h6++
                n++
            }
            END {if (bad || n!=2 || bb!=1 || h6!=1) exit 1}
        ' "$QC"; then
            echo "FAIL outer_${O}/chr_${PAD}: malformed QC or not exactly one BB and one H6 row"
            OUTER_OK=0
            continue
        fi

        while IFS=$'\t' read -r ROW_OUTER TRAIT ROW_CHR EXP TESTED VALID INVALID; do
            [[ "$ROW_OUTER" == "Outer" ]] && continue
            if [[ "$ROW_OUTER" != "$O" || "$ROW_CHR" != "$CHR" ]]; then
                echo "FAIL outer_${O}/chr_${PAD}: QC outer/chromosome mismatch"
                OUTER_OK=0
            fi
            RESULT="$D/${TRAIT}_chr${PAD}.mlma"
            if [[ ! -s "$RESULT" ]]; then
                echo "FAIL outer_${O}/chr_${PAD}: missing $RESULT"
                OUTER_OK=0
                continue
            fi
            ROWS="$(awk 'END{print NR-1}' "$RESULT")"
            if [[ "$ROWS" != "$TESTED" ]]; then
                echo "FAIL outer_${O}/chr_${PAD}/${TRAIT}: result rows=$ROWS, QC tested=$TESTED"
                OUTER_OK=0
            fi
            TESTED_SUM=$((TESTED_SUM + TESTED))
            VALID_SUM=$((VALID_SUM + VALID))
            INVALID_SUM=$((INVALID_SUM + INVALID))
            if [[ "$TRAIT" == "BB" ]]; then
                EXPECTED_SUM=$((EXPECTED_SUM + EXP))
            fi
        done < "$QC"
    done

    # Check the chromosome-level expected counts against the autosomal BIM rows.
    # This deliberately checks chromosomes 1-12 only, matching the MLMA jobs.
    if [[ -s "${BFILE}.bim" ]]; then
        BIM_AUTOSOMAL="$(awk '$1 ~ /^[0-9]+$/ && $1>=1 && $1<=12 {n++} END{print n+0}' "${BFILE}.bim")"
        if [[ "$EXPECTED_SUM" -ne "$BIM_AUTOSOMAL" ]]; then
            echo "FAIL outer_${O}: summed expected chr1-12 SNPs=${EXPECTED_SUM}, BIM chr1-12 SNPs=${BIM_AUTOSOMAL}"
            OUTER_OK=0
        fi
    fi

    if (( OUTER_OK )); then
        echo "PASS outer_${O}: 12 chromosomes, 24 trait-GWAS outputs, expected chr1-12 SNPs=${EXPECTED_SUM}, tested=${TESTED_SUM}, valid_P=${VALID_SUM}, invalid_P=${INVALID_SUM}"
        ((PASS+=1))
    else
        echo "FAIL outer_${O}"
        ((FAIL+=1))
    fi
done

echo
echo "OVERALL 09d: PASS=${PASS}/5 FAIL=${FAIL}/5"
(( FAIL == 0 ))
