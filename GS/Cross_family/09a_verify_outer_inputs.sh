#!/bin/bash
set -Eeuo pipefail
WORK="${1:-/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer}"
ROOT="$WORK/Cross_family/results/09_outer_inputs"
SELECT="$WORK/Cross_family/results/08_inner_qtn_selection/selected_qtn_counts.tsv"
[[ -s "$SELECT" ]] || { echo "ERROR missing selection: $SELECT" >&2; exit 1; }
pass=0;fail=0
for outer in 01 02 03 04 05; do
    dir="$ROOT/outer_${outer}"
    valid=1
    for f in prepare.complete chip_ALL.snps outer_train.keep outer_test.ids BB.train.pheno H6.train.pheno frozen_inner_selection.tsv; do
        [[ -s "$dir/$f" ]] || { echo "FAIL outer_$outer missing $f";valid=0; }
    done
    if [[ -s "$dir/prepare.complete" ]]; then
        for flag in 'status=PASS' 'outer_test_phenotype_read=NO' 'inner_validation_pheno_merged=PASS' 'full_chip_consistency=PASS' 'sample_partition_master=PASS'; do
            grep -qx "$flag" "$dir/prepare.complete" || { echo "FAIL outer_$outer missing flag $flag";valid=0; }
        done
    fi
    if (( valid )); then
        n_chip=$(wc -l < "$dir/chip_ALL.snps")
        n_train=$(wc -l < "$dir/outer_train.keep")
        n_test=$(wc -l < "$dir/outer_test.ids")
        if (( n_chip != 43267 || n_train + n_test != 904 )); then
            echo "FAIL outer_$outer counts chip=$n_chip train=$n_train test=$n_test"
            valid=0
        fi
    fi
    for trait in BB H6; do
        # Compare frozen selections exactly with 08c's count and 06's chip selection.
        frozen="$dir/frozen_inner_selection.tsv"
        if [[ -s "$frozen" ]]; then
            n=$(awk -F '\t' -v t="$trait" 'NR>1&&$2==t{print $4}' "$frozen")
            exp=$(awk -F '\t' -v o="$outer" -v t="$trait" 'NR>1&&$1==o&&$2==t{print $3}' "$SELECT")
            if [[ -z "$n" || -z "$exp" || "$n" != "$exp" ]]; then
                echo "FAIL outer_$outer/$trait frozen count $n != selected $exp";valid=0
            fi
        fi
    done
    if (( valid )); then
        echo "PASS outer_$outer: chip_ALL=$n_chip outer_train=$n_train outer_test=$n_test"
        pass=$((pass+1))
    else
        fail=$((fail+1))
    fi
done
echo "OVERALL: PASS=$pass/5 FAIL=$fail/5"
(( pass==5 && fail==0 ))
