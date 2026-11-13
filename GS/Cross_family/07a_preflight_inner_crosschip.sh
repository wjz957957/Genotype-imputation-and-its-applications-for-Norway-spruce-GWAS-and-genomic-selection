#!/bin/bash
set -Eeuo pipefail
export LC_ALL=C

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
CONFIG="$WORK/Cross_family/config/best_chip_by_outer_trait.tsv"
MASTER="$WORK/08_imputed_GWAS/00_reference/all_904_genotypes"
PHENO_ROOT="$WORK/Cross_family/results/02_inner_phenotypes"
CHIP_ROOT="$WORK/Cross_family/results/03_inner_chip_gwas_clumping"
MERGED_ROOT="$WORK/Cross_family/results/06_inner_imputed_clumping_merged"
MERGE_SUMMARY="$MERGED_ROOT/merge_summary.tsv"

for f in "$CONFIG" "$MERGE_SUMMARY" "${MASTER}.bed" "${MASTER}.bim" "${MASTER}.fam"; do
  [[ -s "$f" ]] || { echo "FAIL: missing/empty $f" >&2; exit 1; }
done

HEADER=$(head -n 1 "$CONFIG")
[[ "$HEADER" == $'Outer\tTrait\tBest_Chip\tMean_Fold_Pearson_PA' ]] || {
  echo "FAIL: unexpected config header: $HEADER" >&2; exit 1;
}

N_PASS=0
for O in 01 02 03 04 05; do
  for I in 01 02 03; do
    FOLD="outer_${O}/inner_${I}"
    for T in BB H6; do
      BEST=$(awk -F '\t' -v o="$O" -v t="$T" 'NR>1 && $1==o && $2==t {print $3}' "$CONFIG")
      N_CONFIG=$(awk -F '\t' -v o="$O" -v t="$T" 'NR>1 && $1==o && $2==t {n++} END{print n+0}' "$CONFIG")
      [[ "$N_CONFIG" -eq 1 ]] || { echo "FAIL $FOLD/$T: expected one config row, got $N_CONFIG" >&2; exit 1; }
      [[ "$BEST" =~ ^chip_(50|100|200|300|400|500|1000|2000|3000|4000|5000|8000|10000|20000|30000|40000)$ ]] || {
        echo "FAIL $FOLD/$T: unsupported selected chip strategy '$BEST'" >&2; exit 1;
      }
      CHIP="$CHIP_ROOT/$FOLD/marker_sets/${T}_${BEST}.snps"
      PHENO="$PHENO_ROOT/$FOLD/${T}.train.pheno"
      RANKED="$MERGED_ROOT/$FOLD/${T}_index_ranked.tsv"
      for f in "$CHIP" "$PHENO" "$RANKED"; do
        [[ -s "$f" ]] || { echo "FAIL $FOLD/$T: missing/empty $f" >&2; exit 1; }
      done
      N_SUMMARY=$(awk -F '\t' -v o="$O" -v i="$I" -v t="$T" '
        NR==1 {if (tolower($1)!="outer" || tolower($2)!="inner" || tolower($3)!="trait" || tolower($4)!="n_candidate_snps") exit 2; next}
        $1==o && $2==i && $3==t {n++; count=$4}
        END {if (n!=1 || count !~ /^[0-9]+$/) exit 3; print count}
      ' "$MERGE_SUMMARY") || { echo "FAIL $FOLD/$T: no unique row in merge_summary.tsv" >&2; exit 1; }
      EXPECTED=${BEST#chip_}
      N_CHIP=$(wc -l < "$CHIP")
      [[ "$N_CHIP" -eq "$EXPECTED" ]] || {
        echo "FAIL $FOLD/$T: $BEST has $N_CHIP SNPs, expected $EXPECTED ($CHIP)" >&2; exit 1;
      }
      RANK_HEADER=$(head -n 1 "$RANKED")
      [[ "$RANK_HEADER" == $'SNP\tP' ]] || { echo "FAIL $FOLD/$T: unexpected ranked header '$RANK_HEADER'" >&2; exit 1; }
      N_RANKED=$(( $(wc -l < "$RANKED") - 1 ))
      [[ "$N_RANKED" -eq "$N_SUMMARY" ]] || { echo "FAIL $FOLD/$T: ranked rows $N_RANKED != merge_summary count $N_SUMMARY" >&2; exit 1; }
      (( N_RANKED > 0 )) || { echo "FAIL $FOLD/$T: no clumped candidate SNPs" >&2; exit 1; }
      awk 'NF!=3 || $1!=$2 {exit 1}' "$PHENO" || {
        echo "FAIL $FOLD/$T: phenotype file must be headerless FID IID Y with FID=IID" >&2; exit 1;
      }
      N_PHENO=$(wc -l < "$PHENO")
      (( N_PHENO >= 100 )) || { echo "FAIL $FOLD/$T: unexpected training phenotype count $N_PHENO" >&2; exit 1; }
      printf 'PASS %s/%s: chip=%s (%s SNPs), training=%s, ranked_candidates=%s\n' \
        "$FOLD" "$T" "$BEST" "$N_CHIP" "$N_PHENO" "$N_RANKED"
      ((N_PASS+=1))
    done
  done
done

[[ "$N_PASS" -eq 30 ]] || { echo "FAIL: validated $N_PASS/30 fold-trait tasks" >&2; exit 1; }
echo "OVERALL PASS: all 30 inner-fold/trait inputs are ready for cross-chip LD filtering."
