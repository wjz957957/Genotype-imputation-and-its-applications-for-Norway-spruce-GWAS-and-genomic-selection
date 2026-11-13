#!/usr/bin/env bash
# Publish the previously completed H6 outer_01 inner_01 LD filtering stage.
# Does not rerun genotype extraction or PLINK LD.
set -Eeuo pipefail

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
PY="$WORK/scripts/14b_filter_H6_inner_qtn_LD_top700.py"
FINAL="$WORK/08_imputed_GWAS/08_crosschip_LD_pilot/outer_01/inner_01/H6/r2_0p8"
STAGE="$FINAL/stages/top9182"

for F in "$PY" "$STAGE/stage.complete" "$STAGE/filtered_ranked_extra_snps.tsv" "$STAGE/filter_summary.tsv"; do
  [[ -s "$F" ]] || { echo "ERROR: prerequisite missing or empty: $F" >&2; exit 1; }
done

if [[ -s "$FINAL/complete.flag" ]]; then
  echo "ERROR: already published: $FINAL; refusing to replace it" >&2
  exit 1
fi

N_RETAIN=$(awk -F '\t' '$1 == "Retained" {print $2}' "$STAGE/filter_summary.tsv")
[[ "$N_RETAIN" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid retained count: $N_RETAIN" >&2; exit 1; }
(( N_RETAIN >= 700 )) || { echo "ERROR: fewer than 700 retained SNPs ($N_RETAIN)" >&2; exit 1; }

echo "[REUSE] H6 outer_01 inner_01 stage pool=9182 retained=$N_RETAIN"
python3 "$PY" publish --stage "$STAGE" --out "$FINAL"

[[ $(wc -l < "$FINAL/top700_extra.snps") -eq 700 ]] || { echo "ERROR: top700 has incorrect count" >&2; exit 1; }
for N in 1 2 3 4 5 10 20 30 40 50 75 100 150 200 300 400 500 600 700; do
  [[ $(wc -l < "$FINAL/marker_sets/extra_${N}.snps") -eq $N ]] || { echo "ERROR: extra_${N}.snps has incorrect count" >&2; exit 1; }
done
[[ -s "$FINAL/complete.flag" ]] || { echo "ERROR: missing completion marker" >&2; exit 1; }

echo "DONE: inner_01 published with new 19-count 1..700 grid (plus 0 baseline)"
