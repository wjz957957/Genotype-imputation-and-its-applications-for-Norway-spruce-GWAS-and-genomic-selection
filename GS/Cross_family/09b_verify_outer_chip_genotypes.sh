#!/bin/bash
set -Eeuo pipefail
WORK="${1:-/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer}"
PY=$(command -v python3)
CHECK="$WORK/Cross_family/scripts/09b_validate_outer_chip_inputs.py"
pass=0
fail=0
for outer in 01 02 03 04 05; do
  d="$WORK/Cross_family/results/09_outer_chip_genotypes/outer_${outer}"
  pref="$d/chip_ALL_outer_train"
  if [[ -s "$d/complete.flag" ]] && grep -qx 'status=PASS' "$d/complete.flag" && \
     "$PY" "$CHECK" --work "$WORK" --outer "$outer" --out "$d" --prefix "$pref" --phase post; then
    pass=$((pass+1))
  else
    echo "FAIL outer_${outer}"
    fail=$((fail+1))
  fi
done
echo "OVERALL PASS=$pass/5 FAIL=$fail/5"
((pass==5 && fail==0))
