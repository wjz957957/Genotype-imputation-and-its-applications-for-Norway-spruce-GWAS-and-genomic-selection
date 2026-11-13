#!/bin/bash
set -Eeuo pipefail
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
VERIFY="$WORK/Cross_family/scripts/09b2_verify_outer_chip_selection.py"
PASS=0
FAIL=0
for O in 01 02 03 04 05; do
    if python3 "$VERIFY" --work "$WORK" --outer "$O" --phase post; then
        ((PASS+=1))
    else
        echo "FAIL 09b2 outer_${O}" >&2
        ((FAIL+=1))
    fi
done
echo "OVERALL 09b2 PASS=${PASS}/5 FAIL=${FAIL}/5"
((FAIL==0 && PASS==5))
