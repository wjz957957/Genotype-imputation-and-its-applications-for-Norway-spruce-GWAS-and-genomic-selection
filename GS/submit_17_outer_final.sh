#!/bin/bash
set -Eeuo pipefail
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
mkdir -p logs
J17A=$(sbatch --parsable scripts/17a_filter_outer_crosschip_LD.slurm)
J17A=${J17A%%;*}
J17B=$(sbatch --parsable --dependency="afterok:${J17A}" scripts/17b_prepare_outer_final_genotypes.slurm)
J17B=${J17B%%;*}
printf 'Submitted 17a=%s (outer_train LD) and 17b=%s (final genotype, afterok:17a)\n' "$J17A" "$J17B"
echo 'Wait until all Step17b array tasks succeed; then: python3 scripts/17c_verify_outer_final.py'
