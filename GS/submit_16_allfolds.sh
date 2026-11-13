#!/usr/bin/env bash
# Submit the full all-fold nested-CV re-screening pipeline with afterok dependencies.
set -Eeuo pipefail
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
mkdir -p logs
for F in \
  scripts/16a_allfolds_crosschip_LD.slurm \
  scripts/16b_allfolds_prepare_RAW.slurm \
  scripts/16c_allfolds_crosschip_RRBLUP.slurm \
  scripts/16d_allfolds_crosschip_summarize.slurm; do
  [[ -s "$F" ]] || { echo "Missing $F" >&2; exit 1; }
done
A=$(sbatch --parsable scripts/16a_allfolds_crosschip_LD.slurm); A=${A%%;*}
B=$(sbatch --parsable --dependency=afterok:${A} scripts/16b_allfolds_prepare_RAW.slurm); B=${B%%;*}
C=$(sbatch --parsable --dependency=afterok:${B} scripts/16c_allfolds_crosschip_RRBLUP.slurm); C=${C%%;*}
D=$(sbatch --parsable --dependency=afterok:${C} scripts/16d_allfolds_crosschip_summarize.slurm); D=${D%%;*}
echo "Submitted: LD=${A} RAW=${B} RRBLUP=${C} SUMMARY=${D}"
echo 'All arrays run 60 tasks (10 outer × 3 inner × 2 traits).'
