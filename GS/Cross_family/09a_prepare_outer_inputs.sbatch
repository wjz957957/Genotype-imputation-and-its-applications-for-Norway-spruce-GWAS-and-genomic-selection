#!/bin/bash
#SBATCH -J cf09_outer_inputs
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 00:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-4%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09a_outer_inputs_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09a_outer_inputs_%A_%a.err
set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 09a at line ${LINENO}, code=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
TASK=${SLURM_ARRAY_TASK_ID:?must run as sbatch array}
((TASK >= 0 && TASK < 5)) || { echo "invalid array task id: $TASK" >&2; exit 1; }
OUTER=$(printf '%02d' "$((TASK+1))")
PYTHON=$(command -v python3)
SCRIPT="$WORK/Cross_family/scripts/09a_prepare_outer_inputs.py"
[[ -s "$SCRIPT" ]] || { echo "missing $SCRIPT" >&2; exit 1; }
OUTPUT="$WORK/Cross_family/results/09_outer_inputs/outer_${OUTER}"
mkdir -p "$OUTPUT"
echo "Running 09a outer_${OUTER} using frozen signed-r QTN decisions"
"$PYTHON" "$SCRIPT" --work "$WORK" --outer "$OUTER" --out "$OUTPUT"
for F in "$OUTPUT/prepare.complete" "$OUTPUT/frozen_inner_selection.tsv" "$OUTPUT/chip_ALL.snps" "$OUTPUT/BB.train.pheno" "$OUTPUT/H6.train.pheno" "$OUTPUT/outer_train.keep" "$OUTPUT/outer_test.ids"; do
    [[ -s "$F" ]] || { echo "missing/empty output: $F" >&2; exit 1; }
done
grep -qx 'status=PASS' "$OUTPUT/prepare.complete"
grep -qx 'outer_test_phenotype_read=NO' "$OUTPUT/prepare.complete"
echo "PASS stage09a outer_${OUTER}"
