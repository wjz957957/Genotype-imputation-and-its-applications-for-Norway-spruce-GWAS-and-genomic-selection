#!/bin/bash
#SBATCH -J cf10a_model_data
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-9%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/10a_model_data_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/10a_model_data_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 10a at line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"

module load bioinfo-tools
ml plink/1.90b4.9
PLINK_BIN="$(command -v plink)"
PYTHON_BIN="$(command -v python3)"

TASK="${SLURM_ARRAY_TASK_ID:?Submit this script as an array}"
(( TASK >= 0 && TASK < 10 )) || { echo "Invalid array task ID: $TASK" >&2; exit 2; }
OUTER=$(printf '%02d' "$((TASK / 2 + 1))")
if (( TASK % 2 == 0 )); then TRAIT=BB; else TRAIT=H6; fi

SCRIPT="$WORK/Cross_family/scripts/10a_prepare_outer_model_data.py"
[[ -s "$SCRIPT" ]] || { echo "Missing script: $SCRIPT" >&2; exit 1; }

echo "===== Stage 10a outer_${OUTER}/${TRAIT} ====="
"$PYTHON_BIN" "$SCRIPT" --work "$WORK" --plink "$PLINK_BIN" \
    --outer "$OUTER" --trait "$TRAIT" --cpus "${SLURM_CPUS_PER_TASK:-4}" --memory-mb 24000

OUT="$WORK/Cross_family/results/10_outer_model_data/outer_${OUTER}/${TRAIT}"
[[ -s "$OUT/complete.flag" && -s "$OUT/final_genotype_manifest.tsv" ]] || {
    echo "ERROR: missing completion outputs under $OUT" >&2; exit 1;
}
grep -qx 'status=PASS' "$OUT/complete.flag"
echo "PASS Stage 10a outer_${OUTER}/${TRAIT}"
