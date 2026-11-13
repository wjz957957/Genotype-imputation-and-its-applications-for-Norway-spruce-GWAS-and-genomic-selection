#!/bin/bash
#SBATCH -J cf09_merge
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-4%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09f_merge_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09f_merge_%A_%a.err
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR 09f merge at line ${LINENO} exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
TASK="${SLURM_ARRAY_TASK_ID:?Run via sbatch}"
((TASK>=0 && TASK<5)) || exit 2
OUTER="$(printf '%02d' "$((TASK+1))")"
python3 "$WORK/Cross_family/scripts/09f_merge_outer_imputed_clumps.py" --work "$WORK" --outer "$OUTER"
