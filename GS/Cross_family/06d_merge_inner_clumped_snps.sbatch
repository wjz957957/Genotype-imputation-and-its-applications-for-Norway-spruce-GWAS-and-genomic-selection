#!/bin/bash
#SBATCH -J cf_qtn_merge
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-14%5
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06d_qtn_merge_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06d_qtn_merge_%A_%a.err
set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
(( TASK >= 0 && TASK < 15 )) || { echo "Invalid task id: $TASK" >&2; exit 2; }
OUTER=$(printf '%02d' "$((TASK / 3 + 1))")
INNER=$(printf '%02d' "$((TASK % 3 + 1))")
python3 "$WORK/Cross_family/scripts/06d_merge_inner_clumped_snps.py" --work "$WORK" --outer "$OUTER" --inner "$INNER"
