#!/usr/bin/env bash
#SBATCH -p shared
#SBATCH -A naiss2025-5-732
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00
set -euo pipefail
ROOT="${ROOT:?ROOT must be provided}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?Submit with --array}"
MANIFEST="$ROOT/config/manifest.tsv"
line=$(sed -n "$((TASK_ID + 1))p" "$MANIFEST")
IFS=$'\t' read -r CHUNK CHROM PART POS BCF <<< "$line"
CHUNK=${CHUNK%$'\r'}; POS=${POS%$'\r'}; BCF=${BCF%$'\r'}
echo "START $(date -Is) chunk=$CHUNK node=$(hostname)"
[[ -r "$BCF" ]] || { echo "Cannot read BCF: $BCF" >&2; exit 1; }
[[ -s "$ROOT/work/$CHUNK.keys.tsv" ]] || { echo "Missing keys for $CHUNK" >&2; exit 1; }
module load bioinfo-tools
ml bcftools/1.20
bcftools --version | head -n 1
TARGETS="$ROOT/work_v4/$CHUNK.targets.tsv"
mkdir -p "$ROOT/work_v4" "$ROOT/per_file_v4" "$ROOT/logs"
python3 "$ROOT/scripts/make_targets.py" --keys "$ROOT/work/$CHUNK.keys.tsv" --targets "$TARGETS"
python3 "$ROOT/scripts/match_bcf_stream.py" \
  --chunk "$CHUNK" --bcf "$BCF" --keys "$ROOT/work/$CHUNK.keys.tsv" \
  --targets "$TARGETS" --output "$ROOT/per_file_v4/$CHUNK.json" \
  --progress-every 1000000
echo "SUCCESS $(date -Is) chunk=$CHUNK"
