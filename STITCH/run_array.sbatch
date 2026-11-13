#!/usr/bin/env bash
#SBATCH -p shared
#SBATCH -A naiss2025-5-732
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00

set -euo pipefail

# Resolve project root
ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MODE="${1:?Usage: sbatch run_array.sbatch position|bcf}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?Submit using --array}"

MANIFEST="$ROOT/config/manifest.tsv"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ "$MODE" == "position" || "$MODE" == "bcf" ]] || \
    fail "Unknown mode: $MODE"

[[ -r "$MANIFEST" ]] || \
    fail "Cannot read manifest: $MANIFEST"

[[ "$TASK_ID" =~ ^[1-9][0-9]*$ ]] || \
    fail "Invalid task ID: $TASK_ID"

line=$(sed -n "$((TASK_ID + 1))p" "$MANIFEST")

[[ -n "$line" ]] || \
    fail "Missing manifest row for task $TASK_ID"

IFS=$'\t' read -r CHUNK CHROM PART POS BCF <<< "$line"

[[ -n "${CHUNK:-}" && -n "${POS:-}" && -n "${BCF:-}" ]] || \
    fail "Malformed manifest row: $line"

echo "========================================"
echo "START: $(date -Is)"
echo "MODE:  $MODE"
echo "CHUNK: $CHUNK"
echo "NODE:  $(hostname)"
echo "ROOT:  $ROOT"
echo "POS:   $POS"
echo "BCF:   $BCF"
echo "========================================"

mkdir -p \
    "$ROOT/position_stats" \
    "$ROOT/per_file" \
    "$ROOT/work"

if [[ "$MODE" == "position" ]]; then

    [[ -f "$POS" ]] || \
        fail "Position file does not exist: $POS"

    [[ -r "$POS" ]] || \
        fail "Position file is not readable: $POS"

    command -v python3 >/dev/null 2>&1 || \
        fail "python3 is not available"

    python3 "$ROOT/scripts/prepare_positions.py" \
        --chunk "$CHUNK" \
        --input "$POS" \
        --keys "$ROOT/work/$CHUNK.keys.tsv" \
        --regions "$ROOT/work/$CHUNK.regions.tsv" \
        --stats "$ROOT/position_stats/$CHUNK.json"

elif [[ "$MODE" == "bcf" ]]; then

    [[ -f "$BCF" ]] || \
        fail "BCF file does not exist: $BCF"

    [[ -r "$BCF" ]] || \
        fail "BCF file is not readable: $BCF"

    [[ -s "$ROOT/work/$CHUNK.keys.tsv" ]] || \
        fail "Missing key file for $CHUNK"

    [[ -s "$ROOT/work/$CHUNK.regions.tsv" ]] || \
        fail "Missing region file for $CHUNK"

    [[ -s "$ROOT/position_stats/$CHUNK.json" ]] || \
        fail "Missing position stats for $CHUNK"

    module load bioinfo-tools
    ml bcftools/1.20

    command -v bcftools >/dev/null 2>&1 || \
        fail "bcftools is unavailable"

    bcftools --version | head -n 1

    python3 "$ROOT/scripts/match_bcf.py" \
        --chunk "$CHUNK" \
        --bcf "$BCF" \
        --keys "$ROOT/work/$CHUNK.keys.tsv" \
        --regions "$ROOT/work/$CHUNK.regions.tsv" \
        --output "$ROOT/per_file/$CHUNK.json"

fi

echo "========================================"
echo "SUCCESS: $(date -Is)"
echo "MODE: $MODE"
echo "CHUNK: $CHUNK"
echo "========================================"

