#!/usr/bin/env bash
#SBATCH -p shared
#SBATCH -A naiss2025-5-732
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
python3 "$ROOT/scripts/merge_stats.py" --root "$ROOT"
