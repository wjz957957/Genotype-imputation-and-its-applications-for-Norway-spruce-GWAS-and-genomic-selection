#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="${ROOT:-$(dirname "$SCRIPT_DIR")}";export ROOT
POSROOT=/cfs/klemming/projects/supr/uppstore2017145/V2/users/chen/GWAS
BCFROOT=/cfs/klemming/projects/supr/uppstore2017145/V2/data/variation/P.abies/bcf
python3 "$ROOT/scripts/build_manifest.py" --root "$ROOT" --pos-root "$POSROOT" --bcf-root "$BCFROOT"
N=$(( $(wc -l < "$ROOT/config/manifest.tsv") - 1 ))
[[ "$N" -gt 0 ]] || exit 1
mkdir -p "$ROOT"/{position_stats,per_file,summary,work,logs}
if [[ "${RUN:-0}" != 1 ]];then
  echo "DRY RUN: $N paired inputs. No tasks submitted. RUN=1 bash $0"
  head -n 3 "$ROOT/config/manifest.tsv"
  exit 0
fi
if find "$ROOT/position_stats" "$ROOT/per_file" "$ROOT/work" -type f | grep -q .;then
 echo 'Existing intermediate outputs found. For reproducibility use a fresh root or manually archive/remove old intermediates.' >&2;exit 1
fi
P=$(sbatch --parsable --array="1-${N}%12" --export=ALL \
  --output="$ROOT/logs/position_%A_%a.out" --error="$ROOT/logs/position_%A_%a.err" "$ROOT/scripts/run_array.sbatch" position)
B=$(sbatch --parsable --dependency="afterok:$P" --array="1-${N}%2" --export=ALL \
  --output="$ROOT/logs/bcf_%A_%a.out" --error="$ROOT/logs/bcf_%A_%a.err" "$ROOT/scripts/run_array.sbatch" bcf)
M=$(sbatch --parsable --dependency="afterok:$B" --export=ALL \
  --output="$ROOT/logs/merge_%j.out" --error="$ROOT/logs/merge_%j.err" "$ROOT/scripts/merge.sbatch")
printf 'Submitted position=%s BCF=%s merge=%s chunks=%s\n' "$P" "$B" "$M" "$N"
