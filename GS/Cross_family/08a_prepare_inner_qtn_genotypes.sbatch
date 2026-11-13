#!/bin/bash
#SBATCH -J cf_qtn_raw
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-29%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/08a_qtn_raw_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/08a_qtn_raw_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9
PLINK_BIN="$(command -v plink)"
PYTHON_BIN="$(command -v python3)"
HELPER="$WORK/Cross_family/scripts/08a_prepare_qtn_rrblup_inputs.py"
MASTER="$WORK/08_imputed_GWAS/00_reference/all_904_genotypes"
OUT_ROOT="$WORK/Cross_family/results/08_inner_qtn_genotypes"
TASK="${SLURM_ARRAY_TASK_ID:?task required}"
[[ "$TASK" =~ ^[0-9]+$ ]] && (( TASK >= 0 && TASK < 30 )) || { echo "bad task ID: $TASK" >&2; exit 2; }
OUTER=$(printf '%02d' "$((TASK/6+1))")
INNER=$(printf '%02d' "$((TASK%6/2+1))")
if (( TASK % 2 == 0 )); then TRAIT=BB; else TRAIT=H6; fi
OUT="$OUT_ROOT/outer_${OUTER}/inner_${INNER}/$TRAIT"
mkdir -p "$OUT"
if [[ -s "$OUT/complete.flag" ]]; then
    grep -qx 'status=PASS' "$OUT/complete.flag" || { echo 'ERROR: invalid existing completion flag' >&2; exit 1; }
    for F in train.raw valid.raw selected_chip.snps extra_700.snps; do
        [[ -s "$OUT/$F" ]] || { echo "ERROR: incomplete flagged result $OUT/$F" >&2; exit 1; }
    done
    "$PYTHON_BIN" "$HELPER" verify --work "$WORK" --outer "$OUTER" --inner "$INNER" --trait "$TRAIT" --out "$OUT"
    echo "SKIP: verified existing outer_${OUTER}/inner_${INNER}/$TRAIT"
    exit 0
fi
rm -f "$OUT/complete.flag"
"$PYTHON_BIN" "$HELPER" prepare --work "$WORK" --outer "$OUTER" --inner "$INNER" --trait "$TRAIT" --out "$OUT"
for PART in train valid; do
    "$PLINK_BIN" --bfile "$MASTER" --keep "$OUT/$PART.keep" \
        --extract "$OUT/union.snps" --keep-allele-order --recode A \
        --threads "${SLURM_CPUS_PER_TASK:-4}" --memory 36000 \
        --out "$OUT/$PART"
    [[ -s "$OUT/$PART.raw" ]] || { echo "ERROR: missing $PART.raw" >&2; exit 1; }
done
"$PYTHON_BIN" "$HELPER" verify --work "$WORK" --outer "$OUTER" --inner "$INNER" --trait "$TRAIT" --out "$OUT"
echo "PASS: stage08a outer_${OUTER}/inner_${INNER}/$TRAIT"
