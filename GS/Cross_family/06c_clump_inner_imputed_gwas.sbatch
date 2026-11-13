#!/bin/bash
#SBATCH -J cf_qtn_clump
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-179%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06c_qtn_clump_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/06c_qtn_clump_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer"
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9
PLINK_BIN="$(command -v plink)"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
INPUT_ROOT="$WORK/Cross_family/results/06_inner_imputed_gwas_inputs"
GWAS_ROOT="$WORK/Cross_family/results/06_inner_imputed_gwas"
OUT_ROOT="$WORK/Cross_family/results/06_inner_imputed_clumping_by_chr"
R2=0.2
WINDOW_KB=1000
TASK="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
(( TASK >= 0 && TASK < 180 )) || { echo "Invalid task id: $TASK" >&2; exit 2; }
FOLD_IDX=$((TASK / 12))
CHR=$((TASK % 12 + 1))
OUTER=$(printf '%02d' "$((FOLD_IDX / 3 + 1))")
INNER=$(printf '%02d' "$((FOLD_IDX % 3 + 1))")
PAD=$(printf '%02d' "$CHR")
FOLD="outer_${OUTER}/inner_${INNER}"
BFILE="$INPUT_ROOT/$FOLD/imputed_train_qc"
GWAS_DIR="$GWAS_ROOT/$FOLD/chr_${PAD}"
OUT="$OUT_ROOT/$FOLD/chr_${PAD}"
mkdir -p "$OUT"
[[ -s "$GWAS_DIR/complete.flag" ]] && grep -qx 'status=PASS' "$GWAS_DIR/complete.flag" || { echo "ERROR: GWAS chromosome stage incomplete: $GWAS_DIR/complete.flag" >&2; exit 1; }
for ext in bed bim fam; do [[ -s "${BFILE}.${ext}" ]] || { echo "ERROR: missing ${BFILE}.${ext}" >&2; exit 1; }; done
N_SAMPLES=$(wc -l < "${BFILE}.fam")
N_EXPECT=$(awk -v chr="$CHR" '$1==chr{n++} END{print n+0}' "${BFILE}.bim")
(( N_SAMPLES >= 100 && N_EXPECT >= 1000 )) || { echo "ERROR: unexpected sample/SNP count" >&2; exit 1; }
rm -f "$OUT/complete.flag" "$OUT/clumping_summary.tsv"
printf 'Outer\tInner\tChr\tTrait\tSamples\tGWAS_SNPs\tClumped_Index_SNPs\tR2\tWindow_kb\n' > "$OUT/clumping_summary.tsv"
for TRAIT in BB H6; do
  GWAS="${GWAS_DIR}/${TRAIT}_chr${PAD}.mlma"
  [[ -s "$GWAS" ]] || { echo "ERROR: missing GWAS $GWAS" >&2; exit 1; }
  PREFIX="$OUT/${TRAIT}_chr${PAD}"
  ASSOC="${PREFIX}.for_clump.tsv"
  RANKED="${PREFIX}.index_ranked.tsv"
  INDEX_LIST="${PREFIX}.index_snps.txt"
  awk -v chr="$CHR" -v expected="$N_EXPECT" '
    BEGIN{OFS="\t"}
    NR==1 {for(i=1;i<=NF;i++){field=tolower($i);if(field=="snp")snpcol=i;if(field=="p")pcol=i;if(field=="chr")chrcol=i}; if(!snpcol||!pcol||!chrcol){print "ERROR: invalid GWAS header" > "/dev/stderr";bad=1;exit 1}; print "SNP","P";next}
    {snp=$(snpcol);p=$(pcol);if($(chrcol)+0!=chr){print "ERROR: wrong chromosome" > "/dev/stderr";bad=1;exit 1};if(seen[snp]++){print "ERROR: duplicate SNP" > "/dev/stderr";bad=1;exit 1};if(p !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || p+0<0 || p+0>1){print "ERROR: invalid P-value",p,"for",snp > "/dev/stderr";bad=1;exit 1};if(p+0==0){p="1e-300";zero_p++};print snp,p;n++}
    END{if(bad)exit 1;if(n!=expected){print "ERROR: GWAS SNP count mismatch",n,expected > "/dev/stderr";exit 1};if(zero_p>0)print "WARNING: P-values floored to 1e-300:",zero_p > "/dev/stderr"}
  ' "$GWAS" > "$ASSOC"
  echo "[CLUMP] outer_${OUTER}/inner_${INNER}/${TRAIT} chr${CHR}: r2=${R2}, window=${WINDOW_KB} kb"
  "$PLINK_BIN" --bfile "$BFILE" --chr "$CHR" --clump "$ASSOC" \
    --clump-snp-field SNP --clump-field P --clump-p1 1 --clump-p2 1 \
    --clump-r2 "$R2" --clump-kb "$WINDOW_KB" --threads "$CPUS" \
    --memory 24000 --out "$PREFIX"
  [[ -s "${PREFIX}.clumped" ]] || { echo "ERROR: missing clumped output ${PREFIX}.clumped" >&2; exit 1; }
  awk -v chr="$CHR" '
    BEGIN{OFS="\t"}
    NR==1 {if($1!="CHR"||$3!="SNP"||$5!="P"){print "ERROR: unexpected .clumped header" > "/dev/stderr";bad=1;exit 1};print "SNP","P";next}
    NF==0{next}
    {if($1+0!=chr||$3=="NONE"||$3=="NA"){print "ERROR: invalid clump index SNP" > "/dev/stderr";bad=1;exit 1};if(seen[$3]++){print "ERROR: duplicate clump index SNP" > "/dev/stderr";bad=1;exit 1};print $3,$5;n++}
    END{if(bad||n<1){print "ERROR: no valid index SNPs" > "/dev/stderr";exit 1}}
  ' "${PREFIX}.clumped" > "$RANKED"
  { head -n 1 "$RANKED"; tail -n +2 "$RANKED" | sort -g -k2,2 -k1,1; } > "${RANKED}.tmp"
  mv "${RANKED}.tmp" "$RANKED"
  awk 'NR>1{print $1}' "$RANKED" > "$INDEX_LIST"
  N_INDEX=$(wc -l < "$INDEX_LIST")
  (( N_INDEX >= 1 && N_INDEX <= N_EXPECT )) || { echo "ERROR: invalid index SNP count $N_INDEX" >&2; exit 1; }
  printf '%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n' "$OUTER" "$INNER" "$CHR" "$TRAIT" "$N_SAMPLES" "$N_EXPECT" "$N_INDEX" "$R2" "$WINDOW_KB" >> "$OUT/clumping_summary.tsv"
done
printf 'status=PASS\nouter=%s\ninner=%s\nchr=%s\nr2=%s\nwindow_kb=%s\n' "$OUTER" "$INNER" "$CHR" "$R2" "$WINDOW_KB" > "$OUT/complete.flag"
echo "PASS clumping outer_${OUTER}/inner_${INNER}/chr_${PAD}"
