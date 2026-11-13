#!/bin/bash
#SBATCH -J cf09_clump
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-59%4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09e_clump_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09e_clump_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR 09e at line ${LINENO} exit=${rc}: ${BASH_COMMAND}" >&2' ERR
WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9
PLINK_BIN="$(command -v plink)"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
TASK="${SLURM_ARRAY_TASK_ID:?Use sbatch --array=0-59}"
(( TASK>=0 && TASK<60 )) || { echo "Invalid array index $TASK" >&2; exit 2; }
OUTER="$(printf '%02d' "$((TASK / 12 + 1))")"
CHR="$((TASK % 12 + 1))"
PAD="$(printf '%02d' "$CHR")"
INPUT="$WORK/Cross_family/results/09_outer_imputed_gwas_inputs/outer_${OUTER}"
BFILE="$INPUT/imputed_train_qc"
GWAS="$WORK/Cross_family/results/09_outer_imputed_gwas/outer_${OUTER}/chr_${PAD}"
OUT="$WORK/Cross_family/results/09_outer_imputed_clumping_by_chr/outer_${OUTER}/chr_${PAD}"
R2=0.2
WINDOW_KB=1000
mkdir -p "$OUT"

for F in "$INPUT/complete.flag" "$GWAS/complete.flag" "${BFILE}.bed" "${BFILE}.bim" "${BFILE}.fam" "$GWAS/gwas_qc.tsv"; do
    [[ -s "$F" ]] || { echo "ERROR: missing/empty input: $F" >&2; exit 1; }
done
for ENTRY in 'status=PASS' 'training_only=PASS' 'outer_test_phenotype_read=NO' "outer=${OUTER}"; do
    grep -qx "$ENTRY" "$INPUT/complete.flag" || { echo "ERROR: invalid 09c flag: $ENTRY" >&2; exit 1; }
done
for ENTRY in 'status=PASS' "outer=${OUTER}" "chromosome=${CHR}" 'outer_test_phenotype_read=NO' 'traits=BB,H6'; do
    grep -qx "$ENTRY" "$GWAS/complete.flag" || { echo "ERROR: invalid 09d flag: $ENTRY" >&2; exit 1; }
done
N_TRAIN="$(wc -l < "${BFILE}.fam")"
N_EXPECT="$(awk -v chr="$CHR" '$1==chr {n++} END {print n+0}' "${BFILE}.bim")"
N_TRAIN_FLAG="$(awk -F= '$1=="n_training" {print $2}' "$INPUT/complete.flag")"
(( N_TRAIN>=100 && N_TRAIN==N_TRAIN_FLAG && N_EXPECT>=1000 )) || { echo "ERROR: unexpected training samples/SNPs" >&2; exit 1; }
# 09d ensures the GWAS result was produced with exactly the same training genotype.
for T in BB H6; do
    QC="$(awk -F '\t' -v tr="$T" 'NR>1 && $2==tr {if(++n>1)exit 3; print $4"\t"$5"\t"$6"\t"$7} END{if(n!=1)exit 4}' "$GWAS/gwas_qc.tsv")" || { echo "ERROR: missing 09d QC for $T" >&2; exit 1; }
    read -r EXPECT TESTED VALID INVALID <<< "$QC"
    (( EXPECT == N_EXPECT && TESTED == N_EXPECT && VALID == N_EXPECT && INVALID == 0 )) || {
        echo "ERROR: 09d QC mismatch for $T chr${PAD}: $QC, expected=$N_EXPECT" >&2; exit 1;
    }
done
if [[ -s "$OUT/complete.flag" ]]; then
    for ENTRY in 'status=PASS' "outer=${OUTER}" "chr=${CHR}" 'r2=0.2' 'window_kb=1000' 'training_only=PASS'; do
        grep -qx "$ENTRY" "$OUT/complete.flag" || { echo "ERROR: existing 09e flag mismatch ($ENTRY)" >&2; exit 1; }
    done
    for T in BB H6; do
        [[ -s "$OUT/${T}_chr${PAD}.index_ranked.tsv" && -s "$OUT/${T}_chr${PAD}.clumped" ]] || {
            echo "ERROR: incomplete pre-existing 09e results" >&2; exit 1;
        }
    done
    echo "[SKIP] outer_${OUTER}/chr_${PAD} 09e already PASS"
    exit 0
fi
rm -f "$OUT/complete.flag" "$OUT/clumping_summary.tsv"
printf 'Outer\tChr\tTrait\tSamples\tGWAS_SNPs\tClumped_Index_SNPs\tR2\tWindow_kb\n' > "$OUT/clumping_summary.tsv"
for TRAIT in BB H6; do
    ASSOC_SRC="$GWAS/${TRAIT}_chr${PAD}.mlma"
    PREFIX="$OUT/${TRAIT}_chr${PAD}"
    ASSOC="${PREFIX}.for_clump.tsv"
    RANK="${PREFIX}.index_ranked.tsv"
    INDEX="${PREFIX}.index_snps.txt"
    [[ -s "$ASSOC_SRC" ]] || { echo "ERROR: missing GWAS: $ASSOC_SRC" >&2; exit 1; }
    # Replicates inner 06c: lowercase GCTA p -> PLINK uppercase P, zero p -> 1e-300.
    awk -v chr="$CHR" -v expected="$N_EXPECT" '
      BEGIN {OFS="\t"}
      NR==1 {
        for(i=1;i<=NF;i++){key=tolower($i); if(key=="snp")snpcol=i; if(key=="p")pcol=i; if(key=="chr")chrcol=i}
        if(!snpcol||!pcol||!chrcol){print "ERROR: GWAS header columns missing" > "/dev/stderr";bad=1;exit 1}
        print "SNP","P"; next
      }
      {
        snp=$(snpcol); p=$(pcol)
        if($(chrcol)+0!=chr || snp=="" || seen[snp]++) {print "ERROR: wrong chromosome or duplicate/empty SNP" > "/dev/stderr";bad=1;exit 1}
        if(p !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || p+0<0 || p+0>1) {
          print "ERROR: invalid GWAS P value",p > "/dev/stderr";bad=1;exit 1
        }
        if(p+0==0) {p="1e-300";zero_p++}
        print snp,p; n++
      }
      END {
        if(bad)exit 1
        if(n!=expected) {print "ERROR: GWAS SNP count mismatch",n,expected > "/dev/stderr";exit 1}
        if(zero_p>0)print "WARNING: P-value zero floored count",zero_p > "/dev/stderr"
      }
    ' "$ASSOC_SRC" > "${ASSOC}.tmp"
    mv "${ASSOC}.tmp" "$ASSOC"
    echo "[CLUMP] outer_${OUTER}/${TRAIT} chr${CHR}: training=$N_TRAIN SNP=$N_EXPECT r2=$R2 kb=$WINDOW_KB"
    rm -f "${PREFIX}.clumped" "${PREFIX}.log"
    "$PLINK_BIN" --bfile "$BFILE" --chr "$CHR" --clump "$ASSOC" \
        --clump-snp-field SNP --clump-field P --clump-p1 1 --clump-p2 1 \
        --clump-r2 "$R2" --clump-kb "$WINDOW_KB" --threads "$CPUS" \
        --memory 24000 --out "$PREFIX"
    [[ -s "${PREFIX}.clumped" ]] || { echo "ERROR: missing PLINK clumped file" >&2; exit 1; }
    awk -v chr="$CHR" '
      BEGIN{OFS="\t"}
      NR==1 {if($1!="CHR"||$3!="SNP"||$5!="P") {print "ERROR: unexpected clumped header" > "/dev/stderr";bad=1;exit 1};print "SNP","P";next}
      NF==0{next}
      {
        if($1+0!=chr||$3=="NONE"||$3=="NA"||seen[$3]++) {print "ERROR: invalid/duplicate index SNP" > "/dev/stderr";bad=1;exit 1}
        p=$5
        if(p !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || p+0<=0 || p+0>1) {
          print "ERROR: invalid clumped P",p > "/dev/stderr";bad=1;exit 1
        }
        print $3,p;n++
      }
      END{if(bad||n<1){print "ERROR: empty/invalid clump index list" > "/dev/stderr";exit 1}}
    ' "${PREFIX}.clumped" > "${RANK}.tmp"
    { head -n 1 "${RANK}.tmp"; tail -n +2 "${RANK}.tmp" | sort -g -k2,2 -k1,1; } > "$RANK"
    rm -f "${RANK}.tmp"
    awk 'NR>1{print $1}' "$RANK" > "$INDEX"
    N_INDEX="$(wc -l < "$INDEX")"
    (( N_INDEX>=1 && N_INDEX<=N_EXPECT )) || { echo "ERROR: invalid clump index count" >&2; exit 1; }
    printf '%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\n' "$OUTER" "$CHR" "$TRAIT" "$N_TRAIN" "$N_EXPECT" "$N_INDEX" "$R2" "$WINDOW_KB" >> "$OUT/clumping_summary.tsv"
done
cat > "$OUT/complete.flag.tmp" <<EOF_FLAG
status=PASS
version=1
outer=${OUTER}
chr=${CHR}
r2=${R2}
window_kb=${WINDOW_KB}
training_only=PASS
outer_test_phenotype_read=NO
traits=BB,H6
EOF_FLAG
mv "$OUT/complete.flag.tmp" "$OUT/complete.flag"
echo "PASS outer_${OUTER}/chr_${PAD} clumping complete"
