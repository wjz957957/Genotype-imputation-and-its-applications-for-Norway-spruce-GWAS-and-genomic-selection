#!/bin/bash
#SBATCH -J cf09_chip_gwas
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-4%2
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09b2_chip_gwas_%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer/Cross_family/logs/09b2_chip_gwas_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C
trap 'rc=$?; echo "ERROR at line ${LINENO}, rc=${rc}: ${BASH_COMMAND}" >&2' ERR

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
cd "$WORK"
module load bioinfo-tools
ml plink/1.90b4.9
ml gcta/1.94.1
PLINK_BIN=$(command -v plink)
GCTA_BIN=$(command -v gcta64)
PYTHON_BIN=$(command -v python3)
VALIDATOR="$WORK/Cross_family/scripts/09b2_verify_outer_chip_selection.py"

TASK=${SLURM_ARRAY_TASK_ID:?Run as Slurm array task 0 through 4}
[[ "$TASK" =~ ^[0-4]$ ]] || { echo "Invalid task id: $TASK" >&2; exit 2; }
OUTER=$(printf '%02d' "$(( TASK + 1 ))")
N_EXPECTED_SNP=43267
THREADS=${SLURM_CPUS_PER_TASK:-4}
SCRIPT_VERSION=1
INPUT="$WORK/Cross_family/results/09_outer_inputs/outer_${OUTER}"
GENO="$WORK/Cross_family/results/09_outer_chip_genotypes/outer_${OUTER}/chip_ALL_outer_train"
OUT="$WORK/Cross_family/results/09_outer_chip_gwas_clumping/outer_${OUTER}"
GRM="$OUT/grm/chip_outer_train_grm"
GWAS_DIR="$OUT/gwas"
CLUMP_DIR="$OUT/clumping"
SET_DIR="$OUT/marker_sets"
QC="$OUT/gwas_qc.tsv"
SET_SUMMARY="$OUT/marker_set_counts.tsv"
FLAG="$OUT/complete.flag"

for F in "$VALIDATOR" "$INPUT/prepare.complete" "$INPUT/frozen_inner_selection.tsv" \
         "$INPUT/outer_train.keep" "$INPUT/outer_test.ids" \
         "$INPUT/BB.train.pheno" "$INPUT/H6.train.pheno" \
         "$INPUT/chip_ALL.snps" \
         "${GENO}.bed" "${GENO}.bim" "${GENO}.fam"; do
    [[ -s "$F" ]] || { echo "ERROR missing/empty prerequisite $F" >&2; exit 1; }
done
[[ -s "$WORK/Cross_family/results/09_outer_chip_genotypes/outer_${OUTER}/complete.flag" ]] || {
    echo "ERROR: 09b-1 not complete for outer_${OUTER}" >&2; exit 1;
}
"$PYTHON_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --phase pre
mkdir -p "$OUT" "$OUT/grm" "$GWAS_DIR" "$CLUMP_DIR" "$SET_DIR"

# Guard reuse against changed outer-training genotypes/phenotypes or frozen choices.
SIG="$OUT/input_sha256.txt"
SIG_NOW="$OUT/input_sha256.now"
sha256sum "${GENO}.bed" "${GENO}.bim" "${GENO}.fam" \
          "$INPUT/outer_train.keep" "$INPUT/chip_ALL.snps" \
          "$INPUT/BB.train.pheno" "$INPUT/H6.train.pheno" \
          "$INPUT/frozen_inner_selection.tsv" > "$SIG_NOW"
if [[ -f "$SIG" ]]; then
    cmp -s "$SIG" "$SIG_NOW" || {
        echo "ERROR: Input identity changed for outer_${OUTER}. Do not reuse the old GRM/GWAS/clumping. Inspect $OUT" >&2
        exit 1
    }
else
    # Refuse to adopt stale intermediate results with unknown input identity.
    if [[ -e "$FLAG" || -e "${GRM}.grm.bin" || -e "$GWAS_DIR/chip_BB_gwas.loco.mlma" || -e "$GWAS_DIR/chip_H6_gwas.loco.mlma" ]]; then
        echo "ERROR: existing 09b2 products lack $SIG; inspect before rerunning" >&2; exit 1
    fi
    mv "$SIG_NOW" "$SIG"
fi
rm -f "$SIG_NOW"
if [[ -s "$FLAG" ]]; then
    "$PYTHON_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --phase post
    echo "[SKIP] 09b2 outer_${OUTER} already complete and re-verified"
    exit 0
fi
rm -f "$FLAG"

N_TRAIN=$(wc -l < "${GENO}.fam")
N_SNP=$(wc -l < "${GENO}.bim")
[[ "$N_SNP" -eq "$N_EXPECTED_SNP" ]] || { echo "ERROR: expected $N_EXPECTED_SNP chip SNPs, found $N_SNP" >&2; exit 1; }

# GCTA --mlma-loco is EXACTLY the original inner Stage 03 chip GWAS method.
if [[ ! -s "${GRM}.grm.bin" || ! -s "${GRM}.grm.N.bin" || ! -s "${GRM}.grm.id" ]]; then
    "$GCTA_BIN" --bfile "$GENO" --make-grm \
        --thread-num "$THREADS" --out "$GRM"
fi
for EXT in grm.bin grm.N.bin grm.id; do
    [[ -s "${GRM}.${EXT}" ]] || { echo "ERROR missing GRM: ${GRM}.${EXT}" >&2; exit 1; }
done
[[ $(wc -l < "${GRM}.grm.id") -eq "$N_TRAIN" ]] || { echo 'ERROR: GRM count mismatch' >&2; exit 1; }
awk '{print $2}' "${GRM}.grm.id" | sort > "$OUT/grm_ids.sorted"
awk '{print $2}' "${GENO}.fam" | sort > "$OUT/train_ids.sorted"
cmp -s "$OUT/grm_ids.sorted" "$OUT/train_ids.sorted" || {
    echo 'ERROR: GRM IDs differ from outer-training FAM' >&2; exit 1;
}

printf 'Trait\tTotal_SNPs\tValid_P\tNonfinite_P\tZero_P_Floored\n' > "$QC"
printf 'Trait\tMarker_Set\tRequested_N\tActual_N\n' > "$SET_SUMMARY"

for TRAIT in BB H6; do
    echo "===== outer_${OUTER} / ${TRAIT} : outer-training-only chip GWAS ====="
    PHENO="$INPUT/${TRAIT}.train.pheno"
    PREFIX="$GWAS_DIR/chip_${TRAIT}_gwas"
    GWAS="${PREFIX}.loco.mlma"
    if [[ ! -s "$GWAS" ]]; then
        "$GCTA_BIN" --mlma-loco --bfile "$GENO" \
            --grm "$GRM" --pheno "$PHENO" \
            --thread-num "$THREADS" --out "$PREFIX"
    else
        echo "[REUSE] $GWAS"
    fi
    [[ -s "$GWAS" ]] || { echo "ERROR: GWAS output not produced: $GWAS" >&2; exit 1; }

    TC="$CLUMP_DIR/${TRAIT}"
    mkdir -p "$TC/clump"
    PFILE="$TC/gwas_p.txt"
    COUNTS="$TC/pvalue_counts.txt"
    # Kept aligned to Stage03: accepts finite [0,1] P; skips NA/nan/inf;
    # floors P=0 to 1e-300, requires all 43267 GWAS rows.
    awk -v expected="$N_EXPECTED_SNP" -v counts="$COUNTS" '
        BEGIN {OFS="\t"; numeric="^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$"; print "SNP","P"}
        NR==1 {
            if (tolower($1)!="chr" || tolower($2)!="snp" || tolower($9)!="p") {
                print "ERROR: unexpected GWAS header" > "/dev/stderr"; bad=1; exit 1
            }
            next
        }
        {
            if (NF<9 || seen[$2]++) {
                print "ERROR: malformed or duplicate GWAS SNP",$2 > "/dev/stderr"; bad=1; exit 1
            }
            total++; p=$9; low=tolower(p)
            if (low ~ /^[-+]?nan(\([^)]*\))?$/ || low ~ /^[-+]?inf(inity)?$/ ||
                low=="na" || low=="n/a" || low=="." || low=="-na") { invalid++; next }
            if (p !~ numeric || p+0<0 || p+0>1) {
                print "ERROR invalid GWAS P",$2,p > "/dev/stderr"; bad=1; exit 1
            }
            if (p+0==0) {p="1e-300"; zero++}
            print $2,p; valid++
        }
        END {
            if (bad || total!=expected || valid<1 || valid+invalid!=total) {
                print "ERROR GWAS P QC",total+0,valid+0,invalid+0,expected > "/dev/stderr"
                exit 1
            }
            printf "%d\t%d\t%d\t%d\n",total,valid,invalid,zero+0 > counts
            close(counts)
        }
    ' "$GWAS" > "$PFILE"
    read -r TOTAL VALID INVALID ZERO < "$COUNTS"
    printf '%s\t%s\t%s\t%s\t%s\n' "$TRAIT" "$TOTAL" "$VALID" "$INVALID" "$ZERO" >> "$QC"

    # EXACT original Stage03 chip clumping: r2=0.2, kb=50, P1=P2=1.
    "$PLINK_BIN" --bfile "$GENO" --clump "$PFILE" \
        --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 50 \
        --allow-extra-chr --threads "$THREADS" --memory 24000 \
        --out "$TC/clump/chip"
    CLUMPED="$TC/clump/chip.clumped"
    [[ -s "$CLUMPED" ]] || { echo "ERROR: empty clumping result: $CLUMPED" >&2; exit 1; }

    RANKED="$TC/ranked_chip_snps.tsv"
    awk '
        NR==1 {
            if (toupper($1)!="CHR" || toupper($3)!="SNP" || toupper($5)!="P") {
                print "ERROR: unexpected .clumped header" > "/dev/stderr"; bad=1; exit 1
            }
            next
        }
        NF==0 {next}
        {
            if ($3=="NONE" || $3=="NA" || seen[$3]++) {
                print "ERROR: invalid or repeated clump index",$3 > "/dev/stderr"; bad=1; exit 1
            }
            print $3 "\t" $5; n++
        }
        END {if (bad || n<1) {print "ERROR: no valid index SNPs" > "/dev/stderr"; exit 1}}
    ' "$CLUMPED" | sort -k2,2g -k1,1 > "$RANKED"
    N_RANKED=$(wc -l < "$RANKED")
    echo "[INDEX] trait=$TRAIT clumped index SNPs=$N_RANKED"

    for N in 50 100 200 300 400 500 1000 2000 3000 4000 5000 8000 10000 20000 30000 40000; do
        head -n "$N" "$RANKED" | awk '{print $1}' > "$SET_DIR/${TRAIT}_chip_${N}.snps"
        ACTUAL=$(wc -l < "$SET_DIR/${TRAIT}_chip_${N}.snps")
        printf '%s\t%s\t%s\t%s\n' "$TRAIT" "chip_${N}" "$N" "$ACTUAL" >> "$SET_SUMMARY"
    done
    awk '{print $1}' "$RANKED" > "$SET_DIR/${TRAIT}_chip_LD_ALL.snps"
    printf '%s\t%s\t%s\t%s\n' "$TRAIT" chip_LD_ALL ALL_CLUMPED "$N_RANKED" >> "$SET_SUMMARY"
    awk '{print $2}' "${GENO}.bim" > "$SET_DIR/${TRAIT}_chip_ALL.snps"
    printf '%s\t%s\t%s\t%s\n' "$TRAIT" chip_ALL "$N_EXPECTED_SNP" "$N_EXPECTED_SNP" >> "$SET_SUMMARY"

    # The previously selected chip STRATEGY is frozen, but actual outer SNPs
    # must be re-selected based on this outer-training-only GWAS.
    BEST_SET=$(awk -F '\t' -v t="$TRAIT" 'NR>1 && $2==t {print $3}' "$INPUT/frozen_inner_selection.tsv")
    [[ "$BEST_SET" =~ ^chip_[0-9]+$ ]] || { echo "ERROR: missing frozen chip strategy" >&2; exit 1; }
    BEST_N=${BEST_SET#chip_}
    BEST_FILE="$SET_DIR/${TRAIT}_${BEST_SET}.snps"
    ACTUAL=$(wc -l < "$BEST_FILE")
    [[ "$ACTUAL" -eq "$BEST_N" ]] || {
        echo "ERROR: outer_${OUTER}/${TRAIT} needs $BEST_N clumped chip SNPs, only $ACTUAL available" >&2
        exit 1
    }
    echo "[FROZEN] ${TRAIT} selects ${BEST_SET}, actual outer SNP count=${ACTUAL}"
done

# Only publish after all files and sample/SNP/QTN freeze checks pass.
"$PYTHON_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --phase post --no-flag-check
cat > "${FLAG}.tmp" <<EOF_FLAG
status=PASS
version=${SCRIPT_VERSION}
outer=${OUTER}
training_only=PASS
outer_test_phenotype_read=NO
n_train=${N_TRAIN}
chip_SNPs=${N_SNP}
chip_GWAS=GCTA_MLMA_LOCO
chip_clumping_r2=0.2
chip_clumping_kb=50
selected_outer_chip_SNPs=PASS
EOF_FLAG
mv "${FLAG}.tmp" "$FLAG"
"$PYTHON_BIN" "$VALIDATOR" --work "$WORK" --outer "$OUTER" --phase post
printf 'OVERALL PASS: 09b2 outer_%s, training=%s, chip=%s, traits=BB,H6\n' "$OUTER" "$N_TRAIN" "$N_SNP"
