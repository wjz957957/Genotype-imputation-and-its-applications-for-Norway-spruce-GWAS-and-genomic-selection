#!/bin/bash
#SBATCH -J 11b_BB01_COJO
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/11b_BB01_COJO_%j.out
#SBATCH -e logs/11b_BB01_COJO_%j.err

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR
export LC_ALL=C

module load bioinfo-tools
ml gcta/1.94.1

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer
INPUT="${WORK}/10_outer_test/02_imputed_inputs/outer_01"
OUT="${WORK}/10_outer_test/06_COJO_pilot/outer_01/BB"
SCRIPTS="${WORK}/scripts"

CHIP_BFILE="${INPUT}/plink/chip_train_qc"
GRM_ALL="${INPUT}/grm/chip_all"
PHENO="${INPUT}/pheno/BB.pheno"
LD_BFILE="${OUT}/selected_120_train_qc"
FINAL="${OUT}/BB_outer01_20QTN_COJO_vs_marginal.tsv"

mkdir -p "${OUT}"
for E in bed bim fam; do
    [[ -s "${CHIP_BFILE}.${E}" ]] || { echo "Missing CHIP genotype: ${CHIP_BFILE}.${E}" >&2; exit 1; }
    [[ -s "${LD_BFILE}.${E}" ]] || { echo "Missing 11a pilot genotype: ${LD_BFILE}.${E}" >&2; exit 1; }
done
for E in grm.bin grm.N.bin grm.id; do
    [[ -s "${GRM_ALL}.${E}" ]] || { echo "Missing GRM: ${GRM_ALL}.${E}" >&2; exit 1; }
done
[[ -s "${PHENO}" ]] || { echo "Missing phenotype: ${PHENO}" >&2; exit 1; }
[[ -s "${OUT}/BB_outer01_QTN_vs_CHIP_max_r2.tsv" ]] || {
    echo 'Run 11a LD pilot before 11b COJO pilot.' >&2; exit 1;
}
[[ -f "${SCRIPTS}/prepare_cojo_bb01.py" ]] || { echo 'Missing prepare_cojo_bb01.py' >&2; exit 1; }
[[ -f "${SCRIPTS}/summarize_cojo_bb01.py" ]] || { echo 'Missing summarize_cojo_bb01.py' >&2; exit 1; }

[[ $(wc -l < "${LD_BFILE}.fam") -eq 813 ]] || { echo 'Unexpected LD sample count' >&2; exit 1; }
[[ $(wc -l < "${PHENO}") -eq 813 ]] || { echo 'Unexpected phenotype sample count' >&2; exit 1; }
[[ $(wc -l < "${LD_BFILE}.bim") -eq 120 ]] || { echo 'Unexpected pilot SNP count' >&2; exit 1; }

# Use exactly the same individuals for the GWAS and COJO LD reference.
# Differences in FID/IID are fatal even if both datasets happen to have 813 rows.
if ! cmp -s \
   <(awk '{print $1 "\t" $2}' "${CHIP_BFILE}.fam" | sort) \
   <(awk '{print $1 "\t" $2}' "${LD_BFILE}.fam" | sort); then
    echo 'GWAS CHIP sample IDs and COJO LD sample IDs differ.' >&2
    exit 1
fi
if ! cmp -s \
   <(awk '{print $1 "\t" $2}' "${CHIP_BFILE}.fam" | sort) \
   <(awk '{print $1 "\t" $2}' "${PHENO}" | sort); then
    echo 'GWAS CHIP sample IDs and phenotype sample IDs differ.' >&2
    exit 1
fi
if ! cmp -s \
   <(awk '{print $1 "\t" $2}' "${CHIP_BFILE}.fam" | sort) \
   <(awk '{print $1 "\t" $2}' "${GRM_ALL}.grm.id" | sort); then
    echo 'GWAS CHIP sample IDs and GRM sample IDs differ.' >&2
    exit 1
fi

# COJO conditioning on the final 100 selected CHIP SNPs, separately by chromosome.
# Matched CHIP association stats are newly computed under exactly the GRM+LOCO
# formulation of step10b2, not taken from step10a's alternative LOCO command.
for CHR in 1 12; do
    PAD=$(printf '%02d' "$CHR")
    CHR_GRM="${INPUT}/grm/chip_chr${PAD}"
    PREFIX="${OUT}/chip_chr${PAD}_matched"
    for E in grm.bin grm.N.bin grm.id; do
        [[ -s "${CHR_GRM}.${E}" ]] || { echo "Missing GRM: ${CHR_GRM}.${E}" >&2; exit 1; }
    done
    if [[ ! -s "${PREFIX}.mlma" ]]; then
        echo "[RUN] Matched CHIP MLMA GWAS, BB outer_01 chr${CHR}"
        gcta64 \
          --mlma \
          --autosome-num 12 \
          --chr "$CHR" \
          --bfile "${CHIP_BFILE}" \
          --grm "${GRM_ALL}" \
          --mlma-subtract-grm "${CHR_GRM}" \
          --pheno "${PHENO}" \
          --thread-num "${SLURM_CPUS_PER_TASK:-4}" \
          --out "${PREFIX}"
        [[ -s "${PREFIX}.mlma" ]] || { echo "Matched chip GWAS failed for chr${CHR}" >&2; exit 1; }
    else
        echo "[REUSE] ${PREFIX}.mlma"
    fi
done

# Build a genome-wide COJO .ma from all existing BB imputed GWAS (chr1-12)
# and matched CHIP GWAS on chr1 and chr12. Do not make a top-20-only .ma:
# GCTA uses broad summary statistics when estimating the phenotypic variance.
echo '[PREP] Genome-wide summary file and condition/extract lists'
python3 "${SCRIPTS}/prepare_cojo_bb01.py" \
    --work "${WORK}" --out "${OUT}"
MA="${OUT}/BB_outer01_genomewide_with_matched_chip.ma"
[[ -s "$MA" && -s "${OUT}/cojo_input_ready.flag" ]] || exit 1

for CHR in 1 12; do
    echo "[COJO] BB outer_01 chr${CHR}"
    PREFIX="${OUT}/cojo_chr${CHR}"
    # With 1,000,000 kb, conditioning also accounts for far-apart CHIP variants
    # on the *same* Norway spruce chromosome, instead of assuming no LD at 10 Mb.
    # All 100 CHIP SNPs are included across the two chromosome-specific runs;
    # CHIP SNPs on other chromosomes are not used as conditioning covariates.
    gcta64 \
      --bfile "${LD_BFILE}" \
      --autosome-num 12 \
      --chr "$CHR" \
      --extract "${OUT}/chr${CHR}_selected.snps" \
      --cojo-file "$MA" \
      --cojo-cond "${OUT}/chr${CHR}_conditional_chip.snps" \
      --cojo-wind 1000000 \
      --thread-num "${SLURM_CPUS_PER_TASK:-4}" \
      --out "$PREFIX"

    if [[ ! -s "${PREFIX}.cma.cojo" && ! -s "${PREFIX}.cma" ]]; then
        echo "Missing COJO cma result for chr${CHR}; inspect ${PREFIX}.log" >&2
        exit 1
    fi
done

python3 "${SCRIPTS}/summarize_cojo_bb01.py" --work "${WORK}" --out "${OUT}"
[[ $(wc -l < "${FINAL}") -eq 21 ]] || { echo "Expected 21 report lines: ${FINAL}" >&2; exit 1; }
echo "COJO pilot finished: ${FINAL}"
