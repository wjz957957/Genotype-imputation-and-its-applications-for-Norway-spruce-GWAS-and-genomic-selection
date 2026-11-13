#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=Chr11_D25_prepare
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-02:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/Chr11_D25_prepare_%j.out
#SBATCH --error=logs/Chr11_D25_prepare_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

CHR=11
CHR_LENGTH=954579047
DENSITY=25
REP=1
SEED=1001

# 25 SNP/Mb × chromosome length in Mb
N_SNPS=$(awk -v bp=${CHR_LENGTH} -v d=${DENSITY} \
    'BEGIN{printf "%.0f", d*bp/1000000}')

RAW=${WORKDIR}/Chr11/base/Chr11.WGS50.raw.vcf.gz

OUTDIR=${WORKDIR}/Chr11/D25/rep1
mkdir -p "${OUTDIR}"

ALLPOS=${OUTDIR}/Chr11.all_positions.tsv
SITES=${OUTDIR}/Chr11.D25.rep1.sites.tsv
TARGET=${OUTDIR}/Chr11.D25.rep1.target.unphased.bcf

echo "========================================"
echo "Chromosome length : ${CHR_LENGTH} bp"
echo "Density           : ${DENSITY} SNP/Mb"
echo "Number of SNPs    : ${N_SNPS}"
echo "Random seed       : ${SEED}"
echo "Input             : ${RAW}"
echo "========================================"

# -------------------------------------------------
# 1. Generate candidate site list
# -------------------------------------------------

echo "[1/4] Generating candidate position list..."

bcftools query \
    -f '%CHROM\t%POS\n' \
    "${RAW}" \
    > "${ALLPOS}"

TOTAL=$(wc -l < "${ALLPOS}")

echo "Total candidate variants: ${TOTAL}"

if [[ "${TOTAL}" -lt "${N_SNPS}" ]]; then
    echo "ERROR: fewer candidate sites than requested."
    exit 1
fi

# -------------------------------------------------
# 2. Fixed-seed reservoir sampling
# -------------------------------------------------

echo "[2/4] Randomly selecting ${N_SNPS} SNPs..."

python3 - "${ALLPOS}" "${SITES}" "${N_SNPS}" "${SEED}" <<'PY'
import sys
import random

infile = sys.argv[1]
outfile = sys.argv[2]
n = int(sys.argv[3])
seed = int(sys.argv[4])

rng = random.Random(seed)
reservoir = []

with open(infile, "r") as f:
    for i, line in enumerate(f):
        if i < n:
            reservoir.append(line)
        else:
            j = rng.randint(0, i)
            if j < n:
                reservoir[j] = line

if len(reservoir) != n:
    raise RuntimeError(
        f"Expected {n} sites but selected {len(reservoir)}"
    )

reservoir.sort(
    key=lambda x: (x.split("\t")[0], int(x.split("\t")[1]))
)

with open(outfile, "w") as out:
    out.writelines(reservoir)

print(f"Selected {len(reservoir)} SNPs with seed {seed}")
PY

N_SELECTED=$(wc -l < "${SITES}")

echo "Selected sites: ${N_SELECTED}"

if [[ "${N_SELECTED}" -ne "${N_SNPS}" ]]; then
    echo "ERROR: incorrect number of selected SNPs."
    exit 1
fi

# -------------------------------------------------
# 3. Extract sparse genotype target
# -------------------------------------------------

echo "[3/4] Extracting sparse target..."

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -R "${SITES}" \
    -Ob \
    -o "${TARGET}" \
    "${RAW}"

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${TARGET}"

# -------------------------------------------------
# 4. QC
# -------------------------------------------------

echo "[4/4] QC..."

N_TARGET_SAMPLES=$(bcftools query -l "${TARGET}" | wc -l)
N_TARGET_SNPS=$(bcftools index -n "${TARGET}")

echo
echo "========== TARGET QC =========="
echo "Samples : ${N_TARGET_SAMPLES}"
echo "SNPs    : ${N_TARGET_SNPS}"

if [[ "${N_TARGET_SAMPLES}" -ne 50 ]]; then
    echo "ERROR: target does not contain 50 samples."
    exit 1
fi

if [[ "${N_TARGET_SNPS}" -ne "${N_SNPS}" ]]; then
    echo "ERROR: target does not contain ${N_SNPS} variants."
    echo "Possible duplicated physical positions should be checked."
    exit 1
fi

echo
echo "Example target genotypes:"
bcftools query -f '[%GT\t]\n' "${TARGET}" \
    | head -5 || true

echo
echo "Preparation successfully completed."
echo "Finished: $(date)"
