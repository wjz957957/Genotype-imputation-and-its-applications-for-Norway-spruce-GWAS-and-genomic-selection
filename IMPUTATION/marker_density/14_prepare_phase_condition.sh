#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=prepare_phase
#SBATCH -t 00-04:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/prepare_%j.out
#SBATCH --error=logs/prepare_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20
ml SHAPEIT/v5.1.1

# ============================================================
# Arguments
# ============================================================

CHR2=$1
DENSITY=$2
REP=$3

CHRNUM=$((10#${CHR2}))

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
MAPDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map

# ============================================================
# Physical chromosome length
# ============================================================

case "${CHR2}" in
    01) CHR_LENGTH=1681021526 ;;
    05) CHR_LENGTH=1478570840 ;;
    06) CHR_LENGTH=1098568278 ;;
    11) CHR_LENGTH=954579047 ;;
    *)
        echo "ERROR: unsupported chromosome ${CHR2}"
        exit 1
        ;;
esac

# ============================================================
# Reproducible seed
#
# e.g.
# Chr11 / D25 / rep1 -> 1100251
# ============================================================

SEED=$((CHRNUM * 100000 + DENSITY * 10 + REP))

# Number of markers at requested SNP/Mb
N_SNPS=$(awk \
    -v bp="${CHR_LENGTH}" \
    -v d="${DENSITY}" \
    'BEGIN {printf "%.0f", bp/1000000*d}')

BASE=${WORKDIR}/Chr${CHR2}/base
OUTDIR=${WORKDIR}/Chr${CHR2}/D${DENSITY}/rep${REP}

mkdir -p "${OUTDIR}"

RAW=${BASE}/Chr${CHR2}.WGS50.raw.vcf.gz
CANDIDATES=${BASE}/Chr${CHR2}.candidate_positions.tsv

SITES=${OUTDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.sites.tsv

TARGET=${OUTDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.target.unphased.bcf

PHASED=${OUTDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.target.phased.bcf

MAP=${MAPDIR}/LD.GenMap.chr${CHRNUM}.txt

META=${OUTDIR}/Chr${CHR2}.D${DENSITY}.rep${REP}.metadata.tsv


echo "============================================================"
echo "Chr            : ${CHR2}"
echo "Density        : ${DENSITY} SNP/Mb"
echo "Replicate      : ${REP}"
echo "Physical bp    : ${CHR_LENGTH}"
echo "Number SNPs    : ${N_SNPS}"
echo "Seed           : ${SEED}"
echo "============================================================"


# ============================================================
# Check input
# ============================================================

for f in "${RAW}" "${CANDIDATES}" "${MAP}"
do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: missing ${f}"
        exit 1
    fi
done


# ============================================================
# 1. Fixed-seed random sampling
# ============================================================

python3 - \
"${CANDIDATES}" \
"${SITES}" \
"${N_SNPS}" \
"${SEED}" <<'PY'

import sys
import random

input_file = sys.argv[1]
output_file = sys.argv[2]
n = int(sys.argv[3])
seed = int(sys.argv[4])

rng = random.Random(seed)

reservoir = []

with open(input_file) as f:

    for i, line in enumerate(f):

        if i < n:
            reservoir.append(line)

        else:
            j = rng.randint(0, i)

            if j < n:
                reservoir[j] = line

if len(reservoir) != n:
    raise RuntimeError(
        f"Requested {n} SNPs but obtained {len(reservoir)}"
    )

reservoir.sort(
    key=lambda x: int(x.split()[1])
)

with open(output_file, "w") as out:
    out.writelines(reservoir)

print(f"Selected {len(reservoir)} SNPs using seed {seed}")

PY


N_SELECTED=$(wc -l < "${SITES}")

if [[ "${N_SELECTED}" -ne "${N_SNPS}" ]]; then
    echo "ERROR: selected SNP number incorrect."
    exit 1
fi


# ============================================================
# 2. Extract sparse unphased target
# ============================================================

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -R "${SITES}" \
    -Ob \
    -o "${TARGET}" \
    "${RAW}"

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${TARGET}"

TARGET_N=$(bcftools index -n "${TARGET}")
TARGET_SAMPLES=$(bcftools query -l "${TARGET}" | wc -l)

echo "Sparse target variants: ${TARGET_N}"
echo "Sparse target samples : ${TARGET_SAMPLES}"

if [[ "${TARGET_N}" -ne "${N_SNPS}" ]]; then
    echo "ERROR: sparse target marker count mismatch."
    exit 1
fi

if [[ "${TARGET_SAMPLES}" -ne 50 ]]; then
    echo "ERROR: target sample count is not 50."
    exit 1
fi


# ============================================================
# 3. SHAPEIT5
# ============================================================

SHAPEIT5_phase_common \
    -T ${SLURM_CPUS_PER_TASK} \
    -I "${TARGET}" \
    -O "${PHASED}" \
    --region "${CHRNUM}" \
    -M "${MAP}" \
    --seed "${SEED}"

bcftools index \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f "${PHASED}"


PHASED_N=$(bcftools index -n "${PHASED}")
PHASED_SAMPLES=$(bcftools query -l "${PHASED}" | wc -l)

if [[ "${PHASED_N}" -ne "${N_SNPS}" ]]; then
    echo "ERROR: phased marker count mismatch."
    exit 1
fi

if [[ "${PHASED_SAMPLES}" -ne 50 ]]; then
    echo "ERROR: phased sample count mismatch."
    exit 1
fi


# ============================================================
# Metadata
# ============================================================

{
    echo -e "Parameter\tValue"
    echo -e "Chromosome\t${CHR2}"
    echo -e "Chromosome_bp\t${CHR_LENGTH}"
    echo -e "Density_SNP_per_Mb\t${DENSITY}"
    echo -e "Target_SNPs\t${N_SNPS}"
    echo -e "Replicate\t${REP}"
    echo -e "Random_seed\t${SEED}"
    echo -e "Samples\t${PHASED_SAMPLES}"
} > "${META}"


echo
echo "Preparation + phasing completed."
echo "Finished: $(date)"
