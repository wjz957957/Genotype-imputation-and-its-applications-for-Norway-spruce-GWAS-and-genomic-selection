#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH --job-name=prep_rand
#SBATCH -t 00-04:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/prep_rand_%j.out
#SBATCH --error=logs/prep_rand_%j.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20
ml SHAPEIT/v5.1.1

CHR2=$1
REP=$2
CHRNUM=$((10#${CHR2}))

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review
MAPDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/Chr-phased-map

QC=${WORKDIR}/Chr${CHR2}/array_real/Chr${CHR2}.array_conform_QC.txt

N_SNPS=$(awk '$1=="Usable_array_markers_after_conform"{print $2}' "${QC}")

if [[ -z "${N_SNPS}" ]]; then
    echo "ERROR: cannot read usable array marker count from ${QC}"
    exit 1
fi

SEED=$((CHRNUM * 1000000 + N_SNPS * 10 + REP))

BASE=${WORKDIR}/Chr${CHR2}/base
OUTDIR=${WORKDIR}/Chr${CHR2}/random_matched/rep${REP}

mkdir -p "${OUTDIR}"
mkdir -p "${WORKDIR}/logs"

RAW=${BASE}/Chr${CHR2}.WGS50.raw.vcf.gz
CANDIDATES=${BASE}/Chr${CHR2}.candidate_positions.tsv
MAP=${MAPDIR}/LD.GenMap.chr${CHRNUM}.txt

SITES=${OUTDIR}/Chr${CHR2}.randomN.rep${REP}.sites.tsv
TARGET=${OUTDIR}/Chr${CHR2}.randomN.rep${REP}.target.unphased.bcf
PHASED=${OUTDIR}/Chr${CHR2}.randomN.rep${REP}.target.phased.bcf
META=${OUTDIR}/Chr${CHR2}.randomN.rep${REP}.metadata.tsv

echo "Chr=${CHR2}"
echo "Rep=${REP}"
echo "Matched N=${N_SNPS}"
echo "Seed=${SEED}"

python3 - "${CANDIDATES}" "${SITES}" "${N_SNPS}" "${SEED}" <<'PY'
import sys, random

inp, outp, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

rng = random.Random(seed)
reservoir = []

with open(inp) as f:
    for i, line in enumerate(f):
        if i < n:
            reservoir.append(line)
        else:
            j = rng.randint(0, i)
            if j < n:
                reservoir[j] = line

if len(reservoir) != n:
    raise RuntimeError(f"Requested {n}, obtained {len(reservoir)}")

reservoir.sort(key=lambda x: int(x.split()[1]))

with open(outp, "w") as out:
    out.writelines(reservoir)
PY

[[ $(wc -l < "${SITES}") -eq "${N_SNPS}" ]] || {
    echo "ERROR: selected marker count mismatch"
    exit 1
}

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -R "${SITES}" \
    -Ob \
    -o "${TARGET}" \
    "${RAW}"

bcftools index -f "${TARGET}"

TARGET_N=$(bcftools index -n "${TARGET}")
TARGET_NS=$(bcftools query -l "${TARGET}" | wc -l)

[[ "${TARGET_N}" -eq "${N_SNPS}" ]] || {
    echo "ERROR: target marker count mismatch"
    exit 1
}

[[ "${TARGET_NS}" -eq 50 ]] || {
    echo "ERROR: target sample count mismatch"
    exit 1
}

SHAPEIT5_phase_common \
    -T "${SLURM_CPUS_PER_TASK}" \
    -I "${TARGET}" \
    -O "${PHASED}" \
    --region "${CHRNUM}" \
    -M "${MAP}" \
    --seed "${SEED}"

bcftools index -f "${PHASED}"

PHASED_N=$(bcftools index -n "${PHASED}")

[[ "${PHASED_N}" -eq "${N_SNPS}" ]] || {
    echo "ERROR: phased marker count mismatch"
    exit 1
}

{
    echo -e "Parameter\tValue"
    echo -e "Chromosome\t${CHR2}"
    echo -e "Matched_actual_array_N\t${N_SNPS}"
    echo -e "Replicate\t${REP}"
    echo -e "Random_seed\t${SEED}"
} > "${META}"

echo "Finished: $(date)"
