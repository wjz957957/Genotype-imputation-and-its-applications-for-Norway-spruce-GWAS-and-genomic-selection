#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=build_candidates
#SBATCH --array=0-3
#SBATCH -t 00-04:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/candidates_%A_%a.out
#SBATCH --error=logs/candidates_%A_%a.err

set -euo pipefail

ml bioinfo-tools
ml bcftools/1.20

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

CHRS=(01 05 06 11)
CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}

RAW=${WORKDIR}/Chr${CHR}/base/Chr${CHR}.WGS50.raw.vcf.gz
OUT=${WORKDIR}/Chr${CHR}/base/Chr${CHR}.candidate_positions.tsv
TMP=${OUT}.tmp

echo "============================================"
echo "Chromosome : Chr${CHR}"
echo "Input      : ${RAW}"
echo "Output     : ${OUT}"
echo "Start      : $(date)"
echo "============================================"

if [[ ! -f "${RAW}" ]]; then
    echo "ERROR: missing ${RAW}"
    exit 1
fi

bcftools query \
    -f '%CHROM\t%POS\n' \
    "${RAW}" \
    > "${TMP}"

N=$(wc -l < "${TMP}")

# Check duplicate CHROM/POS
DUP=$(awk '
NR==1 {
    pc=$1
    pp=$2
    next
}
{
    if ($1==pc && $2==pp) n++
    pc=$1
    pp=$2
}
END {
    print n+0
}' "${TMP}")

echo "Candidate SNPs      : ${N}"
echo "Duplicated positions: ${DUP}"

if [[ "${N}" -eq 0 ]]; then
    echo "ERROR: no variants found."
    exit 1
fi

if [[ "${DUP}" -ne 0 ]]; then
    echo "ERROR: duplicate CHROM/POS detected."
    exit 1
fi

mv "${TMP}" "${OUT}"

echo "Finished: $(date)"
