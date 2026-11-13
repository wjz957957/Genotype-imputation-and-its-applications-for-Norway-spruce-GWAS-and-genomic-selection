#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH --job-name=Chr11_D25_eval
#SBATCH --mail-type=END,FAIL
#SBATCH -t 00-12:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --output=logs/Chr11_D25_eval_%j.out
#SBATCH --error=logs/Chr11_D25_eval_%j.err

set -euo pipefail

WORKDIR=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/SNP-density-review

cd "${WORKDIR}"

TRUTH=${WORKDIR}/Chr11/base/Chr11.WGS50.raw.vcf.gz

IMPUTED=${WORKDIR}/Chr11/D25/rep1/Chr11.D25.rep1.imputed.vcf.gz

OUTPUT=${WORKDIR}/Chr11/D25/rep1/Chr11.D25.rep1.accuracy.tsv

echo "======================================"
echo "Truth   : ${TRUTH}"
echo "Imputed : ${IMPUTED}"
echo "Output  : ${OUTPUT}"
echo "Start   : $(date)"
echo "======================================"

python3 ${WORKDIR}/08_evaluate_imputation.py \
    --truth "${TRUTH}" \
    --imputed "${IMPUTED}" \
    --output "${OUTPUT}" \
    --progress 100000

echo
echo "Finished: $(date)"
