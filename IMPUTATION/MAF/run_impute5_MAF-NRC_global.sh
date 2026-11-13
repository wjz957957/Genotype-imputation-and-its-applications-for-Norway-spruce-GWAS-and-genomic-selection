#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -t 00-10:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -J NRC-global-MAF
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL

python3 maf_NRC_global.py \
  -t Chr11_MAF-347_with_MAF.vcf.gz \
  -i impute5-output-347-maf-ac-an.vcf.gz \
  -o results_MAF-NRC-all.txt \
  --matrix-output results_MAF-NRC-matrices.txt \
  --info-threshold 0.8
