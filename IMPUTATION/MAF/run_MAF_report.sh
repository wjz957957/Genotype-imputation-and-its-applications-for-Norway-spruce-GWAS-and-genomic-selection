#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -t 00-24:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH -J MAF_report
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL


python MAF_stats.py -i Chr11_MAF-347_with_MAF.vcf.gz -o Chr11_MAF-347_with_MAF_report.txt
