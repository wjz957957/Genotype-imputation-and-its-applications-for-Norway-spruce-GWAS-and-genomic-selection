#!/bin/bash

#SBATCH -A naiss2024-5-699
#SBATCH -t 00-05:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH -J GWAS
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH --mail-type=ALL


# Load any necessary modules
module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11
module load bioinfo-tools bcftools/1.20
module load bioinfo-tools vcftools/0.1.16
module load bioinfo-tools plink/1.90b4.9
module load bioinfo-tools gcta/1.94.1


# 调用 bcftools 命令
plink --vcf SSF.exome.allchr-bud-filter.vcf.gz --make-bed --out Bud-origin-plink --threads 10 --double-id

