#!/bin/bash

#SBATCH -A naiss2023-5-58
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 2
#SBATCH -J GWAS
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name


# Load any necessary modules
module load bioinfo-tools gcta/1.94.1

# 调用 bcftools 命令 
gcta64 --bfile $1-$2-plink --grm-sparse $1-$2-plink-G --fastGWA-mlm --pheno Traits_$1.txt --qcovar $1-$2-plink-pca5.eigenvec  --thread-num 4 --out $1-$2-PC5-GWAS
