#!/bin/bash

#SBATCH -A naiss2023-22-925
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 5
#SBATCH -J bref3
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name



# 调用 bcftools 命令 进行bcf转vcf.gz
java -jar conform-gt.24May16.cee.jar gt=snp-array-50-phased-11.vcf.gz ref=Chr11_map_phased.vcf.gz match=POS chrom=11 out=snp-array-50-phased-11-GT
