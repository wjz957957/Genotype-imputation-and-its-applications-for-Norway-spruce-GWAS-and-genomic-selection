#!/bin/bash

#SBATCH -A naiss2023-22-925
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 2
#SBATCH -J conform-GT
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name

# Load any necessary modules
module load bioinfo-tools bcftools/1.20

# 调用 conformgt 命令 进行校正并添加ACAN字段
java -jar /crex/proj/snic2020-6-176/nobackup/private/jingz/conform-gt.24May16.cee.jar gt=exome-11-1303-rename.vcf.gz ref=Chr11_map_phased.vcf.gz match=POS chrom=11 out=exome-11-1302-GT

bcftools index exome-11-1302-GT.vcf.gz

bcftools +fill-tags exome-11-1302-GT.vcf.gz -Oz -o exome-11-1302-UL.vcf.gz -- -t AC,AN

bcftools index exome-11-1302-UL.vcf.gz
