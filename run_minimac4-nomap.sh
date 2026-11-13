#!/bin/bash

#SBATCH -A naiss2025-5-732  
#SBATCH -t 01-00:00:00   
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH -J minimac4
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL

# 加载必要的模块
module load bioinfo-tools bcftools/1.20
ml bioinfo-tools Minimac4/4.1.6

# 调用 minimac4 命令 进行vcf.gz转化msav
minimac4 --compress-reference Chr11_phased_1047-50.vcf.gz > Chr11_phased_1047-50.msav

minimac4 Chr11_phased_1047-50.msav \
  snp-array-50-phased-11-UL.vcf.gz \
  --format GT,DS,GP --min-ratio 1e-5 \
  -O vcf.gz \
  -o minimac4_output_nomap.vcf.gz 

bcftools index minimac4_output_nomap.vcf.gz
