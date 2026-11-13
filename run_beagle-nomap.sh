#!/bin/bash

#SBATCH -A naiss2025-5-732  
#SBATCH -t 01-00:00:00   
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH -J beagle
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL

# 加载必要的模块
module load bioinfo-tools bcftools/1.20

java -jar -Xmx40G /cfs/klemming/projects/supr/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/beagle.22Jul22.46e.jar \
  gt=snp-array-50-phased-11-UL.vcf.gz \
  gp=true \
  ref=Chr11_phased_1047-50.vcf.gz \
  overlap=5 \
  out=beagle_output_nomap
  
bcftools index beagle_output_nomap.vcf.gz
