#!/bin/bash

#SBATCH -A naiss2023-22-925
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 10
#SBATCH -J beagle
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name

# Load any necessary modules
module load bioinfo-tools bcftools/1.20

# 调用 beagle 命令 进行填充
# 执行minimac4命令，如果成功则执行bcftools index
/crex/proj/snic2020-6-176/nobackup/private/jingz/jdk-22/bin/java -jar /crex/proj/snic2020-6-176/nobackup/private/jingz/beagle.27May24.118.jar \
  gt=snp-array-50-phased-11-UL.vcf.gz \
  gp=true \
  ref=Chr11_map_phased_1047-50.bref3 \
  map=LD.GenMap.chr11_plink.txt \
  overlap=5 \
  out=beagle_output_map && \
bcftools index beagle_output_map.vcf.gz
