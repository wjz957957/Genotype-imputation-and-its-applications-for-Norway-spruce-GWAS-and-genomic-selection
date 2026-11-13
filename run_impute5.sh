#!/bin/bash

#SBATCH -A naiss2023-22-925
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 10
#SBATCH -J impute5
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name

# Load any necessary modules
module load bioinfo-tools bcftools/1.20

# 调用 impute5 命令 进行填充
# 执行impute5命令，如果成功则执行bcftools index
/crex/proj/snic2020-6-176/nobackup/private/jingz/impute5_v1.2.1_x86_64_static \
  --h Chr11_map_phased_1047-50.vcf.gz \
  --g snp-array-50-phased-11-UL.vcf.gz \
  --m LD.GenMap.chr11.txt \
  --o impute5_output_map.vcf.gz \
  --r 11:42607-952042565 \
  --buffer-region 11:1-954579047 && \
bcftools index impute5_output_map.vcf.gz
