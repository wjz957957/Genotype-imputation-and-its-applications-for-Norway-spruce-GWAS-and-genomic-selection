#!/bin/bash

#SBATCH -A naiss2023-22-925
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 2
#SBATCH -J minimac4
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name

# Load any necessary modules
module load bioinfo-tools bcftools/1.20

# 调用 minimac4 命令 进行填充
# 执行minimac4命令，如果成功则执行bcftools index
/crex/proj/snic2020-6-176/nobackup/private/jingz/minimac4/bin/minimac4 Chr11_map_phased_1047-50.msav snp-array-50-phased-11-UL.vcf.gz --format GT,DS,GP --min-ratio 1e-5 -m LD.GenMap.chr11_plink.txt -O vcf.gz -o minimac4_output_map.vcf.gz && bcftools index minimac4_output_map.vcf.gz
