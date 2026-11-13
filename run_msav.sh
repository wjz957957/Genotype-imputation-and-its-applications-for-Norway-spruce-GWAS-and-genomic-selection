#!/bin/bash

#SBATCH -A naiss2023-22-925
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 5
#SBATCH -J msav
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name

# Define PATH
cd /crex/proj/snic2020-6-176/nobackup/private/jingz/Chr-phased-map/


# 调用 minimac4 命令 进行vcf.gz转化msav
/crex/proj/snic2020-6-176/nobackup/private/jingz/minimac4/bin/minimac4 --compress-reference Chr${1}_map_phased.vcf.gz > Chr${1}_map_phased.msav
