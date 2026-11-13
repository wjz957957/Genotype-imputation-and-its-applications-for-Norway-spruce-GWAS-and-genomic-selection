#!/bin/bash

#SBATCH -A naiss2023-22-925
#SBATCH -t 05-00:00:00
#SBATCH -p core -n 5
#SBATCH -J impute5-R2
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name

# Load any necessary modules
module load bioinfo-tools bcftools/1.20

# 调用 python 命令 进行R2计算
python for+R2-GT-impute5-bin-2024.07.09.py