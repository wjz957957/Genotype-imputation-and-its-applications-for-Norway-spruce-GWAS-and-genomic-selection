#!/bin/bash

#SBATCH -A naiss2025-5-732
#SBATCH -t 00-08:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -J combine
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH --mail-type=ALL


# Load any necessary modules
module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03


# 调用 bcftools 命令
# Rscript plot_GWAS_manhattan_dataset_specific_thresholds.R
Rscript plot_GWAS_manhattan_dataset_specific_thresholds-new.R
