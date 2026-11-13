#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -t 00-24:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH -J R2-global-impute5
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL



python global_R2_info.py -t Chr11_map_phased_50_rename.vcf.gz -i impute5_output_map.vcf.gz -o results.txt
