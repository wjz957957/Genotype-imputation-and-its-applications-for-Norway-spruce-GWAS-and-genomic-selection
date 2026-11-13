#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -t 00-24:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH -J phasing
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL

module load bioinfo-tools SHAPEIT/v5.1.1

SHAPEIT5_phase_common -T 10 \
	-I Chr11_1047-50.vcf.gz \
	-O Chr11_phased_1047-50.bcf \
	--region 11
