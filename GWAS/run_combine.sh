#!/bin/bash

#SBATCH -A naiss2024-5-699
#SBATCH -t 00-08:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -J combine
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH --mail-type=ALL


# Load any necessary modules
module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11
module load bioinfo-tools bcftools/1.20
module load bioinfo-tools vcftools/0.1.16
module load bioinfo-tools plink/1.90b4.9
module load bioinfo-tools gcta/1.94.1


# 调用 bcftools 命令
cat Bud-imputed-PC0-GWASloco_chr1.mlma | while read line; do echo "$line"; done > Bud-imputed-PC0-GWASloco_chrall.mlma
for file in Bud-imputed-PC0-GWASloco_chr2.mlma Bud-imputed-PC0-GWASloco_chr3.mlma Bud-imputed-PC0-GWASloco_chr4.mlma Bud-imputed-PC0-GWASloco_chr5.mlma Bud-imputed-PC0-GWASloco_chr6.mlma Bud-imputed-PC0-GWASloco_chr7.mlma Bud-imputed-PC0-GWASloco_chr8.mlma Bud-imputed-PC0-GWASloco_chr9.mlma Bud-imputed-PC0-GWASloco_chr10.mlma Bud-imputed-PC0-GWASloco_chr11.mlma Bud-imputed-PC0-GWASloco_chr12.mlma; do awk 'NR==1{if(NR==FNR){next}}1' "$file" >> Bud-imputed-PC0-GWASloco_chrall.mlma; done

