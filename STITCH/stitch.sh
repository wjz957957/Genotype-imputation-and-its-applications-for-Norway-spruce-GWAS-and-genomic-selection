#!/bin/bash

#SBATCH -A snic2022-5-632
#SBATCH -p core -n 10
#SBATCH -t 15:00:00
#SBATCH --mail-user zhiqiang.chen@slu.se

# get tabuate postion information
# bcftools query -f '%CHROM \t %POS\t %REF\t %ALT\n' PA_chr01_1.indel.bia.Qual20.DP.mac5.miss0.3.HWE.recode.vcf.gz > pos.tabulated.txt


/crex/proj/uppstore2017145/V2/users/chen/bin/STITCH/STITCH.R --chr=PA_chr01_1 \
# --regionStart=5000 \
# --regionEnd=50000 \
# --buffer=5000 \
--bamlist=bamlist.KAW1063_new.txt \
--posfile=/proj/uppstore2017145/V2/users/chen/GWAS/PA_chr01_1/pos.tabulated.txt \
--outputdir=/proj/uppstore2017145/V2/users/chen/imputation/ \
--K=40 --nGen=100 --nCores=18
