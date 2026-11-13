#!/bin/bash
#SBATCH -J chip_clumping
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/06_chip_clumping_%A_%a.out
#SBATCH -e logs/06_chip_clumping_%A_%a.err
#SBATCH --array=0-59


set -euo pipefail


####################################
# modules
####################################

module load bioinfo-tools
ml plink/1.90b4.9


####################################
# paths
####################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer


GWAS_DIR=${WORK}/05_chip_GWAS

CHIP_DIR=${WORK}/03_chip_selection

OUTDIR=${WORK}/06_chip_clumping



####################################
# task mapping
####################################

TASK=${SLURM_ARRAY_TASK_ID}


# 10 outer × 3 inner × 2 traits

OUTER=$(( TASK / 6 + 1 ))

REMAIN=$(( TASK % 6 ))

INNER=$(( REMAIN / 2 + 1 ))

TRAIT_ID=$(( REMAIN % 2 ))


OUTER=$(printf "%02d" ${OUTER})
INNER=$(printf "%02d" ${INNER})


if [ ${TRAIT_ID} -eq 0 ]; then

    TRAIT="BB"

else

    TRAIT="H6"

fi



echo "=============================="
echo "Task ${TASK}"
echo "Outer ${OUTER}"
echo "Inner ${INNER}"
echo "Trait ${TRAIT}"
echo "=============================="



####################################
# input
####################################


GWAS=${GWAS_DIR}/outer_${OUTER}/inner_${INNER}/chip_${TRAIT}_gwas.loco.mlma


CHIP_PREFIX=${CHIP_DIR}/outer_${OUTER}/inner_${INNER}/chip_train



####################################
# output
####################################


OUT=${OUTDIR}/outer_${OUTER}/inner_${INNER}/${TRAIT}


mkdir -p ${OUT}/clump

mkdir -p ${OUT}/marker_sets



####################################
# Step 1
# Prepare PLINK clump input
####################################


echo "Preparing GWAS p file"


awk '

BEGIN{
    print "SNP\tP"
}

NR>1 && $2!="" && $9!="" {

    print $2"\t"$9

}

' \
${GWAS} \
> ${OUT}/gwas_p.txt



echo "GWAS p file:"
head ${OUT}/gwas_p.txt

echo "Number of GWAS markers:"
wc -l ${OUT}/gwas_p.txt



####################################
# Step 2
# LD clumping
####################################


echo "Running PLINK clumping"


plink \
--bfile ${CHIP_PREFIX} \
--clump ${OUT}/gwas_p.txt \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.2 \
--clump-kb 50 \
--out ${OUT}/clump/chip



####################################
# Step 3
# Extract clump index SNPs
####################################


echo "Generating ranked SNP list"


awk '

NR>1 && $3!="" && $3!="NONE" && $5!="" {

    print $3"\t"$5

}

' \
${OUT}/clump/chip.clumped \
| sort -k2,2g \
> ${OUT}/ranked_chip_snps.txt



echo "Ranked SNP:"
wc -l ${OUT}/ranked_chip_snps.txt


echo "First ranked SNP:"
head ${OUT}/ranked_chip_snps.txt



####################################
# Step 4
# Generate SNP subsets
####################################


for N in \
50 \
100 \
200 \
300 \
400 \
500 \
1000 \
2000 \
3000 \
4000 \
5000 \
8000 \
10000 \
20000 \
30000 \
40000

do


head -${N} \
${OUT}/ranked_chip_snps.txt \
| awk '$1!="" {print $1}' \
> ${OUT}/marker_sets/chip_${N}.snps


echo \
"chip_${N}: $(wc -l < ${OUT}/marker_sets/chip_${N}.snps)"


done



####################################
# Step 5
# All LD-clumped SNPs
####################################


awk '$1!="" {print $1}' \
${OUT}/ranked_chip_snps.txt \
> ${OUT}/marker_sets/chip_LD_ALL.snps



####################################
# Step 6
# All chip SNPs without GWAS selection
####################################


awk '{print $2}' \
${CHIP_PREFIX}.bim \
> ${OUT}/marker_sets/chip_ALL.snps



####################################
# Final check
####################################


echo "=============================="
echo "FINAL CHECK"
echo "=============================="


echo "GWAS input:"
wc -l ${GWAS}


echo "gwas_p:"
wc -l ${OUT}/gwas_p.txt


echo "ranked:"
wc -l ${OUT}/ranked_chip_snps.txt


echo "marker sets:"

for f in ${OUT}/marker_sets/*.snps

do

echo "$(basename ${f}) $(wc -l < ${f})"

done



echo "DONE ${OUTER} ${INNER} ${TRAIT}"
