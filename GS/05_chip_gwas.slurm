#!/bin/bash
#SBATCH -J chip_gwas
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/05_chip_gwas_%A_%a.out
#SBATCH -e logs/05_chip_gwas_%A_%a.err
#SBATCH --array=0-29


set -euo pipefail


####################################
# modules
####################################

module load bioinfo-tools
ml gcta/1.94.1


####################################
# paths
####################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer


CHIP_DIR=${WORK}/03_chip_selection

CV_DIR=${WORK}/02_CV

PHENO_DIR=${WORK}/00_input


OUTDIR=${WORK}/05_chip_GWAS



####################################
# task mapping
####################################

TASK=${SLURM_ARRAY_TASK_ID}


OUTER=$(( TASK / 3 + 1 ))
INNER=$(( TASK % 3 + 1 ))


OUTER=$(printf "%02d" ${OUTER})
INNER=$(printf "%02d" ${INNER})


echo "=============================="
echo "Task ${TASK}"
echo "Outer ${OUTER}"
echo "Inner ${INNER}"
echo "=============================="


####################################
# input
####################################

CHIP_PREFIX=${CHIP_DIR}/outer_${OUTER}/inner_${INNER}/chip_train


KEEP=${CV_DIR}/outer_${OUTER}/inner_${INNER}/inner_train.samples



####################################
# output
####################################

OUT=${OUTDIR}/outer_${OUTER}/inner_${INNER}

mkdir -p ${OUT}/grm
mkdir -p ${OUT}/pheno



####################################
# phenotype prepare
####################################


awk '
NR==FNR{
    a[$1]=$3
    next
}
{
    if($1 in a)
        print $1,$1,a[$1]
}
' \
${PHENO_DIR}/Trait-Adj_Budburst.txt \
${KEEP} \
> ${OUT}/pheno/BB.pheno



awk '
NR==FNR{
    a[$1]=$3
    next
}
{
    if($1 in a)
        print $1,$1,a[$1]
}
' \
${PHENO_DIR}/Trait-Adj_H6.txt \
${KEEP} \
> ${OUT}/pheno/H6.pheno



echo "Phenotype check"

wc -l ${OUT}/pheno/*.pheno



####################################
# GRM
####################################

echo "Building GRM"


gcta64 \
--bfile ${CHIP_PREFIX} \
--make-grm \
--out ${OUT}/grm/chip_grm



####################################
# BB GWAS
####################################


echo "Running BB MLMA-LOCO"


gcta64 \
--mlma-loco \
--bfile ${CHIP_PREFIX} \
--grm ${OUT}/grm/chip_grm \
--pheno ${OUT}/pheno/BB.pheno \
--out ${OUT}/chip_BB_gwas



####################################
# H6 GWAS
####################################


echo "Running H6 MLMA-LOCO"


gcta64 \
--mlma-loco \
--bfile ${CHIP_PREFIX} \
--grm ${OUT}/grm/chip_grm \
--pheno ${OUT}/pheno/H6.pheno \
--out ${OUT}/chip_H6_gwas



echo "DONE outer_${OUTER} inner_${INNER}"
