#!/bin/bash
#SBATCH -J chip_LD_pruning
#SBATCH -A naiss2025-5-732
#SBATCH -p shared

#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL

#SBATCH -o logs/05_chip_LD_pruning_%A_%a.out
#SBATCH -e logs/05_chip_LD_pruning_%A_%a.err

#SBATCH --array=0-29


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


INPUT=${WORK}/03_chip_selection

OUTPUT=${WORK}/04_chip_pruning



####################################
# map task
####################################

TASK=${SLURM_ARRAY_TASK_ID}


OUTER=$(( TASK / 3 + 1 ))
INNER=$(( TASK % 3 + 1 ))


OUTER=$(printf "%02d" ${OUTER})
INNER=$(printf "%02d" ${INNER})


echo "=============================="
echo "Task : ${TASK}"
echo "Outer: ${OUTER}"
echo "Inner: ${INNER}"
echo "=============================="



####################################
# input
####################################

BFILE=${INPUT}/outer_${OUTER}/inner_${INNER}/chip_train



####################################
# output
####################################

OUT=${OUTPUT}/outer_${OUTER}/inner_${INNER}

mkdir -p ${OUT}


PREFIX=${OUT}/chip_pruned



####################################
# LD pruning
####################################


echo "Running LD pruning..."


plink \
 --bfile ${BFILE} \
 --indep-pairwise 50 5 0.2 \
 --allow-extra-chr \
 --out ${PREFIX}



####################################
# extract LD-pruned markers
####################################


plink \
 --bfile ${BFILE} \
 --extract ${PREFIX}.prune.in \
 --make-bed \
 --allow-extra-chr \
 --out ${PREFIX}



####################################
# save snplist
####################################


cp ${PREFIX}.prune.in ${PREFIX}.snplist



####################################
# check
####################################


echo "=============================="
echo "CHECK"
echo "=============================="

echo "Samples:"
wc -l ${PREFIX}.fam


echo "LD-pruned markers:"
wc -l ${PREFIX}.bim


echo "DONE outer_${OUTER} inner_${INNER}"
