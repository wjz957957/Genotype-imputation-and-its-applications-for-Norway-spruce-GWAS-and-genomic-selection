#!/bin/bash
#SBATCH -J run_chip_rrblup
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-29%5
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/07_run_chip_rrblup_%A_%a.out
#SBATCH -e logs/07_run_chip_rrblup_%A_%a.err

set -euo pipefail

####################################
# Modules
####################################

module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

####################################
# Paths
####################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

SCRIPT=${WORK}/scripts/07_run_chip_rrblup.R

####################################
# Fold mapping
####################################

TASK=${SLURM_ARRAY_TASK_ID}

OUTER=$((TASK / 3 + 1))
INNER=$((TASK % 3 + 1))

OUTER=$(printf "%02d" "${OUTER}")
INNER=$(printf "%02d" "${INNER}")

echo "================================"
echo "Task  : ${TASK}"
echo "Outer : ${OUTER}"
echo "Inner : ${INNER}"
echo "================================"

####################################
# Input checks
####################################

test -s "${SCRIPT}"

GENO_DIR=${WORK}/07_chip_RRBLUP/outer_${OUTER}/inner_${INNER}/genotype

for TRAIT in BB H6
do
    for SET in 50 100 200 300 400 500 1000 2000 3000 \
               4000 5000 8000 10000 20000 30000 40000 \
               LD_ALL ALL
    do
        test -s "${GENO_DIR}/chip_${TRAIT}_${SET}_train.raw"
        test -s "${GENO_DIR}/chip_${TRAIT}_${SET}_valid.raw"
    done
done

####################################
# RRBLUP
####################################

cd "${WORK}"

Rscript "${SCRIPT}" \
    "${WORK}" \
    "${OUTER}" \
    "${INNER}"

echo "DONE Step7 outer_${OUTER} inner_${INNER}"
