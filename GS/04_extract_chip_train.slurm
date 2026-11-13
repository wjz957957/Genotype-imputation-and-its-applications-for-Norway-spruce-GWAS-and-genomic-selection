#!/bin/bash
#SBATCH -J extract_chip_train
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL

#SBATCH -o logs/04_extract_chip_train_%A_%a.out
#SBATCH -e logs/04_extract_chip_train_%A_%a.err

#SBATCH --array=0-29


set -euo pipefail


####################################
# modules
####################################

module load bioinfo-tools
ml bcftools/1.20
ml plink/1.90b4.9


####################################
# paths
####################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer


VCF=${WORK}/01_reference/Imputed_array_combine.ID.vcf.gz

CHIP_MARKER=${WORK}/01_reference/chip.marker

CV_DIR=${WORK}/02_CV

OUTDIR=${WORK}/03_chip_selection



####################################
# map array ID
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
# output directory
####################################

OUT=${OUTDIR}/outer_${OUTER}/inner_${INNER}

mkdir -p ${OUT}


PREFIX=${OUT}/chip_train



####################################
# prepare PLINK keep file
####################################

KEEP_RAW=${CV_DIR}/outer_${OUTER}/inner_${INNER}/inner_train.samples

KEEP=${OUT}/inner_train.keep


if [ ! -f ${KEEP_RAW} ]; then
    echo "ERROR: missing ${KEEP_RAW}"
    exit 1
fi


awk '{print $1"\t"$1}' ${KEEP_RAW} > ${KEEP}



echo "KEEP samples:"
wc -l ${KEEP}


echo "First samples:"
head ${KEEP}



####################################
# extract chip SNPs
####################################

plink \
 --vcf ${VCF} \
 --extract ${CHIP_MARKER} \
 --keep ${KEEP} \
 --make-bed \
 --double-id \
 --allow-extra-chr \
 --out ${PREFIX}



####################################
# quality check
####################################


echo "=============================="
echo "CHECK"
echo "=============================="


echo "Samples:"
wc -l ${PREFIX}.fam


echo "Markers:"
wc -l ${PREFIX}.bim



echo "DONE outer_${OUTER} inner_${INNER}"
