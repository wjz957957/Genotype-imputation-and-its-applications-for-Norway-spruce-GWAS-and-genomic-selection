#!/bin/bash
#SBATCH -J prepare_chip_genotype
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/07_prepare_chip_genotype_%A_%a.out
#SBATCH -e logs/07_prepare_chip_genotype_%A_%a.err
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


CHIP_DIR=${WORK}/03_chip_selection


MARKER_DIR=${WORK}/06_chip_clumping


CV_DIR=${WORK}/02_CV


OUTDIR=${WORK}/07_chip_RRBLUP



####################################
# task mapping
####################################

TASK=${SLURM_ARRAY_TASK_ID}


OUTER=$((TASK / 3 + 1))
INNER=$((TASK % 3 + 1))


OUTER=$(printf "%02d" ${OUTER})
INNER=$(printf "%02d" ${INNER})


echo "======================"
echo "Outer ${OUTER}"
echo "Inner ${INNER}"
echo "======================"



####################################
# input
####################################

TRAIN_PREFIX=${CHIP_DIR}/outer_${OUTER}/inner_${INNER}/chip_train


VALID_SAMPLE=${CV_DIR}/outer_${OUTER}/inner_${INNER}/inner_valid.samples



####################################
# output
####################################

OUT=${OUTDIR}/outer_${OUTER}/inner_${INNER}


mkdir -p ${OUT}/genotype



####################################
# make validation keep
####################################


KEEP=${OUT}/valid.keep


awk '
{
print $1"\t"$1
}
' \
${VALID_SAMPLE} \
> ${KEEP}



echo "Validation keep:"
wc -l ${KEEP}

head ${KEEP}



if [ ! -s ${KEEP} ]
then
    echo "ERROR: empty keep"
    exit 1
fi



####################################
# marker sets
####################################


MARKERS=(
50
100
200
300
400
500
1000
2000
3000
4000
5000
8000
10000
20000
30000
40000
LD_ALL
ALL
)



####################################
# prepare genotype
####################################


for TRAIT in BB H6
do


echo "======================"
echo "Trait ${TRAIT}"
echo "======================"


for M in "${MARKERS[@]}"
do


echo "Preparing ${TRAIT} ${M}"


SNPS=${MARKER_DIR}/outer_${OUTER}/inner_${INNER}/${TRAIT}/marker_sets/chip_${M}.snps


PREFIX=${OUT}/genotype/chip_${TRAIT}_${M}



####################################
# training genotype
####################################


plink \
--bfile ${TRAIN_PREFIX} \
--extract ${SNPS} \
--recode A \
--out ${PREFIX}_train



####################################
# validation genotype
####################################


plink \
--vcf ${VCF} \
--extract ${SNPS} \
--keep ${KEEP} \
--double-id \
--allow-extra-chr \
--recode A \
--out ${PREFIX}_valid



done


done



echo "DONE outer_${OUTER} inner_${INNER}"
