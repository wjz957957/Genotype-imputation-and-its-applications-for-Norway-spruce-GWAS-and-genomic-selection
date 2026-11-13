#!/bin/bash

#SBATCH -J 10b2_outer_imputed_gwas
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-239%6
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/10b2_outer_imputed_gwas_%A_%a.out
#SBATCH -e logs/10b2_outer_imputed_gwas_%A_%a.err


set -Eeuo pipefail

trap 'rc=$?; echo "ERROR line ${LINENO}, exit=${rc}: ${BASH_COMMAND}" >&2' ERR


########################################
# Modules
########################################

module load bioinfo-tools
ml gcta/1.94.1


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}



########################################
# Paths
########################################


WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer


INPUT_ROOT=${WORK}/10_outer_test/02_imputed_inputs

OUTROOT=${WORK}/10_outer_test/03_imputed_GWAS



########################################
# Task mapping
#
# 10 outer
# 12 chromosomes
# 2 traits
#
########################################


TASK=${SLURM_ARRAY_TASK_ID}


TRAIT_ID=$((TASK % 2))

CHR=$(( (TASK / 2) % 12 + 1 ))

OUTER=$(( TASK / 24 + 1 ))


OUTER=$(printf "%02d" ${OUTER})

CHR_PAD=$(printf "%02d" ${CHR})


if [[ ${TRAIT_ID} -eq 0 ]]
then
    TRAIT="BB"
else
    TRAIT="H6"
fi



echo "========================================"
echo "Step10b-2: Outer imputed GWAS"
echo "Outer : ${OUTER}"
echo "Chr   : ${CHR}"
echo "Trait : ${TRAIT}"
echo "========================================"



########################################
# Input paths
########################################


INPUT_DIR=${INPUT_ROOT}/outer_${OUTER}


BFILE=${INPUT_DIR}/plink/imputed_train_qc


GRM_DIR=${INPUT_DIR}/grm


GRM_ALL=${GRM_DIR}/chip_all

GRM_CHR=${GRM_DIR}/chip_chr${CHR_PAD}



OUT=${OUTROOT}/outer_${OUTER}/chr_${CHR_PAD}/${TRAIT}


mkdir -p ${OUT}



########################################
# Input checks
########################################


if [[ ! -s ${INPUT_DIR}/complete.flag ]]
then
    echo "ERROR: Step10b-1 incomplete"
    exit 1
fi



for EXT in bed bim fam
do

    if [[ ! -s ${BFILE}.${EXT} ]]
    then
        echo "ERROR missing genotype ${BFILE}.${EXT}"
        exit 1
    fi

done



N_TRAIN=$(wc -l < ${BFILE}.fam)


echo "Training samples: ${N_TRAIN}"



if [[ ${N_TRAIN} -lt 700 ||
      ${N_TRAIN} -gt 900 ]]
then

    echo "ERROR unexpected sample number"

    exit 1

fi



########################################
# GRM checks
########################################


for PREFIX in \
    ${GRM_ALL} \
    ${GRM_CHR}

do

    for EXT in grm.bin grm.N.bin grm.id
    do

        if [[ ! -s ${PREFIX}.${EXT} ]]
        then

            echo "ERROR missing GRM ${PREFIX}.${EXT}"

            exit 1

        fi

    done

done



if [[ \
$(wc -l < ${GRM_ALL}.grm.id) \
-ne \
${N_TRAIN}
]]
then

    echo "ERROR GRM sample mismatch"

    exit 1

fi



########################################
# Chromosome SNP count
########################################


CHR_COUNT_FILE=${INPUT_DIR}/imp_chr_counts.tsv


if [[ ! -s ${CHR_COUNT_FILE} ]]
then

    echo "ERROR missing chromosome SNP count"

    exit 1

fi



N_CHR_SNP=$(awk \
-v chr=${CHR} '

$1==chr {
    print $2
    found=1
}

END{
    if(!found)
        exit 1
}

' ${CHR_COUNT_FILE})



if [[ -z ${N_CHR_SNP} ||
      ${N_CHR_SNP} -lt 1000 ]]
then

    echo "ERROR too few SNPs"

    exit 1

fi



echo "Chromosome SNPs: ${N_CHR_SNP}"



########################################
# Phenotype
########################################
#
# IMPORTANT:
# phenotype files are under pheno/
#
########################################


PHENO=${INPUT_DIR}/pheno/${TRAIT}.pheno



if [[ ! -s ${PHENO} ]]
then

    echo "ERROR missing phenotype:"
    echo ${PHENO}

    exit 1

fi



if [[ $(wc -l < ${PHENO}) -ne ${N_TRAIN} ]]
then

    echo "ERROR phenotype count mismatch"

    exit 1

fi



########################################
# Run GCTA
########################################


RESULT=${OUT}/${TRAIT}_chr${CHR_PAD}.mlma



if [[ -s ${RESULT} &&
      -s ${OUT}/complete.flag ]]
then

    echo "Already completed"

    exit 0

fi



rm -f ${OUT}/complete.flag



echo "========================================"
echo "Running GCTA MLMA"
echo "========================================"



gcta64 \
--mlma \
--autosome-num 12 \
--chr ${CHR} \
--bfile ${BFILE} \
--grm ${GRM_ALL} \
--mlma-subtract-grm ${GRM_CHR} \
--pheno ${PHENO} \
--thread-num 8 \
--out ${OUT}/${TRAIT}_chr${CHR_PAD}



########################################
# Output validation
########################################


if [[ ! -s ${RESULT} ]]
then

    echo "ERROR missing GWAS result"

    exit 1

fi



awk '

NR==1 {

    if($1!="Chr" ||
       $2!="SNP" ||
       $9!="p")
    {
        print "ERROR header"
        exit 1
    }

    next

}


{

tested++

if($9!="nan" &&
   $9!="NA" &&
   $9!="")
{
valid++
}
else
{
invalid++
}

}


END{

printf "Tested\t%d\nValid\t%d\nInvalid\t%d\n",
tested,valid,invalid

}

' ${RESULT} \
> ${OUT}/gwas_qc.tsv



########################################
# Complete
########################################


echo "DONE Step10b-2 outer_${OUTER} chr_${CHR_PAD} ${TRAIT}" \
> ${OUT}/complete.flag



echo "========================================"
echo "DONE"
echo "Outer : ${OUTER}"
echo "Chr   : ${CHR}"
echo "Trait : ${TRAIT}"
echo "========================================"
