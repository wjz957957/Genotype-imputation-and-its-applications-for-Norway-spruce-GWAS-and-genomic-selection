#!/bin/bash
#SBATCH -J nested_cv
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 02:00:00
#SBATCH -o logs/03_make_nested_cv_%j.out
#SBATCH -e logs/03_make_nested_cv_%j.err
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL


set -euo pipefail


module load bioinfo-tools
ml bcftools/1.20


WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer


VCF=${WORK}/01_reference/Imputed_array_combine.ID.vcf.gz

OUT=${WORK}/02_CV


mkdir -p ${OUT}
mkdir -p ${WORK}/logs


####################################
# extract samples
####################################

bcftools query -l ${VCF} \
> ${OUT}/all.samples


echo "Total samples:"
wc -l ${OUT}/all.samples


####################################
# create folds in R
####################################


module load PDC/26.03
module load R/4.5.2-cpeGNU-26.03


Rscript - <<EOF


set.seed(2025)


samples <- readLines(
"${OUT}/all.samples"
)


n <- length(samples)


cat(
"Samples:",
n,
"\\n"
)


####################################
# outer 10 fold
####################################


outer <- sample(
rep(1:10,length.out=n)
)


for(i in 1:10){

dir.create(
paste0("${OUT}/outer_",sprintf("%02d",i)),
recursive=TRUE,
showWarnings=FALSE
)


test <- samples[outer==i]

train <- samples[outer!=i]


writeLines(
test,
paste0(
"${OUT}/outer_",
sprintf("%02d",i),
"/outer_test.samples"
)
)


writeLines(
train,
paste0(
"${OUT}/outer_",
sprintf("%02d",i),
"/outer_train.samples"
)
)



####################################
# inner 3 fold
####################################


inner <- sample(
rep(1:3,length.out=length(train))
)


for(j in 1:3){


dir.create(
paste0(
"${OUT}/outer_",
sprintf("%02d",i),
"/inner_",
sprintf("%02d",j)
),
recursive=TRUE,
showWarnings=FALSE
)


itrain <- train[inner!=j]

ival <- train[inner==j]


writeLines(
itrain,
paste0(
"${OUT}/outer_",
sprintf("%02d",i),
"/inner_",
sprintf("%02d",j),
"/inner_train.samples"
)
)


writeLines(
ival,
paste0(
"${OUT}/outer_",
sprintf("%02d",i),
"/inner_",
sprintf("%02d",j),
"/inner_valid.samples"
)
)

}

}


EOF


echo "DONE"
