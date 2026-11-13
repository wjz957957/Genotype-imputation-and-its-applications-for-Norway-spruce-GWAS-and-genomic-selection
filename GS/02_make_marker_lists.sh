#!/bin/bash
#SBATCH -J make_marker_lists
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 02:00:00
#SBATCH -o logs/02_make_marker_lists_%j.out
#SBATCH -e logs/02_make_marker_lists_%j.err
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL


set -euo pipefail


############################################
# modules
############################################

module load bioinfo-tools
ml bcftools/1.20


############################################
# paths
############################################

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

VCF=${WORK}/01_reference/Imputed_array_combine.ID.vcf.gz

OUT=${WORK}/01_reference


############################################
# check
############################################

mkdir -p ${OUT}

echo "Input VCF:"
echo ${VCF}


############################################
# chip markers
#
# chip definition:
# INFO field does NOT contain IMP
############################################


bcftools view -H ${VCF} \
| awk '
{
    if(index($8,"IMP")==0)
        print $1":"$2":"$4":"$5
}' \
> ${OUT}/chip.marker


############################################
# imputed markers
#
# imputed definition:
# INFO contains IMP
############################################


bcftools view -H ${VCF} \
| awk '
{
    if(index($8,"IMP")>0)
        print $1":"$2":"$4":"$5
}' \
> ${OUT}/imputed.marker


############################################
# summary
############################################


echo "=============================="
echo "Marker summary"
echo "=============================="

echo "chip markers:"
wc -l ${OUT}/chip.marker


echo "imputed markers:"
wc -l ${OUT}/imputed.marker



############################################
# duplicate check
############################################


echo "Checking duplicate IDs..."

echo "chip duplicates:"
sort ${OUT}/chip.marker | uniq -d | wc -l


echo "imputed duplicates:"
sort ${OUT}/imputed.marker | uniq -d | wc -l



echo "DONE"
