#!/bin/bash

################ JOB AND SCHEDULING SPECS
#SBATCH -t 05-00:00:00     # jobtime
#SBATCH -A naiss2023-22-925   # the account to be charged
#SBATCH -J wjz-all     # a job name
#SBATCH -p node
#SBATCH -C mem256GB
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END     # or END, or others
#SBATCH -o shapeCom_%j.out  # standard output with the job number in the name
#SBATCH -e shapeCom_%j.err  # standard error with the job number in the name

cat $0  # print the script. useful for debugging.

echo -e "\n\n\thello from:" $HOSTNAME
echo -e "\n\n\tJob started on:"
date
################

################ RUNNING COMMANDS

module load bioinfo-tools SHAPEIT/v5.1.1

echo $1		# so we know which chr is being processed

SHAPEIT5_phase_common -T 17 \
	-I ${1}_Miss_missing0.2_maf0.05_hwe2.5e-10_info0.8_plink2dosage.vcf.gz \
	-O ${1}_Miss_missing0.2_maf0.05_hwe2.5e-10_info0.8_plink2dosage_phased_common_map.bcf \
	--region $2 \
        -M LD.GenMap.chr11.txt

################

echo -e "\n\n\tJob finished on:"
date

