#!/bin/bash


# Set the starting piece index to 1 (fourth piece in zero-based indexing)
START_INDEX=1

# 循环生成每个窗口的参数并创建并提交SLURM作业
for (( i=START_INDEX; i<4; i++ )); do
    

    # 输出当前窗口的参数
    echo "Creating SLURM job for impute-population-${i}"

    # 创建一个临时的SLURM作业脚本
    JOB_SCRIPT="job-600_${i}.slurm"
    cat <<EOT > $JOB_SCRIPT
#!/bin/bash

#SBATCH -A naiss2023-5-58
#SBATCH -t 05-00:00:00   
#SBATCH -p core -n 2
#SBATCH -J phase_job_$i
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL


# Load any necessary modules
module load bioinfo-tools bcftools/1.20
module load bioinfo-tools SHAPEIT/v5.1.1


# 调用 bcftools 命令
bcftools view -S impute-population-600-${i}.txt exome-11-1297-reserve.vcf.gz -Oz -o exome-population-600-${i}.vcf.gz

bcftools index exome-population-600-${i}.vcf.gz

bcftools view -S impute-population-600-${i}.txt exome-11-1297-real.vcf.gz -Oz -o exome-population-600-real-${i}.vcf.gz

bcftools index exome-population-600-real-${i}.vcf.gz

SHAPEIT5_phase_common -T 10 \
    -I exome-population-600-${i}.vcf.gz \
    -O exome-population-phased-600-${i}.bcf \
    --region 11 \
    -M LD.GenMap.chr11.txt

bcftools view  exome-population-phased-600-${i}.bcf -Ov -o - | bgzip -c > exome-population-phased-600-${i}.vcf.gz

bcftools index exome-population-phased-600-${i}.vcf.gz


EOT

    # 提交SLURM作业
    sbatch $JOB_SCRIPT

done


