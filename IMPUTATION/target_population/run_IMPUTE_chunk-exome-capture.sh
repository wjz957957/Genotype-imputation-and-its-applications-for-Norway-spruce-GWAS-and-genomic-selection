#!/bin/bash

# 手动指定的区间列表，染色体固定为chr11
region1s=(
    "168954-100000000"
    "100000001-200000000"
	"200000001-300000000"
	"300000001-400000000"
	"400000001-500000000"
	"500000001-600000000"
	"600000001-700000000"
	"700000001-800000000"
	"800000001-900000000"
	"900000001-953660116"
)

region2s=(
    "1-100050000"
    "99950001-200050000"
	"199950001-300050000"
	"299950001-400050000"
	"399950001-500050000"
	"499950001-600050000"
	"599950001-700050000"
	"699950001-800050000"
	"799950001-900050000"
	"899950001-954579047"
)

# 确保两个区间列表长度相同
if [ ${#region1s[@]} -ne ${#region2s[@]} ]; then
    echo "Error: The regions arrays must have the same length."
    exit 1
fi

# 为每对区间创建SLURM作业
for i in ${!region1s[@]}; do
    region1=${region1s[$i]}
    region2=${region2s[$i]}

    # 为作业生成一个基于区间的文件名
    JOB_SCRIPT="job_chr11-1297_${region1}.slurm"
    JOB_NAME="job_chr11-1297_${region1}"

    # 创建SLURM作业脚本
    cat <<EOT > $JOB_SCRIPT
#!/bin/bash

#SBATCH -A naiss2023-5-58  
#SBATCH -t 05-00:00:00   
#SBATCH -p core -n 5
#SBATCH --job-name=$JOB_NAME
#SBATCH --mail-user 253721079@qq.com
#SBATCH --mail-type=ALL

# 加载必要的模块
module load bioinfo-tools bcftools/1.20

##进行IMPUTE填充
/crex/proj/snic2020-6-176/nobackup/private/jingz/impute5_v1.2.1_x86_64_static \
--h Chr11_map_phased.vcf.gz \
--g exome-11-1297-reserve-phased.vcf.gz \
--m LD.GenMap.chr11.txt \
--o impute5_output-1297-${region1}.vcf.gz \
--r 11:${region1} \
--buffer-region 11:${region2} \
--threads 2

bcftools index impute5_output-1297-${region1}.vcf.gz


EOT

    # 提交SLURM作业
    sbatch $JOB_SCRIPT


done
