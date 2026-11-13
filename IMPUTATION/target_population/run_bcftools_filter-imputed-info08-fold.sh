#!/bin/bash
# ============================================================
# 脚本名称: submit_filter_jobs.sh
# 功能: 为每个fold提交独立的SLURM作业进行过滤
# ============================================================

SAMPLE_SIZE=200
REPEAT=1
BASE_DIR="/cfs/klemming/projects/supr/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/target-population-10random/IMPUTE-population/population-${SAMPLE_SIZE}/repeat-${REPEAT}"
INFO_THRESHOLD=0.8

# 创建日志目录
LOG_DIR="${BASE_DIR}/filter_logs"
mkdir -p ${LOG_DIR}

for fold in {1..10}; do
    FOLD_DIR="${BASE_DIR}/fold${fold}"
    JOB_NAME="filter_info_${SAMPLE_SIZE}_${REPEAT}_fold${fold}"
    JOB_SCRIPT="${FOLD_DIR}/filter_job_${JOB_NAME}.slurm"
    
    cat > ${JOB_SCRIPT} << EOF
#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -t 00-4:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH -J ${JOB_NAME}
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL


echo "Job started at: \$(date)"
echo "Processing fold ${fold}"
echo "Sample size: ${SAMPLE_SIZE}"
echo "Repeat: ${REPEAT}"
echo "INFO threshold: > ${INFO_THRESHOLD}"
echo "=========================================="

# 加载模块
module load bioinfo-tools bcftools/1.20

# 进入目录
cd ${FOLD_DIR}

INPUT="impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}.vcf.gz"
OUTPUT="impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}-info08.vcf.gz"
STATS_FILE="filter_stats_fold${fold}.txt"

if [ ! -f "\${INPUT}" ]; then
    echo "Error: Input file not found"
    echo "Status: FAILED - Input file not found" > \${STATS_FILE}
    exit 1
fi

# 统计原始位点数
ORIGINAL=\$(bcftools view -H \${INPUT} | wc -l)
echo "Original SNPs: \${ORIGINAL}"

# 过滤：保留INFO > 阈值的位点
echo "Filtering SNPs with INFO > ${INFO_THRESHOLD}..."
bcftools filter -i "INFO > ${INFO_THRESHOLD}" \${INPUT} -Oz -o \${OUTPUT}

if [ \$? -ne 0 ]; then
    echo "Error: bcftools filter failed"
    echo "Status: FAILED - bcftools filter failed" > \${STATS_FILE}
    exit 1
fi

# 建立索引
bcftools index \${OUTPUT}

# 统计过滤后的位点数
FILTERED=\$(bcftools view -H \${OUTPUT} | wc -l)
REMOVED=\$((ORIGINAL - FILTERED))

if [ \${ORIGINAL} -gt 0 ]; then
    RETENTION=\$(echo "scale=2; \${FILTERED} * 100 / \${ORIGINAL}" | bc)
else
    RETENTION="0"
fi

echo "=========================================="
echo "Filtering Results:"
echo "  Original SNPs: \${ORIGINAL}"
echo "  Filtered SNPs (INFO > ${INFO_THRESHOLD}): \${FILTERED}"
echo "  Removed SNPs: \${REMOVED}"
echo "  Retention rate: \${RETENTION}%"
echo "=========================================="

# 保存统计信息到当前fold文件夹
cat > \${STATS_FILE} << STATS_EOF
# Filter Statistics for Fold ${fold}
# ==========================================
Filter date: \$(date)
Input file: \${INPUT}
Output file: \${OUTPUT}
Sample size: ${SAMPLE_SIZE}
Repeat: ${REPEAT}
INFO threshold: > ${INFO_THRESHOLD}
Original SNPs: \${ORIGINAL}
Filtered SNPs: \${FILTERED}
Removed SNPs: \${REMOVED}
Retention rate: \${RETENTION}%
Status: SUCCESS
Output file size: \$(ls -lh \${OUTPUT} | awk '{print \$5}')
# ==========================================
STATS_EOF

echo "Statistics saved to: \${STATS_FILE}"
echo "Job completed successfully at: \$(date)"
echo "=========================================="

EOF
    
    chmod +x ${JOB_SCRIPT}
    
    echo "Submitting job for fold ${fold}..."
    sbatch ${JOB_SCRIPT}
    
    # 避免同时提交过多作业
    sleep 0.5
done

echo ""
echo "=========================================="
echo "All jobs submitted!"
echo "=========================================="
echo "After jobs complete, statistics will be saved in each fold folder:"
for fold in {1..10}; do
    echo "  fold${fold}/filter_stats_fold${fold}.txt"
done
echo "=========================================="
