#!/bin/bash
# ============================================================
# 批量提交 concat 和提取作业
# ============================================================

# 参数设置
SAMPLE_SIZE=200
REPEAT=1

# 固定路径
BASE_DIR="/cfs/klemming/projects/supr/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/target-population-10random/"
SITE_FOLDS_DIR="${BASE_DIR}/site_fold-no-phasing"
OUTPUT_BASE="${BASE_DIR}/IMPUTE-population/population-${SAMPLE_SIZE}/repeat-${REPEAT}"

# ============================================================
# 提交所有fold的作业
# ============================================================
for fold in {1..10}; do
    echo "=========================================="
    echo "Processing fold ${fold}"
    echo "=========================================="
    
    # 定义作业名称
    JOB_NAME="concat_${SAMPLE_SIZE}_${REPEAT}_fold${fold}"
    
    # 定义输出目录
    FOLD_DIR="${OUTPUT_BASE}/fold${fold}"
    LIST_FILE="${OUTPUT_BASE}/list-${SAMPLE_SIZE}-${REPEAT}-fold${fold}.txt"
    
    # 创建SLURM作业脚本
    JOB_SCRIPT="${FOLD_DIR}/job_${JOB_NAME}.slurm"
    
    # 确保目录存在
    mkdir -p ${FOLD_DIR}
    
    cat <<EOT > ${JOB_SCRIPT}
#!/bin/bash
#SBATCH -A naiss2025-5-732
#SBATCH -t 00-3:00:00
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -J ${JOB_NAME}
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL


echo "=========================================="
echo "Job started at: \$(date)"
echo "Sample size: ${SAMPLE_SIZE}"
echo "Repeat: ${REPEAT}"
echo "Fold: ${fold}"
echo "Working directory: ${FOLD_DIR}"
echo "=========================================="

# 加载模块
module load bioinfo-tools bcftools/1.20

# 进入工作目录
cd ${FOLD_DIR}

# 检查list文件是否存在
if [ ! -f "${LIST_FILE}" ]; then
    echo "Error: List file not found: ${LIST_FILE}"
    exit 1
fi

# ============================================================
# 步骤1: 合并所有region的VCF文件
# ============================================================
echo "Step 1: Concatenating VCF files..."
bcftools concat -n -f ${LIST_FILE} -Oz -o impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}.vcf.gz

if [ \$? -ne 0 ]; then
    echo "Error: bcftools concat failed"
    exit 1
fi

# ============================================================
# 步骤2: 建立索引
# ============================================================
echo "Step 2: Indexing concatenated VCF..."
bcftools index impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}.vcf.gz

# ============================================================
# 步骤3: 提取评估位点（基于SNP ID）
# ============================================================
echo "Step 3: Extracting evaluation SNPs..."
EVAL_FILE="${SITE_FOLDS_DIR}/fold_${fold}_eval_pos.txt"

if [ ! -f "\${EVAL_FILE}" ]; then
    echo "Error: Evaluation SNP list not found: \${EVAL_FILE}"
    exit 1
fi

# 使用 -R 基于SNP ID进行筛选
bcftools view \
    -R "\${EVAL_FILE}" \
    impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}.vcf.gz \
    -Oz -o impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}_real.vcf.gz

if [ \$? -ne 0 ]; then
    echo "Error: bcftools view failed"
    exit 1
fi
# ============================================================
# 步骤4: 建立索引
# ============================================================
echo "Step 4: Indexing extracted VCF..."
bcftools index impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}_real.vcf.gz

# ============================================================
# 步骤5: 统计信息
# ============================================================
echo "Step 5: Statistics..."
echo "Original VCF (concatenated):"
bcftools view -H impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}.vcf.gz | wc -l
echo "Extracted VCF (evaluation SNPs only):"
bcftools view -H impute5_output_${SAMPLE_SIZE}_${REPEAT}_fold${fold}_real.vcf.gz | wc -l

echo "=========================================="
echo "Job completed at: \$(date)"
echo "=========================================="

EOT

    # 提交作业
    echo "Submitting job: ${JOB_NAME}"
    sbatch ${JOB_SCRIPT}
    
    # 等待1秒，避免同时提交过多作业
    sleep 1
    
done

echo ""
echo "=========================================="
echo "All jobs submitted for sample_size=${SAMPLE_SIZE}, repeat=${REPEAT}"
echo "Total jobs: 10"
echo "=========================================="

