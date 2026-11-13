	#!/bin/bash
# ============================================================
# 脚本名称: run_all_folds.sh
# 功能: 为每个fold生成评估脚本并运行
# 运行位置: /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/target-population-10random/IMPUTE-population/population-200/repeat-1
# ============================================================

# 设置参数
SAMPLE_SIZE=200
REPEAT=1
BASE_DIR="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/target-population-10random/IMPUTE-population/population-${SAMPLE_SIZE}/repeat-${REPEAT}"
SITE_FOLDS_DIR="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/target-population-10random/site_fold-no-phasing"

# 切换到基础目录
cd ${BASE_DIR} || exit 1

echo "=========================================="
echo "Running evaluation for all folds"
echo "Base directory: ${BASE_DIR}"
echo "Sample size: ${SAMPLE_SIZE}"
echo "Repeat: ${REPEAT}"
echo "=========================================="

# 循环处理每个fold
for fold in {1..10}; do
    FOLD_DIR="${BASE_DIR}/fold${fold}"
    
    # 检查fold目录是否存在
    if [ ! -d "${FOLD_DIR}" ]; then
        echo "Warning: Directory ${FOLD_DIR} not found, skipping..."
        continue
    fi
    
    echo ""
    echo "=========================================="
    echo "Processing fold ${fold}"
    echo "Directory: ${FOLD_DIR}"
    echo "=========================================="
    
    # ============================================================
    # 1. 生成 Concordance 评估脚本（添加summary功能）
    # ============================================================
    cat > ${FOLD_DIR}/IMPUTE-reference-${SAMPLE_SIZE}-concordance-nophasing-for.py << 'CONCORDANCE_EOF'
#!/usr/bin/env python3
# ============================================================
# 脚本: IMPUTE-reference-200-concordance-for.py
# 功能: 计算填充准确性 (Concordance) - 使用矩阵累加算法
# 用法: python IMPUTE-reference-200-concordance-for.py
# 运行位置: 在 fold1, fold2, ... 目录下运行
# ============================================================

import gzip
import numpy as np
import os
import sys

# ============================================================
# 参数设置
# ============================================================
SAMPLE_SIZE = 200
REPEAT = 1

# 获取当前目录名，确定 fold 编号
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith('fold'):
    FOLD = int(current_dir.replace('fold', ''))
else:
    print("Error: This script must be run in a fold directory (fold1, fold2, ...)")
    sys.exit(1)

# 文件路径
# 填充结果 VCF: GT:DS:GP 格式，位于当前目录
IMPUTED_VCF = f"impute5_output_{SAMPLE_SIZE}_{REPEAT}_fold{FOLD}_real.vcf.gz"
# 真实基因型 VCF: GT:DS:GP 格式
TRUE_VCF = f"exome-population-{SAMPLE_SIZE}-{REPEAT}-fold{FOLD}-eval.vcf.gz"
# 评估位点列表
EVAL_SNPS_FILE = f"/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/target-population-10random/site_fold-no-phasing/fold_{FOLD}_eval.txt"

print("=" * 60)
print(f"Concordance Evaluation - Fold {FOLD}")
print(f"Sample size: {SAMPLE_SIZE}")
print(f"Repeat: {REPEAT}")
print(f"Working directory: {os.getcwd()}")
print("=" * 60)

# ============================================================
# 第一步：读取真实基因型 (GT 格式: 0/0, 0/1, 1/0, 1/1)
# ============================================================
print("\n[1/4] Reading true genotypes from evaluation VCF...")

# 定义真实基因型到数字的映射
true_genotype_map = {
    '0/0': 0,
    '0/1': 1,
    '1/0': 1,
    '1/1': 2
}

true_genotypes = []  # 存储每行的真实基因型数字列表

with gzip.open(TRUE_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        if len(columns) > 9:
            # 提取所有样本的 GT 字段
            gt_values = []
            for sample in columns[9:]:
                gt = sample.split(':')[0]  # 提取 GT 字段
                if gt in true_genotype_map:
                    gt_values.append(true_genotype_map[gt])
                else:
                    gt_values.append('NA')  # 缺失值标记为 'NA'
            true_genotypes.append(gt_values)

print(f"  Found {len(true_genotypes)} SNPs in true VCF")

# ============================================================
# 第二步：读取填充基因型 (GP 字段)
# ============================================================
print("\n[2/4] Reading imputed GP values from imputed VCF...")

imputed_gp = []  # 存储每行的 GP 值列表

with gzip.open(IMPUTED_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        if len(columns) > 9:
            # 提取所有样本的 GP 字段
            gp_values = []
            for sample in columns[9:]:
                parts = sample.split(':')
                if len(parts) >= 3:
                    # 标准格式: GT:DS:GP
                    try:
                        gp = [float(x) for x in parts[2].split(',')]
                        gp_values.append(gp)
                    except:
                        gp_values.append([np.nan, np.nan, np.nan])
                elif len(parts) == 2:
                    # 格式: GT:GP
                    try:
                        gp = [float(x) for x in parts[1].split(',')]
                        gp_values.append(gp)
                    except:
                        gp_values.append([np.nan, np.nan, np.nan])
                else:
                    gp_values.append([np.nan, np.nan, np.nan])
            imputed_gp.append(gp_values)

print(f"  Found {len(imputed_gp)} SNPs in imputed VCF")

# ============================================================
# 第三步：确保两个文件的 SNP 数量一致
# ============================================================
print("\n[3/4] Aligning SNPs...")

min_snps = min(len(true_genotypes), len(imputed_gp))
if len(true_genotypes) != len(imputed_gp):
    print(f"  Warning: True SNPs ({len(true_genotypes)}) != Imputed SNPs ({len(imputed_gp)})")
    print(f"  Using first {min_snps} SNPs for analysis")

true_genotypes = true_genotypes[:min_snps]
imputed_gp = imputed_gp[:min_snps]

print(f"  Processing {min_snps} SNPs")

# ============================================================
# 第四步：计算每个 SNP 的 Concordance (使用矩阵累加算法)
# ============================================================
print("\n[4/4] Calculating concordance using matrix accumulation...")

def create_matrix_from_line(true_gt_list, gp_list):
    """
    根据真实基因型和 GP 概率创建混淆矩阵
    
    矩阵结构:
          预测0   预测1   预测2   总计
    真实0   a       b       c     sum0
    真实1   d       e       f     sum1
    真实2   g       h       i     sum2
    总计   sum0'   sum1'   sum2'  total
    """
    # 初始化 4x4 矩阵
    matrix = [[0.0 for _ in range(4)] for _ in range(4)]
    
    # 遍历每个样本
    for true_gt, gp in zip(true_gt_list, gp_list):
        # 跳过缺失的真实基因型
        if true_gt == 'NA' or true_gt is None:
            continue
        
        true_idx = int(true_gt)  # 真实基因型索引 (0, 1, 2)
        
        # 跳过 GP 缺失的样本
        if np.isnan(gp[0]):
            continue
        
        # 将三个概率值累加到矩阵的对应列
        for pred_idx, prob in enumerate(gp):
            if pred_idx < 3:  # 只处理前三个预测类别
                matrix[pred_idx][true_idx] += prob
    
    # 计算总计行和列
    for i in range(3):
        matrix[i][3] = sum(matrix[i][:3])      # 行总计
        matrix[3][i] = sum(matrix[j][i] for j in range(3))  # 列总计
    matrix[3][3] = sum(matrix[i][3] for i in range(3))  # 总样本数
    
    return matrix

def calculate_concordance(matrix):
    """
    从混淆矩阵计算加权一致性率 (Po)
    Po = 对角线元素之和 / 总样本数
    """
    diagonal_sum = sum(matrix[i][i] for i in range(3))
    total = matrix[3][3]
    
    if total > 0:
        return diagonal_sum / total
    else:
        return np.nan

# 计算每个 SNP 的 concordance
concordance_values = []

for snp_idx in range(min_snps):
    true_gt_row = true_genotypes[snp_idx]
    gp_row = imputed_gp[snp_idx]
    
    # 确保长度一致
    min_len = min(len(true_gt_row), len(gp_row))
    true_gt_row = true_gt_row[:min_len]
    gp_row = gp_row[:min_len]
    
    # 创建混淆矩阵
    matrix = create_matrix_from_line(true_gt_row, gp_row)
    
    # 计算 concordance
    concordance = calculate_concordance(matrix)
    concordance_values.append(concordance)
    
    if (snp_idx + 1) % 200 == 0:
        print(f"  Processed {snp_idx + 1} SNPs...")

print(f"  Calculated {len(concordance_values)} concordance values")

# ============================================================
# 第五步：读取评估位点 ID 并保存结果
# ============================================================
print("\n[5/5] Saving results...")

# 读取评估位点 ID
eval_snp_ids = []
with open(EVAL_SNPS_FILE, 'r') as f:
    for line in f:
        eval_snp_ids.append(line.strip())

# 确保数量一致
if len(eval_snp_ids) < min_snps:
    print(f"  Warning: Evaluation SNP list has {len(eval_snp_ids)} SNPs, but only {min_snps} processed")
    eval_snp_ids = eval_snp_ids[:min_snps]
elif len(eval_snp_ids) > min_snps:
    eval_snp_ids = eval_snp_ids[:min_snps]

# 保存结果
output_file = f"concordance_nophasing-fold{FOLD}.txt"
with open(output_file, 'w') as f:
    for i, snp_id in enumerate(eval_snp_ids):
        f.write(f"{snp_id}\t{concordance_values[i]:.6f}\n")

print(f"  Results saved to: {output_file}")

# 统计信息
valid_values = [v for v in concordance_values if not np.isnan(v)]
if valid_values:
    # 计算统计量
    mean_val = np.mean(valid_values)
    median_val = np.median(valid_values)
    std_val = np.std(valid_values)
    min_val = np.min(valid_values)
    max_val = np.max(valid_values)
    valid_count = len(valid_values)
    total_count = len(concordance_values)
    
    print(f"\nStatistics:")
    print(f"  Mean concordance: {mean_val:.4f}")
    print(f"  Median concordance: {median_val:.4f}")
    print(f"  Std concordance: {std_val:.4f}")
    print(f"  Valid SNPs: {valid_count} / {total_count}")
    print(f"  Min concordance: {min_val:.4f}")
    print(f"  Max concordance: {max_val:.4f}")
    
    # 保存统计信息到summary文件
    summary_file = f"summary_concordance_nophasing-fold{FOLD}.txt"
    with open(summary_file, 'w') as f:
        f.write("# Summary statistics for Concordance evaluation\n")
        f.write("# =================================\n")
        f.write(f"Fold\t{FOLD}\n")
        f.write(f"Sample_Size\t{SAMPLE_SIZE}\n")
        f.write(f"Repeat\t{REPEAT}\n")
        f.write(f"Total_SNPs\t{total_count}\n")
        f.write(f"Valid_SNPs\t{valid_count}\n")
        f.write(f"Mean_Concordance\t{mean_val:.6f}\n")
        f.write(f"Median_Concordance\t{median_val:.6f}\n")
        f.write(f"Std_Concordance\t{std_val:.6f}\n")
        f.write(f"Min_Concordance\t{min_val:.6f}\n")
        f.write(f"Max_Concordance\t{max_val:.6f}\n")
        f.write("# =================================\n")
    
    print(f"  Summary saved to: {summary_file}")
    
else:
    print("\nWarning: No valid concordance values calculated")
    
    # 即使没有有效值，也保存summary文件
    summary_file = f"summary_concordance_nophasing-fold{FOLD}.txt"
    with open(summary_file, 'w') as f:
        f.write("# Summary statistics for Concordance evaluation\n")
        f.write("# =================================\n")
        f.write(f"Fold\t{FOLD}\n")
        f.write(f"Sample_Size\t{SAMPLE_SIZE}\n")
        f.write(f"Repeat\t{REPEAT}\n")
        f.write(f"Total_SNPs\t{len(concordance_values)}\n")
        f.write(f"Valid_SNPs\t0\n")
        f.write(f"Mean_Concordance\tNA\n")
        f.write(f"Median_Concordance\tNA\n")
        f.write(f"Std_Concordance\tNA\n")
        f.write(f"Min_Concordance\tNA\n")
        f.write(f"Max_Concordance\tNA\n")
        f.write("# =================================\n")
    
    print(f"  Summary saved to: {summary_file}")

print("\n" + "=" * 60)
print("Concordance evaluation completed!")
print("=" * 60)
CONCORDANCE_EOF

    # ============================================================
    # 2. 生成 R² 评估脚本（支持两种模式，采用简单顺序对齐）
    # ============================================================
    cat > ${FOLD_DIR}/IMPUTE-reference-${SAMPLE_SIZE}-R2-nophasing-for.py << 'R2_EOF'
#!/usr/bin/env python3
# ============================================================
# 脚本: IMPUTE-reference-200-R2-for.py
# 功能: 计算填充准确性 (R²) - 使用 DS 字段
#       支持两种计算模式：
#       1. 按位点计算R²后取平均（默认）
#       2. 全局R²（所有SNP和样本合并计算）
# 用法: python IMPUTE-reference-200-R2-for.py [--global-mode]
# 运行位置: 在 fold1, fold2, ... 目录下运行
# ============================================================

import gzip
import numpy as np
import os
import sys
import argparse

# ============================================================
# 参数设置
# ============================================================
SAMPLE_SIZE = 200
REPEAT = 1

# 解析命令行参数
parser = argparse.ArgumentParser(description='Calculate R² for imputation accuracy')
parser.add_argument('--global-mode', dest='global_mode', action='store_true', 
                    help='Calculate global R² (default: per-SNP R² then average)')
args = parser.parse_args()
use_global_r2 = args.global_mode

# 获取当前目录名，确定 fold 编号
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith('fold'):
    FOLD = int(current_dir.replace('fold', ''))
else:
    print("Error: This script must be run in a fold directory (fold1, fold2, ...)")
    sys.exit(1)

# 文件路径
IMPUTED_VCF = f"impute5_output_{SAMPLE_SIZE}_{REPEAT}_fold{FOLD}_real.vcf.gz"
TRUE_VCF = f"exome-population-{SAMPLE_SIZE}-{REPEAT}-fold{FOLD}-eval.vcf.gz"
EVAL_SNPS_FILE = f"/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/target-population-10random/site_fold-no-phasing/fold_{FOLD}_eval.txt"

print("=" * 60)
print(f"R² Evaluation - Fold {FOLD}")
print(f"Sample size: {SAMPLE_SIZE}")
print(f"Repeat: {REPEAT}")
print(f"Calculation mode: {'Global R²' if use_global_r2 else 'Per-SNP R² (then average)'}")
print(f"Working directory: {os.getcwd()}")
print("=" * 60)

# ============================================================
# 第一步：提取填充基因型的 DS（剂量）
# ============================================================
print("\n[1/5] Extracting DS from imputed VCF...")

imputed_ds_data = []  # 存储每行的 DS 值（浮点数列表）

with gzip.open(IMPUTED_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        if len(columns) > 9:
            # 提取所有样本的 DS 字段（第二个冒号分隔的值）
            ds_values = []
            for sample in columns[9:]:
                parts = sample.split(':')
                if len(parts) >= 2:
                    try:
                        ds = float(parts[1])  # DS 是第二个字段
                        ds_values.append(ds)
                    except:
                        ds_values.append(np.nan)
                else:
                    ds_values.append(np.nan)
            imputed_ds_data.append(ds_values)

print(f"  Found {len(imputed_ds_data)} SNPs in imputed VCF")

# ============================================================
# 第二步：提取真实基因型并转换为 0/1/2 格式
# ============================================================
print("\n[2/5] Converting true genotypes to 0/1/2 format...")

# 定义真实基因型到数字的映射
def convert_true_genotype(genotype):
    """将 VCF 格式的 GT 转换为 0/1/2/NA"""
    if genotype == './.':
        return 'NA'
    alleles = genotype.split('/')
    if alleles[0] == '0' and alleles[1] == '0':
        return '0'
    elif alleles[0] == '1' and alleles[1] == '1':
        return '2'
    else:
        return '1'

true_genotype_data = []  # 存储每行的真实基因型（0/1/2/NA）

with gzip.open(TRUE_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        if len(columns) > 9:
            # 提取所有样本的 GT 字段
            gt_values = []
            for sample in columns[9:]:
                gt = sample.split(':')[0]  # GT 是第一个字段
                converted = convert_true_genotype(gt)
                gt_values.append(converted)
            true_genotype_data.append(gt_values)

print(f"  Found {len(true_genotype_data)} SNPs in true VCF")

# ============================================================
# 第三步：确保两个文件的 SNP 数量一致（与concordance相同的方式）
# ============================================================
print("\n[3/5] Aligning SNPs...")

min_snps = min(len(imputed_ds_data), len(true_genotype_data))
if len(imputed_ds_data) != len(true_genotype_data):
    print(f"  Warning: Imputed SNPs ({len(imputed_ds_data)}) != True SNPs ({len(true_genotype_data)})")
    print(f"  Using first {min_snps} SNPs for analysis")

imputed_ds_data = imputed_ds_data[:min_snps]
true_genotype_data = true_genotype_data[:min_snps]

print(f"  Processing {min_snps} SNPs")

# ============================================================
# 第四步：根据模式计算 R²
# ============================================================
print(f"\n[4/5] Calculating R² values ({'Global mode' if use_global_r2 else 'Per-SNP mode'})...")

def calculate_pearson_r2(x, y):
    """计算皮尔逊相关系数的平方"""
    if len(x) < 2:
        return np.nan
    
    # 转换为 numpy 数组
    x = np.array(x)
    y = np.array(y)
    
    # 过滤缺失值
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_clean = x[mask]
    y_clean = y[mask]
    
    if len(x_clean) < 2:
        return np.nan
    
    # 计算相关系数的平方
    r = np.corrcoef(x_clean, y_clean)[0, 1]
    return r ** 2

if use_global_r2:
    # 全局R²模式
    print("  Using global R² calculation (all SNPs combined)...")
    
    all_ds_values = []
    all_gt_values = []
    
    for snp_idx in range(min_snps):
        imputed_ds = imputed_ds_data[snp_idx]
        true_gt = true_genotype_data[snp_idx]
        
        # 确保长度一致
        min_len = min(len(imputed_ds), len(true_gt))
        
        for i in range(min_len):
            # 跳过真实基因型缺失的样本
            if true_gt[i] == 'NA' or true_gt[i] is None:
                continue
            # 跳过填充 DS 缺失的样本
            if np.isnan(imputed_ds[i]):
                continue
            
            all_ds_values.append(imputed_ds[i])
            all_gt_values.append(float(true_gt[i]))
    
    print(f"  Total valid data points: {len(all_ds_values)}")
    
    if len(all_ds_values) >= 2:
        global_r2 = calculate_pearson_r2(all_ds_values, all_gt_values)
        print(f"  Global R² calculated: {global_r2:.6f}")
    else:
        global_r2 = np.nan
        print("  Warning: Insufficient valid data points for global R²")
    
else:
    # 按位点模式（与原来的方法相同）
    print("  Using per-SNP R² calculation (then average)...")
    
    r2_values = []
    
    for snp_idx in range(min_snps):
        # 获取该 SNP 的填充 DS 和真实基因型
        imputed_ds = imputed_ds_data[snp_idx]
        true_gt = true_genotype_data[snp_idx]
        
        # 确保长度一致
        min_len = min(len(imputed_ds), len(true_gt))
        if min_len == 0:
            r2_values.append(np.nan)
            continue
        
        # 过滤缺失值（真实基因型为 NA 的样本）
        filtered_imputed = []
        filtered_true = []
        
        for i in range(min_len):
            # 跳过真实基因型缺失的样本
            if true_gt[i] == 'NA':
                continue
            # 跳过填充 DS 缺失的样本
            if np.isnan(imputed_ds[i]):
                continue
            
            filtered_imputed.append(imputed_ds[i])
            filtered_true.append(float(true_gt[i]))
        
        # 计算 R²
        if len(filtered_imputed) >= 2:  # 至少需要2个样本才能计算相关系数
            r2 = calculate_pearson_r2(filtered_imputed, filtered_true)
        else:
            r2 = np.nan
        
        r2_values.append(r2)
        
        if (snp_idx + 1) % 200 == 0:
            print(f"  Processed {snp_idx + 1} SNPs...")
    
    print(f"  Calculated {len(r2_values)} R² values")

# ============================================================
# 第五步：读取评估位点 ID 并保存结果
# ============================================================
print("\n[5/5] Saving results...")

# 读取评估位点 ID
eval_snp_ids = []
with open(EVAL_SNPS_FILE, 'r') as f:
    for line in f:
        eval_snp_ids.append(line.strip())

if use_global_r2:
    # 全局模式：只保存汇总统计
    summary_file = f"summary_r2_global_nophasing-fold{FOLD}.txt"
    
    with open(summary_file, 'w') as f:
        f.write("# Summary statistics for R² evaluation (Global mode)\n")
        f.write("# =================================\n")
        f.write(f"Fold\t{FOLD}\n")
        f.write(f"Sample_Size\t{SAMPLE_SIZE}\n")
        f.write(f"Repeat\t{REPEAT}\n")
        f.write(f"Calculation_Mode\tGlobal\n")
        f.write(f"Total_SNPs\t{min_snps}\n")
        f.write(f"Total_Valid_Data_Points\t{len(all_ds_values) if 'all_ds_values' in dir() else 'N/A'}\n")
        f.write(f"Global_R2\t{global_r2:.6f}\n")
        f.write("# =================================\n")
    
    print(f"  Summary saved to: {summary_file}")
    print(f"\nGlobal R²: {global_r2:.6f}")
    
else:
    # 按位点模式：保存每个SNP的R²和汇总统计
    # 确保数量一致
    if len(eval_snp_ids) < min_snps:
        print(f"  Warning: Evaluation SNP list has {len(eval_snp_ids)} SNPs, but only {min_snps} processed")
        eval_snp_ids = eval_snp_ids[:min_snps]
    elif len(eval_snp_ids) > min_snps:
        eval_snp_ids = eval_snp_ids[:min_snps]
    
    # 保存结果
    output_file = f"r2_nophasing-fold{FOLD}.txt"
    with open(output_file, 'w') as f:
        for i, snp_id in enumerate(eval_snp_ids):
            f.write(f"{snp_id}\t{r2_values[i]:.6f}\n")
    
    print(f"  Results saved to: {output_file}")
    
    # 统计信息
    valid_values = [v for v in r2_values if not np.isnan(v)]
    if valid_values:
        mean_val = np.mean(valid_values)
        median_val = np.median(valid_values)
        std_val = np.std(valid_values)
        min_val = np.min(valid_values)
        max_val = np.max(valid_values)
        valid_count = len(valid_values)
        total_count = len(r2_values)
        
        print(f"\nStatistics (Per-SNP R²):")
        print(f"  Mean R²: {mean_val:.6f}")
        print(f"  Median R²: {median_val:.6f}")
        print(f"  Std R²: {std_val:.6f}")
        print(f"  Valid SNPs: {valid_count} / {total_count}")
        print(f"  Min R²: {min_val:.6f}")
        print(f"  Max R²: {max_val:.6f}")
        
        # 保存统计信息到summary文件
        summary_file = f"summary_r2_nophasing-fold{FOLD}.txt"
        with open(summary_file, 'w') as f:
            f.write("# Summary statistics for R² evaluation (Per-SNP mode)\n")
            f.write("# =================================\n")
            f.write(f"Fold\t{FOLD}\n")
            f.write(f"Sample_Size\t{SAMPLE_SIZE}\n")
            f.write(f"Repeat\t{REPEAT}\n")
            f.write(f"Calculation_Mode\tPer-SNP\n")
            f.write(f"Total_SNPs\t{total_count}\n")
            f.write(f"Valid_SNPs\t{valid_count}\n")
            f.write(f"Mean_R2\t{mean_val:.6f}\n")
            f.write(f"Median_R2\t{median_val:.6f}\n")
            f.write(f"Std_R2\t{std_val:.6f}\n")
            f.write(f"Min_R2\t{min_val:.6f}\n")
            f.write(f"Max_R2\t{max_val:.6f}\n")
            f.write("# =================================\n")
        
        print(f"  Summary saved to: {summary_file}")
        
    else:
        print("\nWarning: No valid R² values calculated")
        
        # 即使没有有效值，也保存summary文件
        summary_file = f"summary_r2_nophasing-fold{FOLD}.txt"
        with open(summary_file, 'w') as f:
            f.write("# Summary statistics for R² evaluation (Per-SNP mode)\n")
            f.write("# =================================\n")
            f.write(f"Fold\t{FOLD}\n")
            f.write(f"Sample_Size\t{SAMPLE_SIZE}\n")
            f.write(f"Repeat\t{REPEAT}\n")
            f.write(f"Calculation_Mode\tPer-SNP\n")
            f.write(f"Total_SNPs\t{len(r2_values)}\n")
            f.write(f"Valid_SNPs\t0\n")
            f.write(f"Mean_R2\tNA\n")
            f.write(f"Median_R2\tNA\n")
            f.write(f"Std_R2\tNA\n")
            f.write(f"Min_R2\tNA\n")
            f.write(f"Max_R2\tNA\n")
            f.write("# =================================\n")
        
        print(f"  Summary saved to: {summary_file}")

print("\n" + "=" * 60)
print("R² evaluation completed!")
print("=" * 60)
R2_EOF

    # ============================================================
    # 3. 赋予执行权限并运行脚本（运行两种模式）
    # ============================================================
    chmod +x ${FOLD_DIR}/IMPUTE-reference-${SAMPLE_SIZE}-concordance-nophasing-for.py
    chmod +x ${FOLD_DIR}/IMPUTE-reference-${SAMPLE_SIZE}-R2-nophasing-for.py
    
    echo ""
    echo "Running Concordance evaluation for fold ${fold}..."
    cd ${FOLD_DIR}
    python IMPUTE-reference-${SAMPLE_SIZE}-concordance-nophasing-for.py
    
    echo ""
    echo "Running R² evaluation for fold ${fold} (Per-SNP mode)..."
    python IMPUTE-reference-${SAMPLE_SIZE}-R2-nophasing-for.py
    
    echo ""
    echo "Running R² evaluation for fold ${fold} (Global mode)..."
    python IMPUTE-reference-${SAMPLE_SIZE}-R2-nophasing-for.py --global-mode
    
    # 返回基础目录
    cd ${BASE_DIR}
    
    echo ""
    echo "Fold ${fold} completed"
    echo "=========================================="

done

echo ""
echo "=========================================="
echo "All folds processed!"
echo "Results saved in each fold directory:"
for fold in {1..10}; do
    echo "  fold${fold}/concordance_nophasing-fold${fold}.txt"
    echo "  fold${fold}/r2_nophasing-fold${fold}.txt (Per-SNP)"
    echo "  fold${fold}/summary_concordance_nophasing-fold${fold}.txt"
    echo "  fold${fold}/summary_r2_nophasing-fold${fold}.txt (Per-SNP stats)"
    echo "  fold${fold}/summary_r2_global_nophasing-fold${fold}.txt (Global R²)"
done
echo "=========================================="

