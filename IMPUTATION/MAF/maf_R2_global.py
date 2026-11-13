#!/usr/bin/env python3
"""
MAF-stratified R² Calculator - Memory Efficient (Incremental Algorithm)
自动匹配真实VCF和填充VCF的位点，按MAF区间计算全局R²
内存占用恒定，可处理任意大小数据
"""

import gzip
import numpy as np
import argparse
import sys
from math import sqrt


class IncrementalCorrelation:
    """增量计算相关系数的类 - 不需要存储任何数据点"""
    
    def __init__(self):
        self.n = 0          # 数据点数量
        self.sum_x = 0.0    # DS值的总和
        self.sum_y = 0.0    # GT值的总和
        self.sum_xy = 0.0   # DS * GT 的总和
        self.sum_x2 = 0.0   # DS² 的总和
        self.sum_y2 = 0.0   # GT² 的总和
    
    def add(self, x, y):
        """添加一个数据点（增量更新）"""
        self.n += 1
        self.sum_x += x
        self.sum_y += y
        self.sum_xy += x * y
        self.sum_x2 += x * x
        self.sum_y2 += y * y
    
    def get_r2(self):
        """计算R²（皮尔逊相关系数的平方）"""
        if self.n < 3:
            return np.nan
        
        # 计算相关系数
        numerator = self.n * self.sum_xy - self.sum_x * self.sum_y
        
        denominator_x = self.n * self.sum_x2 - self.sum_x * self.sum_x
        denominator_y = self.n * self.sum_y2 - self.sum_y * self.sum_y
        
        # 避免除零错误
        if denominator_x <= 0 or denominator_y <= 0:
            return np.nan
        
        denominator = sqrt(denominator_x * denominator_y)
        
        if denominator == 0:
            return np.nan
        
        r = numerator / denominator
        return r * r  # R² = r²
    
    def get_stats(self):
        """返回当前统计信息"""
        return {
            'n': self.n,
            'sum_x': self.sum_x,
            'sum_y': self.sum_y,
            'mean_x': self.sum_x / self.n if self.n > 0 else 0,
            'mean_y': self.sum_y / self.n if self.n > 0 else 0
        }


def get_sample_count(vcf_path):
    """快速获取VCF文件的样本数量"""
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                return len(line.strip().split('\t')) - 9
    return 0


def get_vcf_samples(vcf_path):
    """获取VCF文件的样本列表"""
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                return line.strip().split('\t')[9:]
    return []


def convert_true_genotype(gt):
    """将真实基因型GT转换为数值 0/1/2"""
    if gt == '0/0':
        return 0.0
    elif gt == '0/1' or gt == '1/0':
        return 1.0
    elif gt == '1/1':
        return 2.0
    else:
        return np.nan


def extract_dosage(sample_data):
    """从填充VCF的样本数据中提取DS（剂量）值"""
    parts = sample_data.split(':')
    
    # 尝试多种格式
    if len(parts) >= 2:
        # 格式1: GT:DS:GP 或 GT:DS
        try:
            if 'DS' in parts:
                ds_idx = parts.index('DS') + 1
                if ds_idx < len(parts):
                    return float(parts[ds_idx])
            else:
                # 假设第二个字段是DS
                return float(parts[1])
        except (ValueError, IndexError):
            pass
    
    return np.nan


def extract_maf_from_info(info_str):
    """
    从VCF的INFO字段提取MAF值（次要等位基因频率）
    
    支持的格式：
    - MAF=0.05;AN=100;AC=5  → 直接提取MAF值
    - AF=0.05;AN=100;AC=5   → 提取AF后转换为MAF = min(AF, 1-AF)
    
    Parameters:
    -----------
    info_str : str
        VCF文件中INFO字段的内容，例如 "MAF=0.05;AN=100;AC=5"
    
    Returns:
    --------
    float or None
        返回MAF值（范围0-0.5），如果无法提取则返回None
    """
    if not info_str:
        return None
    
    fields = info_str.split(';')
    for field in fields:
        # 优先查找MAF字段
        if field.startswith('MAF='):
            try:
                maf = float(field.split('=')[1])
                # MAF应该在0-0.5之间
                if 0 <= maf <= 0.5:
                    return maf
                else:
                    return None
            except (ValueError, IndexError):
                continue
        
        # 如果没有MAF字段，尝试AF字段
        elif field.startswith('AF='):
            try:
                af = float(field.split('=')[1])
                # 等位基因频率，转换为MAF
                maf = min(af, 1 - af)
                if 0 <= maf <= 0.5:
                    return maf
                else:
                    return None
            except (ValueError, IndexError):
                continue
    
    return None


class MAFStratifiedCalculator:
    """
    按MAF区间分层计算R²的类
    
    功能：
    1. 根据MAF值将SNP分配到不同的MAF区间
    2. 在每个区间内增量计算相关系数
    3. 统计每个区间的SNP数量和数据点数量
    
    Attributes:
    -----------
    maf_bins : list of tuples
        MAF区间列表，每个元素为(min_maf, max_maf)
    correlations : dict
        存储每个区间的增量计算器对象
    snp_counts : dict
        存储每个区间的SNP数量
    total_snps_processed : int
        总处理的SNP数量
    """
    
    def __init__(self, maf_bins):
        """
        初始化分层计算器
        
        Parameters:
        -----------
        maf_bins : list of tuples
            MAF区间列表，每个元素为(min_maf, max_maf)
            区间规则：左闭右开 [min, max)
            最后一个区间可以设置max=1.0来包含MAF=0.5的情况
        """
        self.maf_bins = maf_bins
        # 为每个MAF区间创建一个独立的增量计算器
        self.correlations = {f"{low:.3f}-{high:.3f}": IncrementalCorrelation() 
                            for low, high in maf_bins}
        # 记录每个区间包含的SNP数量
        self.snp_counts = {f"{low:.3f}-{high:.3f}": 0 for low, high in maf_bins}
        self.total_snps_processed = 0
    
    def add_snp(self, maf, ds_values, gt_values):
        """
        添加一个SNP的数据到对应的MAF区间
        
        Parameters:
        -----------
        maf : float
            该SNP的MAF值
        ds_values : list
            填充数据的DS值列表
        gt_values : list
            真实数据的GT值列表（0/1/2格式）
        
        Returns:
        --------
        int
            该SNP的有效数据点数量（非NaN的样本数）
        """
        # 找到对应的MAF区间
        bin_key = None
        for low, high in self.maf_bins:
            # 左闭右开区间，最后一个区间特殊处理
            if low <= maf < high or (high == 1.0 and maf >= low):
                bin_key = f"{low:.3f}-{high:.3f}"
                break
        
        if bin_key is None:
            return 0  # MAF不在任何区间内
        
        # 添加所有样本数据点到对应的计算器
        corr = self.correlations[bin_key]
        valid_points = 0
        
        for ds, gt in zip(ds_values, gt_values):
            if not np.isnan(ds) and not np.isnan(gt):
                corr.add(ds, gt)
                valid_points += 1
        
        # 只有至少有1个有效数据点才计数这个SNP
        if valid_points > 0:
            self.snp_counts[bin_key] += 1
            self.total_snps_processed += 1
        
        return valid_points
    
    def get_results(self):
        """
        获取所有区间的R²计算结果
        
        Returns:
        --------
        dict
            键：MAF区间名称（如"0.000-0.010"）
            值：包含r2、n_snps、n_data_points、mean_ds、mean_gt的字典
        """
        results = {}
        for bin_key in self.correlations:
            r2 = self.correlations[bin_key].get_r2()
            stats = self.correlations[bin_key].get_stats()
            results[bin_key] = {
                'r2': r2,                      # 该区间的R²值
                'n_snps': self.snp_counts[bin_key],  # 该区间的SNP数量
                'n_data_points': stats['n'],   # 该区间总数据点数量
                'mean_ds': stats['mean_x'],    # 填充DS均值
                'mean_gt': stats['mean_y']     # 真实GT均值
            }
        return results


def calculate_maf_stratified_r2(true_vcf, imputed_vcf, maf_bins, verbose=True):
    """
    计算按MAF分层的全局R² - 自动匹配位点，增量计算
    
    Parameters:
    -----------
    true_vcf : str
        真实基因型VCF文件路径（必须包含MAF或AF信息在INFO字段）
    imputed_vcf : str
        填充VCF文件路径（必须包含DS剂量值）
    maf_bins : list of tuples
        MAF区间列表
    verbose : bool
        是否输出详细信息
    
    Returns:
    --------
    dict : 各MAF区间的R²结果
    """
    
    print("=" * 70)
    print("MAF-Stratified R² Calculator (Incremental Algorithm)")
    print("=" * 70)
    print(f"True VCF:     {true_vcf}")
    print(f"Imputed VCF:  {imputed_vcf}")
    print(f"MAF bins:     {len(maf_bins)} intervals")
    print("=" * 70)
    
    # 获取样本信息
    print("\n[1/5] Checking VCF files...")
    
    true_samples = get_vcf_samples(true_vcf)
    imputed_samples = get_vcf_samples(imputed_vcf)
    
    print(f"  True VCF samples:     {len(true_samples)}")
    print(f"  Imputed VCF samples:  {len(imputed_samples)}")
    
    # 验证样本顺序
    if len(true_samples) != len(imputed_samples):
        print(f"  ERROR: Sample count mismatch! {len(true_samples)} vs {len(imputed_samples)}")
        sys.exit(1)
    
    n_samples = len(true_samples)
    print(f"  Number of samples to analyze: {n_samples}")
    
    # 初始化分层计算器
    print("\n[2/5] Initializing MAF-stratified calculator...")
    calculator = MAFStratifiedCalculator(maf_bins)
    print(f"  Created {len(maf_bins)} MAF bins")
    
    # 打开两个VCF文件
    print("\n[3/5] Streaming VCF files and matching SNPs...")
    
    true_f = gzip.open(true_vcf, 'rt')
    imputed_f = gzip.open(imputed_vcf, 'rt')
    
    # 跳过标题行，找到第一个数据行
    true_line = None
    imputed_line = None
    
    for line in true_f:
        if not line.startswith('#'):
            true_line = line
            break
    
    for line in imputed_f:
        if not line.startswith('#'):
            imputed_line = line
            break
    
    # 统计信息
    total_snps = 0
    matched_snps = 0
    skipped_snps = 0
    total_data_points = 0
    no_maf_snps = 0  # 没有MAF信息的SNP数量
    
    # 逐位点处理
    while true_line and imputed_line:
        total_snps += 1
        
        # 解析真实VCF的当前位置
        true_cols = true_line.strip().split('\t')
        true_chrom = true_cols[0].replace('chr', '')
        true_pos = true_cols[1]
        true_snp_id = f"{true_chrom}:{true_pos}"
        
        # 解析填充VCF的当前位置
        imp_cols = imputed_line.strip().split('\t')
        imp_chrom = imp_cols[0].replace('chr', '')
        imp_pos = imp_cols[1]
        imp_snp_id = f"{imp_chrom}:{imp_pos}"
        
        # 比较位点位置
        if true_snp_id < imp_snp_id:
            # 真实VCF的位点在填充VCF中不存在，跳过
            true_line = true_f.readline()
            skipped_snps += 1
            continue
        elif true_snp_id > imp_snp_id:
            # 填充VCF的位点在真实VCF中不存在，跳过
            imputed_line = imputed_f.readline()
            skipped_snps += 1
            continue
        
        # 位点匹配成功！
        matched_snps += 1
        
        # 提取MAF值（从真实VCF的INFO字段）
        # INFO字段通常是第8列（索引7）
        info_fields = true_cols[7] if len(true_cols) > 7 else ""
        maf = extract_maf_from_info(info_fields)
        
        if maf is None or np.isnan(maf):
            # 如果无法提取MAF，跳过该位点
            no_maf_snps += 1
            true_line = true_f.readline()
            imputed_line = imputed_f.readline()
            continue
        
        # 提取真实基因型（GT字段）
        true_gt_values = []
        for i in range(n_samples):
            if i < len(true_cols[9:]):
                gt = true_cols[9+i].split(':')[0]
                gt_val = convert_true_genotype(gt)
                true_gt_values.append(gt_val)
            else:
                true_gt_values.append(np.nan)
        
        # 提取填充剂量（DS字段）
        imp_ds_values = []
        for i in range(n_samples):
            if i < len(imp_cols[9:]):
                ds_val = extract_dosage(imp_cols[9+i])
                imp_ds_values.append(ds_val)
            else:
                imp_ds_values.append(np.nan)
        
        # 添加到分层计算器
        valid_points = calculator.add_snp(maf, imp_ds_values, true_gt_values)
        total_data_points += valid_points
        
        # 读取下一行
        true_line = true_f.readline()
        imputed_line = imputed_f.readline()
        
        # 进度显示
        if verbose and matched_snps % 10000 == 0:
            print(f"    Processed {matched_snps} matched SNPs, "
                  f"{total_data_points} data points, "
                  f"skipped {skipped_snps} SNPs, "
                  f"no MAF: {no_maf_snps}")
    
    # 关闭文件
    true_f.close()
    imputed_f.close()
    
    # 输出统计信息
    print(f"\n  ✓ Processing complete!")
    print(f"    Total SNPs processed:    {total_snps}")
    print(f"    Matched SNPs:            {matched_snps}")
    print(f"    Unmatched SNPs:          {skipped_snps}")
    print(f"    SNPs without MAF info:   {no_maf_snps}")
    print(f"    Total data points:       {total_data_points}")
    
    # 计算所有区间的R²
    print("\n[4/5] Calculating R² for each MAF bin...")
    results = calculator.get_results()
    
    # 输出结果
    print("\n[5/5] Results")
    print("=" * 70)
    print(f"{'MAF Bin':<20} {'SNPs':<10} {'Data Points':<15} {'R²':<12}")
    print("-" * 70)
    
    for bin_key in sorted(results.keys()):
        res = results[bin_key]
        if not np.isnan(res['r2']):
            print(f"{bin_key:<20} {res['n_snps']:<10} {res['n_data_points']:<15} {res['r2']:.6f}")
        else:
            print(f"{bin_key:<20} {res['n_snps']:<10} {res['n_data_points']:<15} {'N/A':<12}")
    
    print("=" * 70)
    
    return results


def get_default_maf_bins():
    """返回默认的MAF区间（左闭右开）"""
    return [
        (0.00, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04),
        (0.04, 0.05), (0.05, 0.06), (0.06, 0.07), (0.07, 0.08),
        (0.08, 0.09), (0.09, 0.10), (0.10, 0.15), (0.15, 0.20),
        (0.20, 0.25), (0.25, 0.30), (0.30, 0.35), (0.35, 0.40),
        (0.40, 0.45), (0.45, 0.51)  # 0.51确保包含MAF=0.5
    ]


def main():
    parser = argparse.ArgumentParser(
        description='Calculate MAF-stratified R² between true and imputed VCF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 基本使用（使用默认MAF区间）
  python maf_stratified_r2.py -t real.vcf.gz -i imputed.vcf.gz -o results.txt
  
  # 自定义MAF区间
  python maf_stratified_r2.py -t real.vcf.gz -i imputed.vcf.gz \\
      --maf-bins 0.00,0.01 0.01,0.05 0.05,0.10 0.10,0.50 -o results.txt
  
  # 静默模式
  python maf_stratified_r2.py -t real.vcf.gz -i imputed.vcf.gz -q
  
Note:
  真实VCF文件必须在INFO字段中包含MAF或AF信息，例如：
    MAF=0.05;AN=100;AC=5
  或
    AF=0.05;AN=100;AC=5
        """
    )
    
    parser.add_argument('-t', '--true-vcf', required=True,
                        help='True genotype VCF file (must contain MAF or AF in INFO field)')
    parser.add_argument('-i', '--imputed-vcf', required=True,
                        help='Imputed VCF file (must contain DS dosage values)')
    parser.add_argument('-o', '--output', default='maf_stratified_R2.txt',
                        help='Output file (default: maf_stratified_R2.txt)')
    parser.add_argument('--maf-bins', nargs='+', action='append',
                        help='MAF bins in format "min,max" (can specify multiple)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress detailed progress output')
    
    args = parser.parse_args()
    
    # 解析MAF区间
    if args.maf_bins:
        maf_bins = []
        for bin_spec in args.maf_bins[0]:  # 因为action='append'会创建嵌套列表
            parts = bin_spec.split(',')
            if len(parts) == 2:
                maf_bins.append((float(parts[0]), float(parts[1])))
            else:
                print(f"ERROR: Invalid MAF bin format: {bin_spec}")
                sys.exit(1)
    else:
        maf_bins = get_default_maf_bins()
    
    # 计算分层R²
    results = calculate_maf_stratified_r2(
        args.true_vcf, 
        args.imputed_vcf, 
        maf_bins,
        verbose=not args.quiet
    )
    
    if results is None:
        sys.exit(1)
    
    # 保存结果
    with open(args.output, 'w') as f:
        f.write("# MAF-Stratified R² Calculation Results\n")
        f.write("# ========================================\n")
        f.write(f"# True VCF:      {args.true_vcf}\n")
        f.write(f"# Imputed VCF:   {args.imputed_vcf}\n")
        f.write("# ========================================\n")
        f.write("MAF_Bin\tSNPs_Count\tData_Points\tMean_DS\tMean_GT\tR2\n")
        
        for bin_key in sorted(results.keys()):
            res = results[bin_key]
            f.write(f"{bin_key}\t{res['n_snps']}\t{res['n_data_points']}\t"
                   f"{res['mean_ds']:.6f}\t{res['mean_gt']:.6f}\t")
            if not np.isnan(res['r2']):
                f.write(f"{res['r2']:.10f}\n")
            else:
                f.write("NA\n")
        
        f.write("# ========================================\n")
    
    print(f"\n✓ Results saved to: {args.output}")


if __name__ == "__main__":
    main()
