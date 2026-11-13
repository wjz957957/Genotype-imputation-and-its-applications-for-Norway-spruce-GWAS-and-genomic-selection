#!/usr/bin/env python3
"""
MAF-stratified Concordance Calculator - Memory Efficient (Incremental Algorithm)
自动匹配真实VCF和填充VCF的位点，按MAF区间计算基因型一致性
内存占用恒定，可处理任意大小数据
"""

import gzip
import numpy as np
import argparse
import sys
from collections import defaultdict


class IncrementalConcordance:
    """
    增量计算Concordance的类
    每个MAF区间独立维护一个concordance累加器
    """
    
    def __init__(self):
        # 混淆矩阵: 3x3 (真实基因型 0,1,2 vs 预测基因型 0,1,2)
        self.matrix = np.zeros((3, 3), dtype=np.float64)
        self.total_samples = 0
    
    def add_sample(self, true_gt, gp_probs):
        """
        添加一个样本的concordance贡献
        
        Parameters:
        -----------
        true_gt : int
            真实基因型 (0, 1, 2)
        gp_probs : list of float
            基因型概率 [P(0/0), P(0/1), P(1/1)]
        """
        if true_gt == 'NA' or true_gt is None:
            return
        
        true_idx = int(true_gt)
        
        # 对于每个可能的预测基因型，添加概率权重
        for pred_idx in range(3):
            if not np.isnan(gp_probs[pred_idx]):
                self.matrix[pred_idx, true_idx] += gp_probs[pred_idx]
        
        self.total_samples += 1
    
    def get_concordance(self):
        """计算加权一致性"""
        if self.total_samples == 0:
            return np.nan
        
        diagonal_sum = np.sum(np.diag(self.matrix))
        total = np.sum(self.matrix)
        
        if total > 0:
            return diagonal_sum / total
        else:
            return np.nan
    
    def get_matrix(self):
        """返回混淆矩阵"""
        return self.matrix
    
    def get_stats(self):
        """返回统计信息"""
        concordance = self.get_concordance()
        return {
            'concordance': concordance,
            'total_weighted_samples': self.total_samples,
            'matrix': self.matrix
        }


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
        return 0
    elif gt == '0/1' or gt == '1/0':
        return 1
    elif gt == '1/1':
        return 2
    else:
        return 'NA'


def extract_gp(sample_data):
    """从填充VCF的样本数据中提取GP（基因型概率）"""
    parts = sample_data.split(':')
    
    gp = [np.nan, np.nan, np.nan]
    
    # 查找GP字段
    if 'GP' in parts:
        gp_idx = parts.index('GP')
        if gp_idx + 1 < len(parts):
            try:
                gp = [float(x) for x in parts[gp_idx+1].split(',')]
            except:
                pass
    elif len(parts) >= 3:
        # 假设第三部分是GP
        try:
            gp = [float(x) for x in parts[2].split(',')]
        except:
            pass
    elif len(parts) >= 2:
        # 可能只有GT:GP格式
        try:
            gp = [float(x) for x in parts[1].split(',')]
        except:
            pass
    
    return gp


def extract_maf_from_info(info_str):
    """
    从VCF的INFO字段提取MAF值（次要等位基因频率）
    
    支持的格式：
    - MAF=0.05;AN=100;AC=5  → 直接提取MAF值
    - AF=0.05;AN=100;AC=5   → 提取AF后转换为MAF = min(AF, 1-AF)
    """
    if not info_str:
        return None
    
    fields = info_str.split(';')
    for field in fields:
        if field.startswith('MAF='):
            try:
                maf = float(field.split('=')[1])
                if 0 <= maf <= 0.5:
                    return maf
            except:
                continue
        elif field.startswith('AF='):
            try:
                af = float(field.split('=')[1])
                maf = min(af, 1 - af)
                if 0 <= maf <= 0.5:
                    return maf
            except:
                continue
    
    return None


class MAFStratifiedConcordanceCalculator:
    """
    按MAF区间分层计算Concordance的类
    
    功能：
    1. 根据MAF值将SNP分配到不同的MAF区间
    2. 在每个区间内增量计算concordance
    3. 统计每个区间的SNP数量
    """
    
    def __init__(self, maf_bins):
        """
        初始化分层计算器
        
        Parameters:
        -----------
        maf_bins : list of tuples
            MAF区间列表，每个元素为(min_maf, max_maf)
            区间规则：左闭右开 [min, max)
        """
        self.maf_bins = maf_bins
        # 为每个MAF区间创建一个独立的concordance计算器
        self.concordance_calcs = {f"{low:.3f}-{high:.3f}": IncrementalConcordance() 
                                   for low, high in maf_bins}
        # 记录每个区间包含的SNP数量
        self.snp_counts = {f"{low:.3f}-{high:.3f}": 0 for low, high in maf_bins}
        self.total_snps_processed = 0
    
    def add_snp(self, maf, gp_values, gt_values):
        """
        添加一个SNP的数据到对应的MAF区间
        
        Parameters:
        -----------
        maf : float
            该SNP的MAF值
        gp_values : list
            填充数据的GP值列表，每个元素为[P0, P1, P2]
        gt_values : list
            真实数据的GT值列表（0/1/2格式）
        
        Returns:
        --------
        int
            该SNP的有效样本数量
        """
        # 找到对应的MAF区间
        bin_key = None
        for low, high in self.maf_bins:
            if low <= maf < high or (high == 1.0 and maf >= low):
                bin_key = f"{low:.3f}-{high:.3f}"
                break
        
        if bin_key is None:
            return 0  # MAF不在任何区间内
        
        # 添加所有样本数据点到对应的计算器
        conc_calc = self.concordance_calcs[bin_key]
        valid_samples = 0
        
        for gp, gt in zip(gp_values, gt_values):
            if gt != 'NA' and not np.isnan(gp[0]):
                conc_calc.add_sample(gt, gp)
                valid_samples += 1
        
        # 只要有至少1个有效样本就计数这个SNP
        if valid_samples > 0:
            self.snp_counts[bin_key] += 1
            self.total_snps_processed += 1
        
        return valid_samples
    
    def get_results(self):
        """
        获取所有区间的concordance计算结果
        
        Returns:
        --------
        dict
            键：MAF区间名称
            值：包含concordance、n_snps、n_samples、matrix的字典
        """
        results = {}
        for bin_key in self.concordance_calcs:
            stats = self.concordance_calcs[bin_key].get_stats()
            results[bin_key] = {
                'concordance': stats['concordance'],
                'n_snps': self.snp_counts[bin_key],
                'n_weighted_samples': stats['total_weighted_samples'],
                'confusion_matrix': stats['matrix']
            }
        return results


def calculate_maf_stratified_concordance(true_vcf, imputed_vcf, maf_bins, 
                                          region=None, verbose=True):
    """
    计算按MAF分层的Concordance - 自动匹配位点，增量计算
    
    Parameters:
    -----------
    true_vcf : str
        真实基因型VCF文件路径（必须包含MAF或AF信息）
    imputed_vcf : str
        填充VCF文件路径（必须包含GP概率值）
    maf_bins : list of tuples
        MAF区间列表
    region : str
        可选，分析区域（如 "11:83036-100000000"）
    verbose : bool
        是否输出详细信息
    
    Returns:
    --------
    dict : 各MAF区间的concordance结果
    """
    
    print("=" * 70)
    print("MAF-Stratified Concordance Calculator (Incremental Algorithm)")
    print("=" * 70)
    print(f"True VCF:     {true_vcf}")
    print(f"Imputed VCF:  {imputed_vcf}")
    print(f"MAF bins:     {len(maf_bins)} intervals")
    if region:
        print(f"Region:       {region}")
    print("=" * 70)
    
    # 获取样本信息
    print("\n[1/5] Checking VCF files...")
    
    true_samples = get_vcf_samples(true_vcf)
    imputed_samples = get_vcf_samples(imputed_vcf)
    
    print(f"  True VCF samples:     {len(true_samples)}")
    print(f"  Imputed VCF samples:  {len(imputed_samples)}")
    
    # 验证样本数量
    if len(true_samples) != len(imputed_samples):
        print(f"  ERROR: Sample count mismatch!")
        sys.exit(1)
    
    n_samples = len(true_samples)
    print(f"  Number of samples:    {n_samples}")
    
    # 解析区域过滤
    region_chrom = None
    region_start = None
    region_end = None
    if region:
        try:
            if ':' in region:
                region_chrom, pos_range = region.split(':')
                if '-' in pos_range:
                    region_start, region_end = map(int, pos_range.split('-'))
                else:
                    region_start = int(pos_range)
                    region_end = int(pos_range)
            else:
                region_chrom = region
            print(f"  Region filter:        {region_chrom}:{region_start}-{region_end}" if region_start else f"  Region filter:        {region_chrom}")
        except:
            print(f"  Warning: Invalid region format '{region}', ignoring filter")
    
    # 初始化分层计算器
    print("\n[2/5] Initializing MAF-stratified calculator...")
    calculator = MAFStratifiedConcordanceCalculator(maf_bins)
    print(f"  Created {len(maf_bins)} MAF bins")
    
    # 打开两个VCF文件
    print("\n[3/5] Streaming VCF files and matching SNPs...")
    
    true_f = gzip.open(true_vcf, 'rt')
    imputed_f = gzip.open(imputed_vcf, 'rt')
    
    # 跳过标题行
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
    filtered_snps = 0
    total_samples_processed = 0
    no_maf_snps = 0
    
    # 逐位点处理
    while true_line and imputed_line:
        total_snps += 1
        
        # 解析真实VCF
        true_cols = true_line.strip().split('\t')
        true_chrom = true_cols[0].replace('chr', '')
        true_pos = true_cols[1]
        true_snp_id = f"{true_chrom}:{true_pos}"
        
        # 解析填充VCF
        imp_cols = imputed_line.strip().split('\t')
        imp_chrom = imp_cols[0].replace('chr', '')
        imp_pos = imp_cols[1]
        imp_snp_id = f"{imp_chrom}:{imp_pos}"
        
        # 区域过滤
        if region_chrom:
            if true_chrom != region_chrom:
                true_line = true_f.readline()
                filtered_snps += 1
                continue
            if region_start and (int(true_pos) < region_start or int(true_pos) > region_end):
                true_line = true_f.readline()
                filtered_snps += 1
                continue
        
        # 比较位点位置
        if true_snp_id < imp_snp_id:
            true_line = true_f.readline()
            skipped_snps += 1
            continue
        elif true_snp_id > imp_snp_id:
            imputed_line = imputed_f.readline()
            skipped_snps += 1
            continue
        
        # 位点匹配成功
        matched_snps += 1
        
        # 提取MAF值
        info_fields = true_cols[7] if len(true_cols) > 7 else ""
        maf = extract_maf_from_info(info_fields)
        
        if maf is None or np.isnan(maf):
            no_maf_snps += 1
            true_line = true_f.readline()
            imputed_line = imputed_f.readline()
            continue
        
        # 提取真实基因型
        true_gt_values = []
        for i in range(n_samples):
            if i < len(true_cols[9:]):
                gt = true_cols[9+i].split(':')[0]
                gt_val = convert_true_genotype(gt)
                true_gt_values.append(gt_val)
            else:
                true_gt_values.append('NA')
        
        # 提取填充GP
        imp_gp_values = []
        for i in range(n_samples):
            if i < len(imp_cols[9:]):
                gp = extract_gp(imp_cols[9+i])
                imp_gp_values.append(gp)
            else:
                imp_gp_values.append([np.nan, np.nan, np.nan])
        
        # 添加到分层计算器
        valid_samples = calculator.add_snp(maf, imp_gp_values, true_gt_values)
        total_samples_processed += valid_samples
        
        # 读取下一行
        true_line = true_f.readline()
        imputed_line = imputed_f.readline()
        
        # 进度显示
        if verbose and matched_snps % 10000 == 0:
            print(f"    Processed {matched_snps} matched SNPs, "
                  f"{total_samples_processed} samples, "
                  f"skipped {skipped_snps}, filtered {filtered_snps}")
    
    # 关闭文件
    true_f.close()
    imputed_f.close()
    
    # 输出统计信息
    print(f"\n  ✓ Processing complete!")
    print(f"    Total SNPs scanned:      {total_snps}")
    print(f"    Matched SNPs:            {matched_snps}")
    print(f"    Unmatched SNPs:          {skipped_snps}")
    print(f"    Filtered by region:      {filtered_snps}")
    print(f"    SNPs without MAF info:   {no_maf_snps}")
    print(f"    Total weighted samples:  {total_samples_processed}")
    
    # 计算所有区间的Concordance
    print("\n[4/5] Calculating concordance for each MAF bin...")
    results = calculator.get_results()
    
    # 输出结果
    print("\n[5/5] Results")
    print("=" * 80)
    print(f"{'MAF Bin':<20} {'SNPs':<10} {'Weighted Samples':<18} {'Concordance':<12}")
    print("-" * 80)
    
    for bin_key in sorted(results.keys()):
        res = results[bin_key]
        if not np.isnan(res['concordance']):
            print(f"{bin_key:<20} {res['n_snps']:<10,} {res['n_weighted_samples']:<18,} {res['concordance']:.6f}")
        else:
            print(f"{bin_key:<20} {res['n_snps']:<10,} {res['n_weighted_samples']:<18,} {'N/A':<12}")
    
    print("=" * 80)
    
    return results


def get_default_maf_bins():
    """返回默认的MAF区间（左闭右开）"""
    return [
        (0.00, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04),
        (0.04, 0.05), (0.05, 0.06), (0.06, 0.07), (0.07, 0.08),
        (0.08, 0.09), (0.09, 0.10), (0.10, 0.15), (0.15, 0.20),
        (0.20, 0.25), (0.25, 0.30), (0.30, 0.35), (0.35, 0.40),
        (0.40, 0.45), (0.45, 0.51)
    ]


def main():
    parser = argparse.ArgumentParser(
        description='Calculate MAF-stratified genotype concordance between true and imputed VCF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 基本使用（使用默认MAF区间）
  python maf_concordance.py -t real.vcf.gz -i imputed.vcf.gz -o results.txt
  
  # 指定区域
  python maf_concordance.py -t real.vcf.gz -i imputed.vcf.gz -r "11:83036-100000000"
  
  # 自定义MAF区间
  python maf_concordance.py -t real.vcf.gz -i imputed.vcf.gz \\
      --maf-bins 0.00,0.01 0.01,0.05 0.05,0.10 0.10,0.50 -o results.txt
  
  # 静默模式
  python maf_concordance.py -t real.vcf.gz -i imputed.vcf.gz -q

Note:
  真实VCF文件必须在INFO字段中包含MAF或AF信息
  填充VCF文件必须包含GP（基因型概率）字段
        """
    )
    
    parser.add_argument('-t', '--true-vcf', required=True,
                        help='True genotype VCF file (must contain MAF or AF in INFO field)')
    parser.add_argument('-i', '--imputed-vcf', required=True,
                        help='Imputed VCF file (must contain GP probability values)')
    parser.add_argument('-o', '--output', default='maf_stratified_concordance.txt',
                        help='Output file (default: maf_stratified_concordance.txt)')
    parser.add_argument('-r', '--region',
                        help='Region to analyze (e.g., "11:83036-100000000")')
    parser.add_argument('--maf-bins', nargs='+', action='append',
                        help='MAF bins in format "min,max" (can specify multiple)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress detailed progress output')
    
    args = parser.parse_args()
    
    # 解析MAF区间
    if args.maf_bins:
        maf_bins = []
        for bin_spec in args.maf_bins[0]:
            parts = bin_spec.split(',')
            if len(parts) == 2:
                maf_bins.append((float(parts[0]), float(parts[1])))
            else:
                print(f"ERROR: Invalid MAF bin format: {bin_spec}")
                sys.exit(1)
    else:
        maf_bins = get_default_maf_bins()
    
    # 计算分层Concordance
    results = calculate_maf_stratified_concordance(
        args.true_vcf, 
        args.imputed_vcf, 
        maf_bins,
        region=args.region,
        verbose=not args.quiet
    )
    
    if results is None:
        sys.exit(1)
    
    # 保存结果
    with open(args.output, 'w') as f:
        f.write("# MAF-Stratified Concordance Calculation Results\n")
        f.write("# ================================================\n")
        f.write(f"# True VCF:      {args.true_vcf}\n")
        f.write(f"# Imputed VCF:   {args.imputed_vcf}\n")
        if args.region:
            f.write(f"# Region:        {args.region}\n")
        f.write("# ================================================\n")
        f.write("MAF_Bin\tSNPs_Count\tWeighted_Samples\tConcordance\n")
        
        for bin_key in sorted(results.keys()):
            res = results[bin_key]
            f.write(f"{bin_key}\t{res['n_snps']}\t{res['n_weighted_samples']}\t")
            if not np.isnan(res['concordance']):
                f.write(f"{res['concordance']:.10f}\n")
            else:
                f.write("NA\n")
        
        # 可选：保存混淆矩阵
        matrix_file = args.output.replace('.txt', '_matrices.txt')
        with open(matrix_file, 'w') as mf:
            mf.write("# Confusion Matrices for Each MAF Bin\n")
            mf.write("# Format: [Predicted_0, Predicted_1, Predicted_2] vs True\n")
            for bin_key in sorted(results.keys()):
                res = results[bin_key]
                matrix = res['confusion_matrix']
                mf.write(f"\n# MAF Bin: {bin_key}\n")
                mf.write(f"# Concordance: {res['concordance']:.6f}\n")
                mf.write("# Confusion Matrix (rows=predicted, cols=true):\n")
                mf.write("#          True_0    True_1    True_2\n")
                for i, row in enumerate(matrix):
                    mf.write(f"Pred_{i}   {row[0]:.2f}    {row[1]:.2f}    {row[2]:.2f}\n")
        
        print(f"\n✓ Confusion matrices saved to: {matrix_file}")
    
    print(f"\n✓ Results saved to: {args.output}")


if __name__ == "__main__":
    main()
