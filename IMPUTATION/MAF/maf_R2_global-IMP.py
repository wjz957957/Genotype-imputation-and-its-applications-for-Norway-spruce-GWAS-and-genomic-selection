#!/usr/bin/env python3
"""
MAF-stratified R² Calculator - Memory Efficient (Incremental Algorithm)
自动匹配真实VCF和填充VCF的位点，按MAF区间计算全局R²
只计算填充VCF中带有IMP标记的位点
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
        
        numerator = self.n * self.sum_xy - self.sum_x * self.sum_y
        
        denominator_x = self.n * self.sum_x2 - self.sum_x * self.sum_x
        denominator_y = self.n * self.sum_y2 - self.sum_y * self.sum_y
        
        if denominator_x <= 0 or denominator_y <= 0:
            return np.nan
        
        denominator = sqrt(denominator_x * denominator_y)
        
        if denominator == 0:
            return np.nan
        
        r = numerator / denominator
        return r * r
    
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
    
    if len(parts) >= 2:
        try:
            if 'DS' in parts:
                ds_idx = parts.index('DS') + 1
                if ds_idx < len(parts):
                    return float(parts[ds_idx])
            else:
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
            except (ValueError, IndexError):
                continue
        elif field.startswith('AF='):
            try:
                af = float(field.split('=')[1])
                maf = min(af, 1 - af)
                if 0 <= maf <= 0.5:
                    return maf
            except (ValueError, IndexError):
                continue
    
    return None


# <-- NEW: 检查IMP标记的函数
def has_imp_tag(info_str):
    """
    检查INFO字段是否包含IMP标记（表示该位点是经过填充的）
    
    Parameters:
    -----------
    info_str : str
        VCF文件中INFO字段的内容
    
    Returns:
    --------
    bool
        True表示该位点带有IMP标记（是填充位点），False表示没有IMP标记
    """
    if not info_str:
        return False
    
    fields = info_str.split(';')
    return 'IMP' in fields


class MAFStratifiedCalculator:
    """按MAF区间分层计算R²的类"""
    
    def __init__(self, maf_bins):
        self.maf_bins = maf_bins
        self.correlations = {f"{low:.3f}-{high:.3f}": IncrementalCorrelation() 
                            for low, high in maf_bins}
        self.snp_counts = {f"{low:.3f}-{high:.3f}": 0 for low, high in maf_bins}
        self.total_snps_processed = 0
    
    def add_snp(self, maf, ds_values, gt_values):
        bin_key = None
        for low, high in self.maf_bins:
            if low <= maf < high or (high == 1.0 and maf >= low):
                bin_key = f"{low:.3f}-{high:.3f}"
                break
        
        if bin_key is None:
            return 0
        
        corr = self.correlations[bin_key]
        valid_points = 0
        
        for ds, gt in zip(ds_values, gt_values):
            if not np.isnan(ds) and not np.isnan(gt):
                corr.add(ds, gt)
                valid_points += 1
        
        if valid_points > 0:
            self.snp_counts[bin_key] += 1
            self.total_snps_processed += 1
        
        return valid_points
    
    def get_results(self):
        results = {}
        for bin_key in self.correlations:
            r2 = self.correlations[bin_key].get_r2()
            stats = self.correlations[bin_key].get_stats()
            results[bin_key] = {
                'r2': r2,
                'n_snps': self.snp_counts[bin_key],
                'n_data_points': stats['n'],
                'mean_ds': stats['mean_x'],
                'mean_gt': stats['mean_y']
            }
        return results


def calculate_maf_stratified_r2(true_vcf, imputed_vcf, maf_bins, verbose=True):
    """计算按MAF分层的全局R² - 只计算带有IMP标记的填充位点"""
    
    print("=" * 70)
    print("MAF-Stratified R² Calculator (Incremental Algorithm) - IMP-only mode")
    print("=" * 70)
    print(f"True VCF:     {true_vcf}")
    print(f"Imputed VCF:  {imputed_vcf}")
    print(f"MAF bins:     {len(maf_bins)} intervals")
    print("NOTE: Only SNPs with IMP tag in imputed VCF will be included")
    print("=" * 70)
    
    true_samples = get_vcf_samples(true_vcf)
    imputed_samples = get_vcf_samples(imputed_vcf)
    
    print(f"\n[1/5] Checking VCF files...")
    print(f"  True VCF samples:     {len(true_samples)}")
    print(f"  Imputed VCF samples:  {len(imputed_samples)}")
    
    if len(true_samples) != len(imputed_samples):
        print(f"  ERROR: Sample count mismatch! {len(true_samples)} vs {len(imputed_samples)}")
        sys.exit(1)
    
    n_samples = len(true_samples)
    print(f"  Number of samples to analyze: {n_samples}")
    
    print("\n[2/5] Initializing MAF-stratified calculator...")
    calculator = MAFStratifiedCalculator(maf_bins)
    print(f"  Created {len(maf_bins)} MAF bins")
    
    print("\n[3/5] Streaming VCF files and matching SNPs...")
    
    true_f = gzip.open(true_vcf, 'rt')
    imputed_f = gzip.open(imputed_vcf, 'rt')
    
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
    
    total_snps = 0
    matched_snps = 0
    skipped_snps = 0
    total_data_points = 0
    no_maf_snps = 0
    no_imp_snps = 0  # <-- NEW: 统计没有IMP标记的SNP数量
    
    while true_line and imputed_line:
        total_snps += 1
        
        true_cols = true_line.strip().split('\t')
        true_chrom = true_cols[0].replace('chr', '')
        true_pos = true_cols[1]
        true_snp_id = f"{true_chrom}:{true_pos}"
        
        imp_cols = imputed_line.strip().split('\t')
        imp_chrom = imp_cols[0].replace('chr', '')
        imp_pos = imp_cols[1]
        imp_snp_id = f"{imp_chrom}:{imp_pos}"
        
        if true_snp_id < imp_snp_id:
            true_line = true_f.readline()
            skipped_snps += 1
            continue
        elif true_snp_id > imp_snp_id:
            imputed_line = imputed_f.readline()
            skipped_snps += 1
            continue
        
        matched_snps += 1
        
        # <-- NEW: 检查填充VCF的INFO字段是否包含IMP标记
        imp_info = imp_cols[7] if len(imp_cols) > 7 else ""
        if not has_imp_tag(imp_info):
            # 该位点没有IMP标记，跳过（不属于填充位点）
            no_imp_snps += 1
            true_line = true_f.readline()
            imputed_line = imputed_f.readline()
            continue
        
        info_fields = true_cols[7] if len(true_cols) > 7 else ""
        maf = extract_maf_from_info(info_fields)
        
        if maf is None or np.isnan(maf):
            no_maf_snps += 1
            true_line = true_f.readline()
            imputed_line = imputed_f.readline()
            continue
        
        true_gt_values = []
        for i in range(n_samples):
            if i < len(true_cols[9:]):
                gt = true_cols[9+i].split(':')[0]
                gt_val = convert_true_genotype(gt)
                true_gt_values.append(gt_val)
            else:
                true_gt_values.append(np.nan)
        
        imp_ds_values = []
        for i in range(n_samples):
            if i < len(imp_cols[9:]):
                ds_val = extract_dosage(imp_cols[9+i])
                imp_ds_values.append(ds_val)
            else:
                imp_ds_values.append(np.nan)
        
        valid_points = calculator.add_snp(maf, imp_ds_values, true_gt_values)
        total_data_points += valid_points
        
        true_line = true_f.readline()
        imputed_line = imputed_f.readline()
        
        if verbose and matched_snps % 10000 == 0:
            print(f"    Processed {matched_snps} matched SNPs, "
                  f"{total_data_points} data points, "
                  f"skipped {skipped_snps}, "
                  f"no IMP: {no_imp_snps}, "  # <-- NEW
                  f"no MAF: {no_maf_snps}")
    
    true_f.close()
    imputed_f.close()
    
    print(f"\n  ✓ Processing complete!")
    print(f"    Total SNPs processed:    {total_snps}")
    print(f"    Matched SNPs:            {matched_snps}")
    print(f"    Unmatched SNPs:          {skipped_snps}")
    print(f"    SNPs without IMP tag:    {no_imp_snps}")  # <-- NEW
    print(f"    SNPs without MAF info:   {no_maf_snps}")
    print(f"    Total data points:       {total_data_points}")
    
    print("\n[4/5] Calculating R² for each MAF bin...")
    results = calculator.get_results()
    
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
    return [
        (0.00, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04),
        (0.04, 0.05), (0.05, 0.06), (0.06, 0.07), (0.07, 0.08),
        (0.08, 0.09), (0.09, 0.10), (0.10, 0.15), (0.15, 0.20),
        (0.20, 0.25), (0.25, 0.30), (0.30, 0.35), (0.35, 0.40),
        (0.40, 0.45), (0.45, 0.51)
    ]


def main():
    parser = argparse.ArgumentParser(
        description='Calculate MAF-stratified R² between true and imputed VCF files (IMP-only mode)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python maf_R2_global.py -t real.vcf.gz -i imputed.vcf.gz -o results.txt
  
Note:
  - Only SNPs with IMP tag in the imputed VCF will be included in calculation
  - True VCF must contain MAF or AF in INFO field
        """
    )
    
    parser.add_argument('-t', '--true-vcf', required=True,
                        help='True genotype VCF file (must contain MAF or AF in INFO field)')
    parser.add_argument('-i', '--imputed-vcf', required=True,
                        help='Imputed VCF file (must contain DS dosage values and IMP tag)')
    parser.add_argument('-o', '--output', default='maf_stratified_R2.txt',
                        help='Output file (default: maf_stratified_R2.txt)')
    parser.add_argument('--maf-bins', nargs='+', action='append',
                        help='MAF bins in format "min,max" (can specify multiple)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress detailed progress output')
    
    args = parser.parse_args()
    
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
    
    results = calculate_maf_stratified_r2(
        args.true_vcf, 
        args.imputed_vcf, 
        maf_bins,
        verbose=not args.quiet
    )
    
    if results is None:
        sys.exit(1)
    
    with open(args.output, 'w') as f:
        f.write("# MAF-Stratified R² Calculation Results (IMP-only mode)\n")
        f.write("# ======================================================\n")
        f.write(f"# True VCF:      {args.true_vcf}\n")
        f.write(f"# Imputed VCF:   {args.imputed_vcf}\n")
        f.write("# NOTE: Only SNPs with IMP tag were included\n")
        f.write("# ======================================================\n")
        f.write("MAF_Bin\tSNPs_Count\tData_Points\tMean_DS\tMean_GT\tR2\n")
        
        for bin_key in sorted(results.keys()):
            res = results[bin_key]
            f.write(f"{bin_key}\t{res['n_snps']}\t{res['n_data_points']}\t"
                   f"{res['mean_ds']:.6f}\t{res['mean_gt']:.6f}\t")
            if not np.isnan(res['r2']):
                f.write(f"{res['r2']:.10f}\n")
            else:
                f.write("NA\n")
        
        f.write("# ======================================================\n")
    
    print(f"\n✓ Results saved to: {args.output}")


if __name__ == "__main__":
    main()