#!/usr/bin/env python3
"""
Global R² Calculator with INFO Filtering - Memory Efficient (Incremental Algorithm)
按INFO评分分层过滤，计算不同质量阈值下的全局R²
自动匹配真实VCF和填充VCF的位点，计算全局R²
内存占用恒定 < 50MB，可处理任意大小数据
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
        """返回当前统计信息（用于调试）"""
        return {
            'n': self.n,
            'sum_x': self.sum_x,
            'sum_y': self.sum_y,
            'sum_xy': self.sum_xy,
            'sum_x2': self.sum_x2,
            'sum_y2': self.sum_y2
        }
    
    def reset(self):
        """重置计数器（用于新的INFO阈值）"""
        self.n = 0
        self.sum_x = 0.0
        self.sum_y = 0.0
        self.sum_xy = 0.0
        self.sum_x2 = 0.0
        self.sum_y2 = 0.0


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
    if gt == '0|0':
        return 0.0
    elif gt == '0|1' or gt == '1|0':
        return 1.0
    elif gt == '1|1':
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
                return float(parts[2])
        except (ValueError, IndexError):
            pass
    
    return np.nan


def extract_info_score(info_field):
    """从INFO列提取INFO评分"""
    # INFO列格式: IMP;AC=50;AN=100;AF=0.83941;INFO=0.431171
    try:
        for item in info_field.split(';'):
            if item.startswith('R2='):
                return float(item.split('=')[1])
    except:
        pass
    return np.nan


def calculate_r2_for_threshold(true_vcf, imputed_vcf, info_threshold, verbose=True):
    """
    计算特定INFO阈值下的全局R²
    """
    if verbose:
        print(f"\n  Processing INFO ≥ {info_threshold}...")
    
    # 获取样本信息
    true_samples = get_vcf_samples(true_vcf)
    imputed_samples = get_vcf_samples(imputed_vcf)
    n_samples = len(true_samples)
    
    # 验证样本数量
    if len(true_samples) != len(imputed_samples):
        print(f"  ERROR: Sample count mismatch!")
        return None
    
    # 初始化增量计算器
    corr = IncrementalCorrelation()
    
    # 打开两个VCF文件
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
    filtered_snps = 0
    skipped_snps = 0
    total_data_points = 0
    
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
        
        # 检查INFO评分
        if len(imp_cols) > 7:
            info_score = extract_info_score(imp_cols[7])
            if np.isnan(info_score) or info_score < info_threshold:
                # INFO评分低于阈值，跳过此位点
                # 需要读取两个文件到下一个位点
                # 但需要保持同步
                if true_snp_id < imp_snp_id:
                    true_line = true_f.readline()
                elif true_snp_id > imp_snp_id:
                    imputed_line = imputed_f.readline()
                else:
                    true_line = true_f.readline()
                    imputed_line = imputed_f.readline()
                filtered_snps += 1
                continue
        
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
        else:
            # 位点匹配成功！
            matched_snps += 1
            
            # 提取真实基因型
            true_gt_values = []
            for i in range(n_samples):
                if i < len(true_cols[9:]):
                    gt = true_cols[9+i].split(':')[0]
                    gt_val = convert_true_genotype(gt)
                    true_gt_values.append(gt_val)
                else:
                    true_gt_values.append(np.nan)
            
            # 提取填充剂量
            imp_ds_values = []
            for i in range(n_samples):
                if i < len(imp_cols[9:]):
                    ds_val = extract_dosage(imp_cols[9+i])
                    imp_ds_values.append(ds_val)
                else:
                    imp_ds_values.append(np.nan)
            
            # 增量更新统计量
            for gt, ds in zip(true_gt_values, imp_ds_values):
                if not np.isnan(gt) and not np.isnan(ds):
                    corr.add(ds, gt)
                    total_data_points += 1
            
            # 读取下一行
            true_line = true_f.readline()
            imputed_line = imputed_f.readline()
    
    # 关闭文件
    true_f.close()
    imputed_f.close()
    
    if total_data_points < 3:
        if verbose:
            print(f"    WARNING: Insufficient data points ({total_data_points})")
        return None
    
    r2 = corr.get_r2()
    stats = corr.get_stats()
    
    return {
        'r2': r2,
        'matched_snps': matched_snps,
        'filtered_snps': filtered_snps,
        'skipped_snps': skipped_snps,
        'total_data_points': total_data_points,
        'mean_ds': stats['sum_x'] / stats['n'],
        'mean_gt': stats['sum_y'] / stats['n']
    }


def calculate_global_r2_with_info_filtering(true_vcf, imputed_vcf, info_thresholds, verbose=True):
    """
    计算多个INFO阈值下的全局R²
    """
    
    print("=" * 80)
    print("Global R² Calculator with INFO Filtering (Incremental Algorithm)")
    print("=" * 80)
    print(f"True VCF:     {true_vcf}")
    print(f"Imputed VCF:  {imputed_vcf}")
    print(f"INFO thresholds: {info_thresholds}")
    print("=" * 80)
    
    # 获取基本样本信息
    print("\n[1/5] Checking VCF files...")
    
    true_samples = get_vcf_samples(true_vcf)
    imputed_samples = get_vcf_samples(imputed_vcf)
    
    print(f"  True VCF samples:     {len(true_samples)}")
    print(f"  Imputed VCF samples:  {len(imputed_samples)}")
    
    if len(true_samples) != len(imputed_samples):
        print("  ERROR: Sample count mismatch!")
        sys.exit(1)
    
    n_samples = len(true_samples)
    print(f"  Number of samples: {n_samples}")
    
    print("\n[2/5] Processing INFO thresholds...")
    print("  This will scan the VCF files multiple times (once per threshold)")
    print(f"  Total thresholds: {len(info_thresholds)}")
    
    results = []
    
    for i, threshold in enumerate(info_thresholds, 1):
        print(f"\n  [{i}/{len(info_thresholds)}] INFO ≥ {threshold}")
        
        result = calculate_r2_for_threshold(true_vcf, imputed_vcf, threshold, verbose)
        
        if result:
            results.append({
                'threshold': threshold,
                **result
            })
            print(f"    ✓ R² = {result['r2']:.6f} (based on {result['matched_snps']:,} SNPs, "
                  f"{result['total_data_points']:,} data points)")
        else:
            print(f"    ✗ Failed to calculate R² (insufficient data)")
            results.append({
                'threshold': threshold,
                'r2': np.nan,
                'matched_snps': 0,
                'total_data_points': 0
            })
    
    # 输出汇总结果
    print("\n" + "=" * 80)
    print("SUMMARY RESULTS")
    print("=" * 80)
    print(f"{'INFO ≥':<12} {'R²':<12} {'Matched SNPs':<15} {'Data Points':<15} {'Mean DS':<12} {'Mean GT':<12}")
    print("-" * 80)
    
    for result in results:
        if not np.isnan(result['r2']):
            print(f"{result['threshold']:<12.2f} {result['r2']:<12.6f} "
                  f"{result['matched_snps']:<15,} {result['total_data_points']:<15,} "
                  f"{result.get('mean_ds', 0):<12.4f} {result.get('mean_gt', 0):<12.4f}")
        else:
            print(f"{result['threshold']:<12.2f} {'NA':<12} "
                  f"{result['matched_snps']:<15,} {result['total_data_points']:<15,} "
                  f"{'NA':<12} {'NA':<12}")
    
    print("=" * 80)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Calculate global R² with INFO filtering for imputed VCF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 使用默认INFO阈值（0.0到0.9，步长0.1）
  python global_r2_info.py -t real.vcf.gz -i imputed.vcf.gz -o results.txt
  
  # 自定义INFO阈值
  python global_r2_info.py -t real.vcf.gz -i imputed.vcf.gz --thresholds 0.0 0.2 0.5 0.8
  
  # 详细阈值列表
  python global_r2_info.py -t real.vcf.gz -i imputed.vcf.gz --thresholds 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
  
  # 静默模式
  python global_r2_info.py -t real.vcf.gz -i imputed.vcf.gz -q
        """
    )
    
    parser.add_argument('-t', '--true-vcf', required=True,
                        help='True genotype VCF file')
    parser.add_argument('-i', '--imputed-vcf', required=True,
                        help='Imputed VCF file')
    parser.add_argument('-o', '--output', default='global_R2_by_INFO.txt',
                        help='Output file (default: global_R2_by_INFO.txt)')
    parser.add_argument('--thresholds', nargs='+', type=float,
                        default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                        help='INFO thresholds (default: 0.0 0.1 0.2 ... 0.9)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress detailed progress output')
    
    args = parser.parse_args()
    
    # 排序阈值
    thresholds = sorted(args.thresholds)
    
    # 计算全局R²
    results = calculate_global_r2_with_info_filtering(
        args.true_vcf, 
        args.imputed_vcf, 
        thresholds,
        verbose=not args.quiet
    )
    
    # 保存结果
    with open(args.output, 'w') as f:
        f.write("# Global R² Calculation Results with INFO Filtering\n")
        f.write("# ==================================================\n")
        f.write(f"# True VCF:      {args.true_vcf}\n")
        f.write(f"# Imputed VCF:   {args.imputed_vcf}\n")
        f.write(f"# INFO thresholds: {', '.join([str(t) for t in thresholds])}\n")
        f.write("# ==================================================\n")
        f.write("INFO_threshold\tR²\tMatched_SNPs\tTotal_Data_Points\tMean_DS\tMean_GT\n")
        
        for result in results:
            if not np.isnan(result['r2']):
                f.write(f"{result['threshold']:.2f}\t{result['r2']:.10f}\t"
                       f"{result['matched_snps']}\t{result['total_data_points']}\t"
                       f"{result.get('mean_ds', 0):.6f}\t{result.get('mean_gt', 0):.6f}\n")
            else:
                f.write(f"{result['threshold']:.2f}\tNA\t{result['matched_snps']}\t{result['total_data_points']}\tNA\tNA\n")
    
    print(f"\n✓ Results saved to: {args.output}")
    
    # 额外生成一个R²随INFO阈值变化的汇总表
    summary_file = args.output.replace('.txt', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("# R² by INFO threshold\n")
        f.write("# INFO_threshold\tR²\n")
        for result in results:
            if not np.isnan(result['r2']):
                f.write(f"{result['threshold']:.2f}\t{result['r2']:.10f}\n")
    
    print(f"✓ Summary saved to: {summary_file}")


if __name__ == "__main__":
    main()
