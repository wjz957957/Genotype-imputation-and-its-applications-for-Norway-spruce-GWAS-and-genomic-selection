#!/usr/bin/env python3
"""
统计VCF文件中MAF的分布
"""

import gzip
import numpy as np
import argparse
import sys

def extract_maf_from_info(info_str):
    """从INFO字段提取MAF"""
    if not info_str or info_str == '.':
        return None
    
    fields = info_str.split(';')
    for field in fields:
        if field.startswith('MAF='):
            try:
                return float(field.split('=')[1])
            except:
                return None
    return None

def analyze_maf_distribution(vcf_path, output_file=None):
    """
    分析VCF中MAF的分布
    """
    print("=" * 70)
    print("MAF Distribution Analysis")
    print("=" * 70)
    print(f"VCF file: {vcf_path}")
    print("=" * 70)
    
    # 定义MAF区间（左闭右开）
    maf_bins = [
        (0.00, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04), (0.04, 0.05),
        (0.05, 0.06), (0.06, 0.07), (0.07, 0.08), (0.08, 0.09), (0.09, 0.10),
        (0.10, 0.15), (0.15, 0.20), (0.20, 0.25), (0.25, 0.30),
        (0.30, 0.35), (0.35, 0.40), (0.40, 0.45), (0.45, 0.51)
    ]
    
    # 初始化统计
    bin_counts = {f"{low:.3f}-{high:.3f}": 0 for low, high in maf_bins}
    all_mafs = []
    missing_maf = 0
    total_snps = 0
    
    print("\nReading VCF file...")
    
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            total_snps += 1
            cols = line.strip().split('\t')
            info = cols[7] if len(cols) > 7 else ""
            
            maf = extract_maf_from_info(info)
            
            if maf is None:
                missing_maf += 1
                continue
            
            all_mafs.append(maf)
            
            # 分配到区间
            for low, high in maf_bins:
                if low <= maf < high or (high == 1.0 and maf >= low):
                    bin_key = f"{low:.3f}-{high:.3f}"
                    bin_counts[bin_key] += 1
                    break
    
    # 计算统计量
    all_mafs = np.array(all_mafs)
    
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"Total SNPs analyzed:     {total_snps:,}")
    print(f"SNPs with MAF info:      {len(all_mafs):,}")
    print(f"SNPs without MAF:        {missing_maf:,}")
    print()
    
    if len(all_mafs) > 0:
        print(f"MAF Statistics:")
        print(f"  Minimum MAF:           {np.min(all_mafs):.6f}")
        print(f"  Maximum MAF:           {np.max(all_mafs):.6f}")
        print(f"  Mean MAF:              {np.mean(all_mafs):.6f}")
        print(f"  Median MAF:            {np.median(all_mafs):.6f}")
        print(f"  Standard deviation:    {np.std(all_mafs):.6f}")
        print(f"  25th percentile:       {np.percentile(all_mafs, 25):.6f}")
        print(f"  75th percentile:       {np.percentile(all_mafs, 75):.6f}")
        print()
    
    print("MAF Distribution by Bin:")
    print("-" * 70)
    print(f"{'MAF Bin':<20} {'Count':<12} {'Percentage':<12} {'Cumulative':<12}")
    print("-" * 70)
    
    cumulative = 0
    for low, high in maf_bins:
        bin_key = f"{low:.3f}-{high:.3f}"
        count = bin_counts[bin_key]
        percentage = (count / len(all_mafs) * 100) if len(all_mafs) > 0 else 0
        cumulative += percentage
        bar = '█' * int(percentage / 2) if percentage > 0 else ''
        print(f"{bin_key:<20} {count:<12,} {percentage:<11.2f}% {cumulative:<11.2f}% {bar}")
    
    print("-" * 70)
    print(f"{'Total':<20} {len(all_mafs):<12,} {'100.00':<11}%")
    
    # 特殊检查：是否有MAF < 0.05的位点
    low_maf_count = sum(count for (low, high), count in bin_counts.items() if high <= 0.05)
    if low_maf_count == 0:
        print("\n⚠️  WARNING: No SNPs found with MAF < 0.05")
        print("   This VCF has been filtered for MAF >= 0.05")
    elif low_maf_count > 0:
        print(f"\n✓ Found {low_maf_count:,} SNPs with MAF < 0.05 ({low_maf_count/len(all_mafs)*100:.2f}%)")
    
    # 保存结果到文件
    if output_file:
        with open(output_file, 'w') as f:
            f.write("# MAF Distribution Report\n")
            f.write(f"# VCF file: {vcf_path}\n")
            f.write(f"# Total SNPs: {total_snps}\n")
            f.write(f"# SNPs with MAF: {len(all_mafs)}\n\n")
            
            if len(all_mafs) > 0:
                f.write("MAF Statistics:\n")
                f.write(f"  Min: {np.min(all_mafs):.6f}\n")
                f.write(f"  Max: {np.max(all_mafs):.6f}\n")
                f.write(f"  Mean: {np.mean(all_mafs):.6f}\n")
                f.write(f"  Median: {np.median(all_mafs):.6f}\n")
                f.write(f"  Std: {np.std(all_mafs):.6f}\n\n")
            
            f.write("MAF_Bin\tCount\tPercentage\tCumulative\n")
            cumulative_pct = 0
            for low, high in maf_bins:
                bin_key = f"{low:.3f}-{high:.3f}"
                count = bin_counts[bin_key]
                pct = (count / len(all_mafs) * 100) if len(all_mafs) > 0 else 0
                cumulative_pct += pct
                f.write(f"{bin_key}\t{count}\t{pct:.2f}%\t{cumulative_pct:.2f}%\n")
        
        print(f"\n✓ Detailed report saved to: {output_file}")
    
    return all_mafs, bin_counts

def main():
    parser = argparse.ArgumentParser(description='Analyze MAF distribution in VCF file')
    parser.add_argument('-i', '--input', required=True, help='Input VCF.gz file with MAF info')
    parser.add_argument('-o', '--output', help='Output file for distribution report')
    
    args = parser.parse_args()
    
    analyze_maf_distribution(args.input, args.output)

if __name__ == "__main__":
    main()
