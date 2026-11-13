#!/usr/bin/env python3
"""
MAF-stratified Concordance Calculator - Memory Efficient (Incremental Algorithm)
自动匹配真实VCF和填充VCF的位点，按MAF区间计算基因型一致性
只计算填充VCF中带有IMP标记的位点
"""

import gzip
import numpy as np
import argparse
import sys
from collections import defaultdict


class IncrementalConcordance:
    """增量计算Concordance的类"""
    
    def __init__(self):
        self.matrix = np.zeros((3, 3), dtype=np.float64)
        self.total_samples = 0
    
    def add_sample(self, true_gt, gp_probs):
        if true_gt == 'NA' or true_gt is None:
            return
        
        true_idx = int(true_gt)
        
        for pred_idx in range(3):
            if not np.isnan(gp_probs[pred_idx]):
                self.matrix[pred_idx, true_idx] += gp_probs[pred_idx]
        
        self.total_samples += 1
    
    def get_concordance(self):
        if self.total_samples == 0:
            return np.nan
        
        diagonal_sum = np.sum(np.diag(self.matrix))
        total = np.sum(self.matrix)
        
        if total > 0:
            return diagonal_sum / total
        else:
            return np.nan
    
    def get_matrix(self):
        return self.matrix
    
    def get_stats(self):
        concordance = self.get_concordance()
        return {
            'concordance': concordance,
            'total_weighted_samples': self.total_samples,
            'matrix': self.matrix
        }


def get_vcf_samples(vcf_path):
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                return line.strip().split('\t')[9:]
    return []


def convert_true_genotype(gt):
    if gt == '0/0':
        return 0
    elif gt == '0/1' or gt == '1/0':
        return 1
    elif gt == '1/1':
        return 2
    else:
        return 'NA'


def extract_gp(sample_data):
    parts = sample_data.split(':')
    
    gp = [np.nan, np.nan, np.nan]
    
    if 'GP' in parts:
        gp_idx = parts.index('GP')
        if gp_idx + 1 < len(parts):
            try:
                gp = [float(x) for x in parts[gp_idx+1].split(',')]
            except:
                pass
    elif len(parts) >= 3:
        try:
            gp = [float(x) for x in parts[2].split(',')]
        except:
            pass
    elif len(parts) >= 2:
        try:
            gp = [float(x) for x in parts[1].split(',')]
        except:
            pass
    
    return gp


def extract_maf_from_info(info_str):
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


# <-- NEW: 检查IMP标记的函数
def has_imp_tag(info_str):
    """
    检查INFO字段是否包含IMP标记（表示该位点是经过填充的）
    """
    if not info_str:
        return False
    
    fields = info_str.split(';')
    return 'IMP' in fields


class MAFStratifiedConcordanceCalculator:
    """按MAF区间分层计算Concordance的类"""
    
    def __init__(self, maf_bins):
        self.maf_bins = maf_bins
        self.concordance_calcs = {f"{low:.3f}-{high:.3f}": IncrementalConcordance() 
                                   for low, high in maf_bins}
        self.snp_counts = {f"{low:.3f}-{high:.3f}": 0 for low, high in maf_bins}
        self.total_snps_processed = 0
    
    def add_snp(self, maf, gp_values, gt_values):
        bin_key = None
        for low, high in self.maf_bins:
            if low <= maf < high or (high == 1.0 and maf >= low):
                bin_key = f"{low:.3f}-{high:.3f}"
                break
        
        if bin_key is None:
            return 0
        
        conc_calc = self.concordance_calcs[bin_key]
        valid_samples = 0
        
        for gp, gt in zip(gp_values, gt_values):
            if gt != 'NA' and not np.isnan(gp[0]):
                conc_calc.add_sample(gt, gp)
                valid_samples += 1
        
        if valid_samples > 0:
            self.snp_counts[bin_key] += 1
            self.total_snps_processed += 1
        
        return valid_samples
    
    def get_results(self):
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
    """计算按MAF分层的Concordance - 只计算带有IMP标记的填充位点"""
    
    print("=" * 70)
    print("MAF-Stratified Concordance Calculator (Incremental Algorithm) - IMP-only mode")
    print("=" * 70)
    print(f"True VCF:     {true_vcf}")
    print(f"Imputed VCF:  {imputed_vcf}")
    print(f"MAF bins:     {len(maf_bins)} intervals")
    print("NOTE: Only SNPs with IMP tag in imputed VCF will be included")
    if region:
        print(f"Region:       {region}")
    print("=" * 70)
    
    true_samples = get_vcf_samples(true_vcf)
    imputed_samples = get_vcf_samples(imputed_vcf)
    
    print("\n[1/5] Checking VCF files...")
    print(f"  True VCF samples:     {len(true_samples)}")
    print(f"  Imputed VCF samples:  {len(imputed_samples)}")
    
    if len(true_samples) != len(imputed_samples):
        print(f"  ERROR: Sample count mismatch!")
        sys.exit(1)
    
    n_samples = len(true_samples)
    print(f"  Number of samples:    {n_samples}")
    
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
    
    print("\n[2/5] Initializing MAF-stratified calculator...")
    calculator = MAFStratifiedConcordanceCalculator(maf_bins)
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
    filtered_snps = 0
    total_samples_processed = 0
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
        
        if region_chrom:
            if true_chrom != region_chrom:
                true_line = true_f.readline()
                filtered_snps += 1
                continue
            if region_start and (int(true_pos) < region_start or int(true_pos) > region_end):
                true_line = true_f.readline()
                filtered_snps += 1
                continue
        
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
                true_gt_values.append('NA')
        
        imp_gp_values = []
        for i in range(n_samples):
            if i < len(imp_cols[9:]):
                gp = extract_gp(imp_cols[9+i])
                imp_gp_values.append(gp)
            else:
                imp_gp_values.append([np.nan, np.nan, np.nan])
        
        valid_samples = calculator.add_snp(maf, imp_gp_values, true_gt_values)
        total_samples_processed += valid_samples
        
        true_line = true_f.readline()
        imputed_line = imputed_f.readline()
        
        if verbose and matched_snps % 10000 == 0:
            print(f"    Processed {matched_snps} matched SNPs, "
                  f"{total_samples_processed} samples, "
                  f"skipped {skipped_snps}, filtered {filtered_snps}, "
                  f"no IMP: {no_imp_snps}")  # <-- NEW
    
    true_f.close()
    imputed_f.close()
    
    print(f"\n  ✓ Processing complete!")
    print(f"    Total SNPs scanned:      {total_snps}")
    print(f"    Matched SNPs:            {matched_snps}")
    print(f"    Unmatched SNPs:          {skipped_snps}")
    print(f"    Filtered by region:      {filtered_snps}")
    print(f"    SNPs without IMP tag:    {no_imp_snps}")  # <-- NEW
    print(f"    SNPs without MAF info:   {no_maf_snps}")
    print(f"    Total weighted samples:  {total_samples_processed}")
    
    print("\n[4/5] Calculating concordance for each MAF bin...")
    results = calculator.get_results()
    
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
    return [
        (0.00, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04),
        (0.04, 0.05), (0.05, 0.06), (0.06, 0.07), (0.07, 0.08),
        (0.08, 0.09), (0.09, 0.10), (0.10, 0.15), (0.15, 0.20),
        (0.20, 0.25), (0.25, 0.30), (0.30, 0.35), (0.35, 0.40),
        (0.40, 0.45), (0.45, 0.51)
    ]


def main():
    parser = argparse.ArgumentParser(
        description='Calculate MAF-stratified genotype concordance between true and imputed VCF files (IMP-only mode)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""..."""
    )
    
    parser.add_argument('-t', '--true-vcf', required=True,
                        help='True genotype VCF file')
    parser.add_argument('-i', '--imputed-vcf', required=True,
                        help='Imputed VCF file (must contain GP and IMP tag)')
    parser.add_argument('-o', '--output', default='maf_stratified_concordance.txt',
                        help='Output file')
    parser.add_argument('--maf-bins', nargs='+', action='append',
                        help='MAF bins in format "min,max"')
    parser.add_argument('-r', '--region', default=None,
                        help='Region to analyze (e.g., "chr1" or "chr1:1000-2000")')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress detailed output')
    
    args = parser.parse_args()
    
    # 处理MAF bins
    if args.maf_bins:
        maf_bins = []
        for bin_spec in args.maf_bins[0]:
            parts = bin_spec.split(',')
            if len(parts) == 2:
                maf_bins.append((float(parts[0]), float(parts[1])))
    else:
        maf_bins = get_default_maf_bins()
    
    # 计算
    results = calculate_maf_stratified_concordance(
        args.true_vcf,
        args.imputed_vcf,
        maf_bins,
        region=args.region,
        verbose=not args.quiet
    )
    
    # 保存结果
    with open(args.output, 'w') as f:
        f.write("# MAF-Stratified Concordance Results (IMP-only mode)\n")
        f.write("# ==================================================\n")
        f.write(f"# True VCF:    {args.true_vcf}\n")
        f.write(f"# Imputed VCF: {args.imputed_vcf}\n")
        f.write("# NOTE: Only SNPs with IMP tag were included\n")
        f.write("# ==================================================\n")
        f.write("MAF_Bin\tSNPs\tWeighted_Samples\tConcordance\n")
        
        for bin_key in sorted(results.keys()):
            res = results[bin_key]
            f.write(f"{bin_key}\t{res['n_snps']}\t{res['n_weighted_samples']}\t")
            if not np.isnan(res['concordance']):
                f.write(f"{res['concordance']:.10f}\n")
            else:
                f.write("NA\n")
    
    print(f"\n✓ Results saved to: {args.output}")

if __name__ == "__main__":
    main()