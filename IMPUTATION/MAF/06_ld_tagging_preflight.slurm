#!/bin/bash -l

#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -J chr11_ld_preflight
#SBATCH -t 02:00:00
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH --output=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.err

set -euo pipefail

BASE=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5

WORKDIR="$BASE/05_scaffold_frequency_analysis"
OUT="$WORKDIR/ld_tagging/preflight"

REF="$BASE/reference-population-phased-700-1.vcf.gz"

TARGET="$BASE/Chr11_MAF-347.vcf.gz"
TARGET_MAF="$BASE/Chr11_MAF-347_with_MAF.vcf.gz"

TARGET_4K="$BASE/Chr11_MAF-347-4K-maf-ac-an.vcf.gz"

SCAFFOLD="$BASE/chrom_pos-4k-1.txt"


IMPUTED="$BASE/impute5-output-347-maf-ac-an.vcf.gz"

mkdir -p "$OUT"

exec > >(tee "$OUT/preflight_report.txt") 2>&1

module load bioinfo-tools
module load bcftools/1.20

echo "===== SOFTWARE ====="
bcftools --version | sed -n '1p'

echo
echo "===== INPUT FILE INVENTORY ====="

for F in \
    "$REF" \
    "$TARGET" \
    "$TARGET_MAF" \
    "$TARGET_4K" \
    "$IMPUTED"
do

    echo
    echo "FILE: $F"

    if [[ ! -s "$F" ]]; then
        echo "ERROR: File missing or empty"
        continue
    fi

    echo "Size:"
    ls -lh "$F"

    echo -n "Samples: "
    bcftools query -l "$F" | wc -l

    if [[ -f "${F}.csi" || -f "${F}.tbi" ]]; then
        echo -n "Indexed variants: "
        bcftools index -n "$F"
    else
        echo "Variant count: not indexed"
    fi

done

echo
echo "===== STEP 1: EXTRACT TARGET VARIANT KEYS ====="

TARGET_KEYS="$OUT/target_all_biallelic_snp_keys.tsv"

bcftools view -m2 -M2 -v snps -Ou "$TARGET" \
    | bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    > "$TARGET_KEYS"

echo -n "Target biallelic SNPs: "
wc -l < "$TARGET_KEYS"

echo
echo "===== STEP 2: IDENTIFY VALIDATION SNPs ====="

VALIDATION_KEYS="$OUT/validation_snp_keys.tsv"
SCAFFOLD_KEYS="$OUT/scaffold_snp_keys.tsv"

awk '
NR==FNR {
    scaffold[$1 SUBSEP $2]=1
    next
}

{
    key=$1 SUBSEP $2

    if (key in scaffold)
        print $0 > scaffold_file
    else
        print $0 > validation_file
}
' scaffold_file="$SCAFFOLD_KEYS" \
  validation_file="$VALIDATION_KEYS" \
  "$SCAFFOLD" "$TARGET_KEYS"

echo -n "Scaffold SNPs present in target: "
wc -l < "$SCAFFOLD_KEYS"

echo -n "Validation SNPs: "
wc -l < "$VALIDATION_KEYS"

echo
echo "===== STEP 3: CHECK MISSING SCAFFOLD POSITIONS ====="

awk '
NR==FNR {
    target[$1 SUBSEP $2]=1
    next
}

!(($1 SUBSEP $2) in target) {
    print $0
}
' "$TARGET_KEYS" "$SCAFFOLD" \
    > "$OUT/scaffold_positions_missing_in_target.txt"

echo -n "Missing scaffold positions: "
wc -l < "$OUT/scaffold_positions_missing_in_target.txt"

echo
echo "===== STEP 4: VALIDATION SNP MAF DISTRIBUTION ====="

if bcftools view -h "$TARGET_MAF" \
    | grep -q '^##INFO=<ID=MAF,'; then

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\n' \
        "$TARGET_MAF" \
        > "$OUT/target_all_MAF.tsv"

    awk '
    NR==FNR {
        scaffold[$1 SUBSEP $2]=1
        next
    }

    !(($1 SUBSEP $2) in scaffold) {
        print $0
    }
    ' "$SCAFFOLD" "$OUT/target_all_MAF.tsv" \
        > "$OUT/validation_MAF.tsv"

    awk -F '\t' '
    BEGIN {
        OFS="\t"
        print "MAF_bin","N_SNPs","Percentage"
    }

    $5!="." && $5!="" {

        m=$5+0

        if (m<0.01) b="0-0.01"
        else if (m<0.05) b="0.01-0.05"
        else if (m<0.10) b="0.05-0.10"
        else if (m<0.15) b="0.10-0.15"
        else if (m<0.20) b="0.15-0.20"
        else b="0.20-0.50"

        n[b]++
        total++

        if (m>0.15)
            high++
    }

    END {

        bins[1]="0-0.01"
        bins[2]="0.01-0.05"
        bins[3]="0.05-0.10"
        bins[4]="0.10-0.15"
        bins[5]="0.15-0.20"
        bins[6]="0.20-0.50"

        for (i=1;i<=6;i++) {
            b=bins[i]
            print b,n[b]+0,100*(n[b]+0)/total
        }

        print "TOTAL",total,100
        print "MAF_GT_0.15",high+0,100*(high+0)/total
    }
    ' "$OUT/validation_MAF.tsv" \
        > "$OUT/validation_MAF_distribution.tsv"

    cat "$OUT/validation_MAF_distribution.tsv"

else

    echo "WARNING: INFO/MAF not found in TARGET_MAF"

fi

echo
echo "===== STEP 5: CHECK 22 INCOMPLETE SCAFFOLD SNPs ====="

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/MAF\n' \
    "$TARGET_4K" \
    | awk -F '\t' '
    BEGIN {
        OFS="\t"
        print "CHROM","POS","REF","ALT","AC","AN","MAF","Missing_alleles"
    }

    $6<694 {
        print $1,$2,$3,$4,$5,$6,$7,694-$6
        n++
    }

    END {
        print "Total incomplete scaffold variants:",n+0 > "/dev/stderr"
    }
    ' > "$OUT/incomplete_scaffold_variants.tsv"

cat "$OUT/incomplete_scaffold_variants.tsv"

echo
echo "===== STEP 6: EXTRACT REFERENCE GENOTYPES AT TARGET SITES ====="

TARGET_POS="$OUT/target_all_positions.tsv"

cut -f1,2 "$TARGET_KEYS" \
    | awk '!seen[$1 SUBSEP $2]++' \
    > "$TARGET_POS"

REF_SUBSET="$OUT/reference_at_target_sites.vcf.gz"

bcftools view \
    -T "$TARGET_POS" \
    -m2 -M2 -v snps \
    -Oz \
    -o "$REF_SUBSET" \
    "$REF"

bcftools index -f -c "$REF_SUBSET"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "$REF_SUBSET" \
    > "$OUT/reference_target_site_keys.tsv"

echo -n "Reference variants extracted: "
wc -l < "$OUT/reference_target_site_keys.tsv"

echo
echo "===== STEP 7: EXACT ALLELE MATCHING ====="

awk '
BEGIN {
    FS=OFS="\t"
}

FILENAME==ARGV[1] {
    ref[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4]=1
    next
}

{
    total++

    key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4

    if (key in ref)
        matched++
    else
        missing++
}

END {
    print "Metric","Value"
    print "Target_biallelic_SNPs",total
    print "Exact_reference_matches",matched+0
    print "Missing_or_allele_mismatched",missing+0
    print "Match_percentage",100*(matched+0)/total
}
' "$OUT/reference_target_site_keys.tsv" "$TARGET_KEYS" \
    > "$OUT/reference_target_site_match.tsv"

cat "$OUT/reference_target_site_match.tsv"

echo
echo "===== PREFLIGHT COMPLETED ====="

echo "Output directory: $OUT"

