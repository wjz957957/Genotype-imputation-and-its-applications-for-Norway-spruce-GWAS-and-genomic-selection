#!/bin/bash
#SBATCH -p shared
#SBATCH -A naiss2025-5-732
#SBATCH -J compare_scaffold_AF
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.err

set -euo pipefail

WORKDIR="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis"

TARGET="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/Chr11_MAF-347-4K.vcf.gz"

REFERENCE="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/reference-population-phased-700-1.vcf.gz"

TARGET_MAF="$WORKDIR/results/target347_scaffold_MAF.tsv"
REFERENCE_MAF="$WORKDIR/results/reference700_scaffold_MAF.tsv"

mkdir -p "$WORKDIR/results" "$WORKDIR/summary"

module load PDC/26.03
module load bioinfo-tools
ml bcftools/1.20

echo "===== SOFTWARE ====="
bcftools --version | head -n 1

TMPDIR=$(mktemp -d "$WORKDIR/results/.compare.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

echo "===== CHECK INPUT FILES ====="

for F in "$TARGET_MAF" "$REFERENCE_MAF" "$TARGET" "$REFERENCE"
do
    if [[ ! -r "$F" ]]; then
        echo "ERROR: Missing input: $F"
        exit 1
    fi
done

for F in "$TARGET_MAF" "$REFERENCE_MAF"
do
    N=$(awk 'END {print NR-1}' "$F")

    if [[ "$N" -ne 4000 ]]; then
        echo "ERROR: Expected 4000 records in $F"
        exit 1
    fi
done

echo "===== CHECK SAMPLE OVERLAP ====="

bcftools query -l "$TARGET" | sort -u > "$TMPDIR/target.samples"
bcftools query -l "$REFERENCE" | sort -u > "$TMPDIR/reference.samples"

N_TARGET=$(wc -l < "$TMPDIR/target.samples")
N_REFERENCE=$(wc -l < "$TMPDIR/reference.samples")

N_OVERLAP=$(comm -12 \
    "$TMPDIR/target.samples" \
    "$TMPDIR/reference.samples" | wc -l)

{
    printf "Metric\tValue\n"
    printf "Target_samples\t%d\n" "$N_TARGET"
    printf "Reference_samples\t%d\n" "$N_REFERENCE"
    printf "Shared_samples\t%d\n" "$N_OVERLAP"
} > "$WORKDIR/summary/sample_overlap.tsv"

echo "Target samples: $N_TARGET"
echo "Reference samples: $N_REFERENCE"
echo "Shared samples: $N_OVERLAP"

if [[ "$N_TARGET" -ne 347 || "$N_REFERENCE" -ne 700 ]]; then
    echo "ERROR: Unexpected sample numbers"
    exit 1
fi

echo "===== JOIN BOTH DATASETS ====="

JOIN="$TMPDIR/joined.tsv"

awk -F '\t' -v OFS='\t' '
function numeric(x) {
    return x ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
}

function mafbin(x) {
    if (x < 0.01) return "0-0.01"
    if (x < 0.05) return "0.01-0.05"
    if (x < 0.10) return "0.05-0.10"
    if (x <= 0.15) return "0.10-0.15"
    if (x <= 0.20) return ">0.15-0.20"
    return ">0.20-0.50"
}

BEGIN {
    print "CHROM","POS","REF","ALT",
          "Target_AC","Target_AN","Target_AF","Target_MAF",
          "Reference_AC","Reference_AN","Reference_AF","Reference_MAF",
          "Delta_AF_abs","Delta_AF_target_minus_reference",
          "Delta_MAF_abs","Target_MAF_bin","Reference_MAF_bin"
}

NR==FNR {

    if (FNR==1) next

    if (NF!=8 || index($4,",")>0) {
        print "ERROR: Invalid or multiallelic reference record at line " FNR > "/dev/stderr"
        bad=1
        next
    }

    key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4

    if (key in ref_seen) {
        print "ERROR: Duplicate reference variant: " $1 ":" $2 > "/dev/stderr"
        bad=1
        next
    }

    for (i=5;i<=8;i++) {
        if (!numeric($i)) {
            print "ERROR: Invalid reference frequency field" > "/dev/stderr"
            bad=1
        }
    }

    ref_seen[key]=1

    rac[key]=$5
    ran[key]=$6
    raf[key]=$7
    rmaf[key]=$8

    next
}

FNR==1 {next}

{
    if (NF!=8 || index($4,",")>0) {
        print "ERROR: Invalid or multiallelic target record at line " FNR > "/dev/stderr"
        bad=1
        next
    }

    key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4

    if (!(key in ref_seen)) {
        print "ERROR: Target variant absent from reference: " $1 ":" $2 > "/dev/stderr"
        bad=1
        next
    }

    if (key in used) {
        print "ERROR: Duplicate target variant: " $1 ":" $2 > "/dev/stderr"
        bad=1
        next
    }

    for (i=5;i<=8;i++) {
        if (!numeric($i)) {
            print "ERROR: Invalid target frequency field" > "/dev/stderr"
            bad=1
        }
    }

    if (bad) next

    if ($7<0 || $7>1 || $8<0 || $8>0.5 ||
        raf[key]<0 || raf[key]>1 || rmaf[key]<0 || rmaf[key]>0.5) {
        print "ERROR: Frequency outside valid range" > "/dev/stderr"
        bad=1
        next
    }

    used[key]=1
    n++

    delta_af=$7-raf[key]
    delta_maf=$8-rmaf[key]

    if (delta_af<0) abs_af=-delta_af
    else abs_af=delta_af

    if (delta_maf<0) abs_maf=-delta_maf
    else abs_maf=delta_maf

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.8f\t%.8f\t%.8f\t%s\t%s\n", \
        $1,$2,$3,$4,$5,$6,$7,$8, \
        rac[key],ran[key],raf[key],rmaf[key], \
        abs_af,delta_af,abs_maf,mafbin($8),mafbin(rmaf[key])
}

END {

    missing=0

    for (k in ref_seen) {
        if (!(k in used)) missing++
    }

    if (bad || n!=4000 || missing>0) {
        print "ERROR: Dataset matching failed" > "/dev/stderr"
        print "Matched:",n,"Missing reference:",missing > "/dev/stderr"
        exit 1
    }
}
' "$REFERENCE_MAF" "$TARGET_MAF" > "$JOIN"

mv "$JOIN" "$WORKDIR/results/scaffold_frequency_comparison.tsv"

JOIN="$WORKDIR/results/scaffold_frequency_comparison.tsv"

echo "Matched variants:"
awk 'END {print NR-1}' "$JOIN"

echo "===== MAF DISTRIBUTION ====="

awk -F '\t' -v OFS='\t' '

function bin(x) {
    if (x < 0.01) return 1
    if (x < 0.05) return 2
    if (x < 0.10) return 3
    if (x <= 0.15) return 4
    if (x <= 0.20) return 5
    return 6
}

BEGIN {
    label[1]="0-0.01"
    label[2]="0.01-0.05"
    label[3]="0.05-0.10"
    label[4]="0.10-0.15"
    label[5]=">0.15-0.20"
    label[6]=">0.20-0.50"
}

NR==1 {next}

{
    target[bin($8)]++
    reference[bin($12)]++
    n++
}

END {

    print "Population","MAF_bin","N_SNPs","Percentage"

    for (i=1;i<=6;i++) {
        printf "Target\t%s\t%d\t%.2f\n",
            label[i],target[i]+0,100*(target[i]+0)/n
    }

    for (i=1;i<=6;i++) {
        printf "Reference\t%s\t%d\t%.2f\n",
            label[i],reference[i]+0,100*(reference[i]+0)/n
    }
}
' "$JOIN" > "$WORKDIR/summary/MAF_distribution.tsv"

echo "===== HIGH MAF SUMMARY ====="

awk -F '\t' -v OFS='\t' '

NR==1 {next}

{
    n++

    d=$13+0

    if ($8>0.15) {
        target_high++
        target_delta+=d
    }

    if ($12>0.15) {
        reference_high++
        reference_delta+=d
    }

    if ($8>0.15 && $12>0.15) both++
    if ($8>0.15 && $12<=0.15) target_only++
    if ($8<=0.15 && $12>0.15) reference_only++
}

END {

    print "Metric","Value"

    printf "Target_MAF_gt_0.15_count\t%d\n",target_high
    printf "Target_MAF_gt_0.15_percent\t%.2f\n",100*target_high/n

    printf "Reference_MAF_gt_0.15_count\t%d\n",reference_high
    printf "Reference_MAF_gt_0.15_percent\t%.2f\n",100*reference_high/n

    printf "Both_groups_MAF_gt_0.15\t%d\n",both
    printf "Target_only_MAF_gt_0.15\t%d\n",target_only
    printf "Reference_only_MAF_gt_0.15\t%d\n",reference_only

    printf "Mean_abs_Delta_AF_target_MAF_gt_0.15\t"
    if (target_high>0) printf "%.8f\n",target_delta/target_high
    else print "NA"

    printf "Mean_abs_Delta_AF_reference_MAF_gt_0.15\t"
    if (reference_high>0) printf "%.8f\n",reference_delta/reference_high
    else print "NA"
}
' "$JOIN" > "$WORKDIR/summary/high_MAF_summary.tsv"

echo "===== AF MISMATCH SUMMARY ====="

awk -F '\t' '

NR==1 {next}

{
    n++
    d=$13+0

    sum+=d
    if (d>max) max=d

    if (d>=0.05) n05++
    if (d>=0.10) n10++

    if ($8>0.15) {
        nhigh++
        high_sum+=d
    }
}

END {

    print "Metric\tValue"

    printf "N_variants\t%d\n",n
    printf "Mean_abs_Delta_AF\t%.8f\n",sum/n
    printf "Max_abs_Delta_AF\t%.8f\n",max

    printf "N_Delta_AF_ge_0.05\t%d\n",n05
    printf "Percentage_Delta_AF_ge_0.05\t%.2f\n",100*n05/n

    printf "N_Delta_AF_ge_0.10\t%d\n",n10
    printf "Percentage_Delta_AF_ge_0.10\t%.2f\n",100*n10/n

    printf "Mean_abs_Delta_AF_target_MAF_gt_0.15\t"
    if (nhigh>0) printf "%.8f\n",high_sum/nhigh
    else print "NA"
}
' "$JOIN" > "$WORKDIR/summary/AF_mismatch_summary.tsv"

awk -F '\t' 'NR>1 {print $13}' "$JOIN" \
    | sort -g > "$TMPDIR/delta_AF_sorted.txt"

read -r MEDIAN P90 < <(
    awk '
    {
        v[NR]=$1
    }

    END {
        n=NR

        if (n%2)
            median=v[(n+1)/2]
        else
            median=(v[n/2]+v[n/2+1])/2

        rank=int(0.90*n)
        if (rank<0.90*n) rank++

        printf "%.8f %.8f\n",median,v[rank]
    }
    ' "$TMPDIR/delta_AF_sorted.txt"
)

{
    printf "Median_abs_Delta_AF\t%s\n" "$MEDIAN"
    printf "P90_abs_Delta_AF\t%s\n" "$P90"
} >> "$WORKDIR/summary/AF_mismatch_summary.tsv"

echo "===== MAF BIN TRANSITIONS ====="

awk -F '\t' -v OFS='\t' '

BEGIN {
    label[1]="0-0.01"
    label[2]="0.01-0.05"
    label[3]="0.05-0.10"
    label[4]="0.10-0.15"
    label[5]=">0.15-0.20"
    label[6]=">0.20-0.50"
}

NR==1 {next}

{
    key=$16 SUBSEP $17
    count[key]++
}

END {

    print "Target_MAF_bin","Reference_MAF_bin","N_SNPs"

    for (i=1;i<=6;i++) {
        for (j=1;j<=6;j++) {
            key=label[i] SUBSEP label[j]

            if (count[key]>0)
                print label[i],label[j],count[key]
        }
    }
}
' "$JOIN" > "$WORKDIR/summary/MAF_bin_transitions.tsv"

echo "===== GENOTYPE CALL COMPLETENESS ====="

awk -F '\t' '

NR==1 {next}

{
    n++

    if ($6<694) target_incomplete++
    if ($10<1400) reference_incomplete++

    if (n==1 || $6<min_target) min_target=$6
    if (n==1 || $10<min_reference) min_reference=$10
}

END {

    print "Metric\tValue"

    printf "Target_variants_AN_below_694\t%d\n",target_incomplete
    printf "Reference_variants_AN_below_1400\t%d\n",reference_incomplete

    printf "Minimum_Target_AN\t%d\n",min_target
    printf "Minimum_Reference_AN\t%d\n",min_reference
}
' "$JOIN" > "$WORKDIR/summary/genotype_call_completeness.tsv"

echo "===== ANALYSIS FINISHED ====="

ls -lh "$WORKDIR/results/scaffold_frequency_comparison.tsv"
ls -lh "$WORKDIR/summary/"

