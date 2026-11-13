#!/bin/bash
#SBATCH -p shared
#SBATCH -A naiss2025-5-732
#SBATCH -J reference700_MAF
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.err

set -euo pipefail

WORKDIR="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis"

INPUT="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/reference-population-phased-700-1.vcf.gz"

SITES="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/chrom_pos-4k-1.txt"

OUT="$WORKDIR/results/reference700_scaffold_MAF.tsv"

mkdir -p "$WORKDIR/results"

echo "===== LOAD SOFTWARE ====="

module load PDC/26.03
module load bioinfo-tools
ml bcftools/1.20

bcftools --version | head -n 1

echo "===== CHECK INPUT ====="

if [[ ! -r "$INPUT" || ! -r "$SITES" ]]; then
    echo "ERROR: Input file missing"
    exit 1
fi

N_SAMPLES=$(bcftools query -l "$INPUT" | awk 'END{print NR+0}')

echo "Reference samples: $N_SAMPLES"

if [[ "$N_SAMPLES" -ne 700 ]]; then
    echo "ERROR: Expected 700 reference individuals"
    exit 1
fi

echo "===== VALIDATE 4000 SITES ====="

REGIONS_TMP=$(mktemp "$WORKDIR/results/.regions.XXXXXX")
TMP=$(mktemp "$WORKDIR/results/.reference700.XXXXXX")

trap 'rm -f "$REGIONS_TMP" "$TMP"' EXIT

awk '
BEGIN {OFS="\t"}

{
    if (NF != 2 || $2 !~ /^[0-9]+$/) {
        bad=1
        next
    }

    key=$1 FS $2

    if (seen[key]++) {
        duplicate=1
    }

    print $1,$2
}

END {
    if (NR != 4000 || bad || duplicate) {
        print "ERROR: Invalid, duplicated or incomplete site list" > "/dev/stderr"
        exit 1
    }
}
' "$SITES" > "$REGIONS_TMP"

echo "Validated 4000 unique positions"

echo "===== EXTRACT REFERENCE SCAFFOLD SNPs ====="

{
    printf "CHROM\tPOS\tREF\tALT\tAC\tAN\tAF\tMAF\n"

    bcftools view \
        -R "$REGIONS_TMP" \
        -Ou "$INPUT" \
    | bcftools +fill-tags \
        -Ou -- -t AC,AN,AF,MAF \
    | bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/MAF\n'

} > "$TMP"

N_OUT=$(awk 'END{print NR-1}' "$TMP")

echo "Extracted variants: $N_OUT"

if [[ "$N_OUT" -ne 4000 ]]; then
    echo "ERROR: Expected 4000 variants but obtained $N_OUT"
    echo "Check missing positions or multiple records per position"
    exit 1
fi

awk -F '\t' '
NR>1 {
    if (NF != 8) bad=1
    if ($8 != "." && ($8 < 0 || $8 > 0.5)) bad=1
}
END {exit bad}
' "$TMP"

mv "$TMP" "$OUT"

echo "===== REFERENCE RESULTS ====="
echo "Output: $OUT"

head -n 5 "$OUT"

echo "===== REFERENCE ANALYSIS FINISHED ====="
