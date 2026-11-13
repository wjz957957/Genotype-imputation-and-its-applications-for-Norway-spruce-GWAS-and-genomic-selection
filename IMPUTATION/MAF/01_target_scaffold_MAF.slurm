#!/bin/bash
#SBATCH -p shared
#SBATCH -A naiss2025-5-732
#SBATCH -J target347_MAF
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis/logs/%x_%j.err

set -euo pipefail

WORKDIR="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/05_scaffold_frequency_analysis"

INPUT="/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/IMPUTATION/MAF-impute5/Chr11_MAF-347-4K.vcf.gz"

OUT="$WORKDIR/results/target347_scaffold_MAF.tsv"

mkdir -p "$WORKDIR/results"

echo "===== LOAD SOFTWARE ====="

module load PDC/26.03
module load bioinfo-tools
ml bcftools/1.20

bcftools --version | head -n 1

echo "===== CHECK INPUT ====="

if [[ ! -r "$INPUT" ]]; then
    echo "ERROR: Input VCF not found"
    exit 1
fi

N_SAMPLES=$(bcftools query -l "$INPUT" | awk 'END{print NR+0}')
N_VARIANTS=$(bcftools view -H "$INPUT" | awk 'END{print NR+0}')

echo "Samples: $N_SAMPLES"
echo "Variants: $N_VARIANTS"

if [[ "$N_SAMPLES" -ne 347 ]]; then
    echo "ERROR: Expected 347 samples"
    exit 1
fi

if [[ "$N_VARIANTS" -ne 4000 ]]; then
    echo "ERROR: Expected 4000 variants"
    exit 1
fi

echo "===== CALCULATE AC, AN, AF AND MAF ====="

TMP=$(mktemp "$WORKDIR/results/.target347.XXXXXX")

trap 'rm -f "$TMP"' EXIT

{
    printf "CHROM\tPOS\tREF\tALT\tAC\tAN\tAF\tMAF\n"

    bcftools +fill-tags "$INPUT" \
        -Ou -- -t AC,AN,AF,MAF \
    | bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/MAF\n'

} > "$TMP"

N_OUT=$(awk 'END{print NR-1}' "$TMP")

if [[ "$N_OUT" -ne 4000 ]]; then
    echo "ERROR: Output does not contain 4000 variants"
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

echo "===== TARGET RESULTS ====="
echo "Output: $OUT"
echo "Variants: $N_OUT"

head -n 5 "$OUT"

echo "===== TARGET ANALYSIS FINISHED ====="
