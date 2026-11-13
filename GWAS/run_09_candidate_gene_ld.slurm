#!/bin/bash
#SBATCH -p shared
#SBATCH -A naiss2025-5-732
#SBATCH --mail-user=253721079@qq.com
#SBATCH --job-name=GWAS_genebody_LD
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=GWAS_genebody_LD_%j.out
#SBATCH --error=GWAS_genebody_LD_%j.err
set -euo pipefail

ml bioinfo-tools plink/1.90b4.9
ROOT=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer
NEW="$ROOT/GWAS_EDITOR_RESPONSE/03_LD_analysis/reanalysis_v1"
HERE="$NEW/scripts"
OUT="$NEW/03_ld_gene_body/09_gene_lead_LD_final"
BB_BFILE="${BB_BFILE:-$ROOT/GWAS/GWAS-Bud/Bud-imputed-plink}"

# DBH prefix must be verified. Attempt discovery only if unambiguous.
if [[ -z "${DBH_BFILE:-}" ]]; then
    mapfile -t candidates < <(find "$ROOT/GWAS" -maxdepth 6 -type f -iname '*.bim' \
      | grep -Ei 'dbh|diameter' | grep -Ei 'imput' || true)
    if [[ "${#candidates[@]}" -ne 1 ]]; then
        echo "ERROR: Cannot unambiguously identify the imputed DBH PLINK prefix." >&2
        printf 'DBH candidates found:\n%s\n' "${candidates[@]}" >&2
        echo 'Submit with: sbatch --export=ALL,DBH_BFILE=/full/verified/prefix scripts/run_09_candidate_gene_ld.slurm' >&2
        exit 2
    fi
    DBH_BFILE="${candidates[0]%.bim}"
fi

command -v python3 >/dev/null || { echo 'python3 not in PATH' >&2; exit 2; }
command -v plink >/dev/null || { echo 'PLINK 1.9 not in PATH; load the same PLINK module used previously' >&2; exit 2; }

mkdir -p "$OUT" "$NEW/04_functional_annotation"
cmd=(
  python3 "$HERE/run_09_candidate_gene_ld.py"
  --seed "$HERE/Supplementary_Table_S16_input.tsv"
  --candidates "$NEW/02_gene_annotation_fresh/imputed_specific_genes.tsv"
  --regions "$NEW/01_regions/all_index_snps.tsv"
  --bb-bfile "$BB_BFILE"
  --dbh-bfile "$DBH_BFILE"
  --outdir "$OUT"
  --plink plink
  --threads "${SLURM_CPUS_PER_TASK:-4}"
  --memory 16000
)
if [[ -n "${BB_ORIGINAL_BIM:-}" ]]; then
    cmd+=(--bb-original-bim "$BB_ORIGINAL_BIM")
fi
if [[ -n "${DBH_ORIGINAL_BIM:-}" ]]; then
    cmd+=(--dbh-original-bim "$DBH_ORIGINAL_BIM")
fi
"${cmd[@]}"

cp "$OUT/Supplementary_Table_S16_with_LD.tsv" "$NEW/04_functional_annotation/"
cp "$OUT/Manuscript_LD_insertions.txt" "$NEW/04_functional_annotation/"

echo 'All nine lead SNP to gene-body LD comparisons completed; S16 TSV and draft text updated.'
