Figure 4 additional baselines (ABLUP, RRBLUP-array / GBLUP-K, RRBLUP-imputed / GBLUP-K)
========================================================================================

What we verified from uploads
-----------------------------
- phe-ag-blup.csv: 904 individuals, each has BB=Adj_Budburst and H6=Adj_H6.
- pedigree.csv: 953 individuals, comprising 49 founder parents and 904 offspring.
- The 904 phenotype IDs match exactly the 904 offspring IDs in the pedigree.
- The G matrices were not uploaded and could NOT be validated here; the script will validate them on the server.

Important G matrix definitions
------------------------------
- Snparray-plink-G-all.csv must be based on the un-imputed SNP array marker set.
- Combine-plink-G-all.csv MUST be checked to verify whether it contains ONLY the imputed marker set,
  or array+imputed SNPs combined, and whether it is made using test-sample genotypes.
- This analysis fits a K-matrix version of GBLUP (the relationship-kernel equivalent of standard
  ridge-regression BLUP under compatible SNP coding/scaling). Manuscript naming should be clear.
- Do not reconstruct outer folds! Use exactly the same original 02_CV/outer_01..outer_10
  outer_train.samples and outer_test.samples as the nested-CV analysis.

Run on your new server (example)
--------------------------------
cd ~/GS-904-final
# Copy these files into this folder:
# 20_run_matrix_baselines_nestedcv.R
# phe-ag-blup.csv
# pedigree.csv
# Snparray-plink-G-all.csv
# Combine-plink-G-all.csv
# The original 02_CV/outer_01..outer_10 sample lists (or CV/outer_XX/)

# 1) First verify the input files and outer folds (will not fit models):
Rscript 20_run_matrix_baselines_nestedcv.R --work "$PWD" --validate-only

# 2) Run all baseline models:
mkdir -p logs
nohup Rscript 20_run_matrix_baselines_nestedcv.R --work "$PWD" \
  > logs/matrix_baselines_outer10.log 2>&1 &
tail -n 35 logs/matrix_baselines_outer10.log

# 3) Read final output:
column -t matrix_baselines_outer10/matrix_baselines_summary.tsv
column -t matrix_baselines_outer10/matrix_baselines_outer_metrics.tsv | head -n 15

# If G matrices are not available yet, ABLUP can run by itself:
Rscript 20_run_matrix_baselines_nestedcv.R --work "$PWD" --models ABLUP

# If sample lists have another root, add --cv-dir /absolute/path/to/02_CV
# If the matrix files live elsewhere, add --g-chip /path/to/Snparray-plink-G-all.csv \
#   --g-imputed /path/to/Combine-plink-G-all.csv
# If using an already-computed pedigree A matrix, pass --a /path/to/A-matrix.csv
# If inputs change after a partial run, use --force to recompute rather than reuse existing outputs.

Output description
------------------
- matrix_baselines_outer10/matrix_baselines_summary.tsv: 2 traits x 3 baseline models, 10-fold means/SD
- matrix_baselines_outer10/matrix_baselines_outer_metrics.tsv: all 60 outer-fold metric records
- matrix_baselines_outer10/per_fold/{BB,H6}/outer_XX/*_predictions.tsv: outer-test predictions
- matrix_baselines_outer10/per_fold/{BB,H6}/outer_XX/*_metrics.tsv: fold metrics
- RMSE_raw/MAE_raw/Bias_raw correspond to unstandardized adjusted phenotypes.
- RMSE_trainSD/MAE_trainSD/Bias_trainSD divide error by the TRAINING fold's phenotype SD.
- Earlier nested-CV DeepGS/RR-BLUP RMSE may have used a different/global phenotype Z scale,
  so do not compare absolute errors with those outputs without matching the phenotype scale.
- Pearson/Spearman correlations are invariant to a positive linear phenotype scale change.
- Predictions use raw phenotype units; if needed, later convert the prediction outputs to the
  exact scale used by the other Figure 4 analyses before comparing absolute-error metrics.

Performance method
------------------
For each model, each trait, each outer fold:
1) fit rrBLUP::mixed.solve() using y_train and K[train, train] only;
2) estimate REML variance ratio Ve/Vg in the outer training individuals;
3) predict test phenotypes with K[test, train] cross-covariances;
4) compute all metrics on unseen test phenotypes.
This is a standard kernel BLUP predictor with no inner hyperparameter selection.

Input validation
----------------
- Test samples across outer_01..outer_10 must be a disjoint partition of the 904 phenotype IDs.
- Every train and test pair must partition the full 904 IDs.
- G matrices are read as 905x905 CSV text INCLUDING one ID header row and one ID column.
  Numerically they must be 904x904; row and column IDs must match the phenotype sample IDs.
- G and pedigree A matrices must be symmetric, finite, with positive diagonal; G matrices
  will be tested for positive semidefiniteness up to minor rounding error.

Important methodological caveat
------------------------------
The script reads a precomputed, cohort-wide G matrix. It is acceptable for a *genotype-only*
transductive genomic prediction benchmark when the genotype information on test subjects is
available at prediction time; HOWEVER its precomputation might use allele frequencies, genotype
imputation, filtering or other transformations estimated with test genotypes. This is not the same
as rebuilding every data-dependent preprocessing step from the training subset. State this clearly
in Methods, and where necessary compute training-fold-specific G matrices for strictest comparison.
