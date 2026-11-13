#!/usr/bin/env Rscript
# Cross-family 5 x 3 inner CV; selected-chip + GWAS-prioritized QTN (ordinary markers).
# Fit algorithm is intentionally matched to Cross_family/scripts/05_run_inner_rrblup.R.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript 08b_run_inner_qtn_rrblup.R WORK OUTER INNER TRAIT")
WORK <- normalizePath(args[1], mustWork = TRUE)
OUTER <- args[2]; INNER <- args[3]; TRAIT <- args[4]
if (!OUTER %in% sprintf('%02d', 1:5) || !INNER %in% sprintf('%02d', 1:3) ||
    !TRAIT %in% c('BB','H6')) stop('Invalid outer / inner / trait')
if (!requireNamespace('rrBLUP', quietly = TRUE)) stop('rrBLUP package is unavailable')
COUNT <- c(0L,1L,2L,3L,4L,5L,10L,20L,30L,40L,50L,75L,100L,150L,200L,300L,400L,500L,600L,700L)
FOLD <- file.path(paste0('outer_', OUTER), paste0('inner_', INNER))
IN_DIR <- file.path(WORK, 'Cross_family/results/08_inner_qtn_genotypes', FOLD, TRAIT)
PHENO_DIR <- file.path(WORK, 'Cross_family/results/02_inner_phenotypes', FOLD)
QTN_DIR <- file.path(WORK, 'Cross_family/results/07_inner_crosschip_ld_filtering', FOLD, TRAIT, 'marker_sets')
CHIP_BASE <- file.path(WORK, 'Cross_family/results/05_inner_rrblup', FOLD, 'metrics.tsv')
OUT_DIR <- file.path(WORK, 'Cross_family/results/08_inner_qtn_rrblup', FOLD, TRAIT)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FLAG <- file.path(OUT_DIR, 'complete.flag')
if (file.exists(FLAG)) file.remove(FLAG)

require_nonempty <- function(p) {
  if (!file.exists(p) || is.na(file.info(p)$size) || file.info(p)$size <= 0L)
    stop('Missing or empty: ', p)
}
read_ids <- function(p, allow_empty = FALSE) {
  if (!file.exists(p)) stop('Missing: ', p)
  x <- trimws(readLines(p, warn = FALSE))
  if (any(!nzchar(x))) stop('Blank row in ', p)
  if (!allow_empty && length(x) == 0L) stop('Empty ID list ', p)
  if (anyDuplicated(x)) stop('Duplicate IDs in ', p)
  x
}
read_pheno <- function(p) {
  require_nonempty(p)
  x <- read.table(p, header = FALSE, sep = '', colClasses = 'character',
                  quote = '', comment.char = '', stringsAsFactors = FALSE)
  if (ncol(x) != 3L || nrow(x) == 0L || any(x[[1]] != x[[2]]) ||
      anyDuplicated(x[[2]])) stop('Invalid phenotype file ', p)
  y <- suppressWarnings(as.numeric(x[[3]]))
  if (any(!is.finite(y))) stop('Non-finite phenotype in ', p)
  data.frame(id = x[[2]], y = y, stringsAsFactors = FALSE)
}
read_raw <- function(p) {
  require_nonempty(p)
  if (requireNamespace('data.table', quietly = TRUE)) {
    x <- data.table::fread(p, header=TRUE, check.names=FALSE, data.table=FALSE,
                            na.strings = c('NA', 'NaN'))
  } else {
    x <- read.table(p, header = TRUE, check.names=FALSE, sep='',
                    na.strings=c('NA','NaN'), quote='', comment.char='', fill=FALSE)
  }
  if (ncol(x) < 7L || !identical(names(x)[1:2], c('FID', 'IID')) ||
      any(as.character(x[[1]]) != as.character(x[[2]])) || anyDuplicated(x[[2]]))
    stop('Invalid PLINK RAW: ', p)
  col_heads <- names(x)[-(1:6)]
  ids <- sub('_[^_]+$', '', col_heads)
  if (anyDuplicated(ids)) stop('Duplicate SNP IDs in RAW ', p)
  G <- as.matrix(x[, -(1:6), drop=FALSE]); storage.mode(G) <- 'double'
  if (any(is.infinite(G)) || any(G < 0 | G > 2, na.rm=TRUE))
    stop('Genotype dosage outside [0,2] or infinite in ', p)
  rownames(G) <- as.character(x[[2]])
  list(ids=as.character(x[[2]]), headers=col_heads, snps=ids, G=G)
}

a <- readLines(file.path(IN_DIR, 'complete.flag'), warn = FALSE)
if (!all(c('status=PASS', 'train_valid_snp_header_match=PASS',
           'outer_test_excluded=PASS', 'training_only_QTN_selection=PASS') %in% a))
  stop('Stage08a completion flag invalid')
chip <- read_ids(file.path(IN_DIR, 'selected_chip.snps'))
top <- read_ids(file.path(IN_DIR, 'extra_700.snps'))
if (length(top) != 700L || length(intersect(chip, top)) > 0L)
  stop('Invalid extra 700 or chip/QTN overlap')
for (k in COUNT) {
  fk <- file.path(QTN_DIR, paste0('extra_', k, '.snps'))
  ids <- read_ids(fk, allow_empty = (k == 0L))
  if (!identical(ids, head(top, k))) stop('Nested QTN set does not match top700: ', fk)
}

TRAIN <- read_raw(file.path(IN_DIR, 'train.raw'))
VALID <- read_raw(file.path(IN_DIR, 'valid.raw'))
if (!identical(TRAIN$headers, VALID$headers)) stop('Training/validation allele headers differ')
if (!setequal(TRAIN$snps, c(chip, top)) || length(TRAIN$snps) != length(chip) + 700L)
  stop('RAW marker set differs from selected chip plus QTN top700')
ytr <- read_pheno(file.path(PHENO_DIR, paste0(TRAIT,'.train.pheno')))
yva <- read_pheno(file.path(PHENO_DIR, paste0(TRAIT,'.valid.pheno')))
if (!setequal(TRAIN$ids, ytr$id) || !setequal(VALID$ids, yva$id) ||
    anyDuplicated(ytr$id) || anyDuplicated(yva$id)) stop('Phenotype/genotype sample mismatch')
TRAIN_INDEX <- match(ytr$id, TRAIN$ids)
VALID_INDEX <- match(yva$id, VALID$ids)
if (anyNA(TRAIN_INDEX) || anyNA(VALID_INDEX)) stop('Failed to align genotype rows')
# Reorder all variants in every model deterministically: chip, then QTN rank order.
IDX <- match(c(chip, top), TRAIN$snps)
if (anyNA(IDX)) stop('Missing markers in RAW matrix')
Xtr <- TRAIN$G[TRAIN_INDEX, IDX, drop=FALSE]
Xva <- VALID$G[VALID_INDEX, IDX, drop=FALSE]
rm(TRAIN, VALID); invisible(gc())

# Strict baseline validation against the original Stage05 chip-only fit.
require_nonempty(CHIP_BASE)
old <- read.delim(CHIP_BASE, sep='\t', header=TRUE, check.names=FALSE,
                  stringsAsFactors=FALSE)
base_config <- read.delim(file.path(WORK, 'Cross_family/config/best_chip_by_outer_trait.tsv'),
                          sep='\t', header=TRUE, stringsAsFactors=FALSE,
                          colClasses='character')
best <- base_config[base_config$Outer == OUTER & base_config$Trait == TRAIT, 'Best_Chip']
if (length(best) != 1L || !identical(as.integer(sub('^chip_', '', best)), length(chip)))
  stop('Selected chip strategy is inconsistent with chip SNP count')
old_base <- old[old$trait==TRAIT & old$marker_set==best, ,drop=FALSE]
if (nrow(old_base) != 1L || old_base$n_train != length(ytr$y) ||
    old_base$n_valid != length(yva$y)) stop('Stage05 baseline row missing/mismatched')

metrics <- vector('list', length(COUNT)); preds <- vector('list', length(COUNT))
base_r <- base_rmse <- NA_real_
for (i in seq_along(COUNT)) {
  k <- COUNT[i]
  index <- seq_len(length(chip) + k)
  xt <- Xtr[, index, drop=FALSE]
  xv <- Xva[, index, drop=FALSE]
  xt[!is.finite(xt)] <- NA_real_
  xv[!is.finite(xv)] <- NA_real_
  # Exactly follow original 05_run_inner_rrblup.R: train-mean imputation,
  # drop training-monomorphic markers, then train-centered dosage matrix.
  means <- colMeans(xt, na.rm = TRUE)
  keep <- is.finite(means)
  xt <- xt[, keep, drop=FALSE]
  xv <- xv[, keep, drop=FALSE]
  means <- means[keep]
  if (ncol(xt) == 0L) stop('No finite training means for k=', k)
  for (j in seq_len(ncol(xt))) {
    m1 <- !is.finite(xt[,j]); m2 <- !is.finite(xv[,j])
    if (any(m1)) xt[m1,j] <- means[j]
    if (any(m2)) xv[m2,j] <- means[j]
  }
  std <- apply(xt, 2L, sd)
  keep_var <- is.finite(std) & std > 0
  # Determine chip/QTN retained columns in the original ordered union.
  used_full <- which(keep)[keep_var]
  n_chip_used <- sum(used_full <= length(chip))
  n_extra_used <- sum(used_full > length(chip))
  xt <- xt[, keep_var, drop=FALSE]
  xv <- xv[, keep_var, drop=FALSE]
  if (ncol(xt) == 0L) stop('All training markers monomorphic for k=',k)
  center <- colMeans(xt)
  xt <- sweep(xt, 2L, center, '-')
  xv <- sweep(xv, 2L, center, '-')
  cat('Fitting ',TRAIT,' extra_',k,' | requested ',length(chip)+k,
      ' markers | used ',ncol(xt),'\n',sep='')
  fit <- rrBLUP::mixed.solve(y=ytr$y, Z=xt,
                             X=matrix(1,nrow=length(ytr$y),ncol=1L))
  prediction <- as.numeric(fit$beta[1L] + xv %*% fit$u)
  if (any(!is.finite(prediction))) stop('Invalid predictions for k=',k)
  r <- if (length(yva$y)>1L && sd(yva$y)>0 && sd(prediction)>0)
          as.numeric(cor(yva$y,prediction,method='pearson')) else NA_real_
  if (!is.finite(r)) stop('Undefined Pearson r for k=',k)
  rmse <- sqrt(mean((prediction - yva$y)^2))
  mae <- mean(abs(prediction - yva$y))
  bias <- mean(prediction - yva$y)
  vu <- as.numeric(fit$Vu); ve <- as.numeric(fit$Ve)
  lambda <- if (is.finite(vu) && vu>0) ve/vu else NA_real_
  if (k == 0L) {
    if (ncol(xt) != old_base$n_markers_used ||
        abs(r-old_base$pearson_r) > 1e-4 ||
        abs(rmse-old_base$RMSE) > 1e-4 * max(1, abs(old_base$RMSE))) {
      stop(sprintf('Stage05 chip-only baseline mismatch: new r %.8f old r %.8f; new RMSE %.8f old %.8f; used %d old %d',
                   r,old_base$pearson_r,rmse,old_base$RMSE,ncol(xt),old_base$n_markers_used))
    }
    base_r <- r; base_rmse <- rmse
    cat('PASS: extra_0 reproduces Stage05 chip-only baseline\n')
  }
  metrics[[i]] <- data.frame(outer=OUTER, inner=INNER, trait=TRAIT, chip_strategy=best,
      n_qtn=k, n_train=length(ytr$y), n_valid=length(yva$y),
      n_chip_requested=length(chip), n_markers_requested=length(chip)+k,
      n_markers_used=ncol(xt), n_chip_used=n_chip_used,
      n_qtn_used=n_extra_used, pearson_r=r, RMSE=rmse,
      MAE=mae, bias=bias, Vu=vu, Ve=ve, lambda=lambda,
      delta_r=r-base_r, delta_RMSE=rmse-base_rmse, status='PASS')
  preds[[i]] <- data.frame(outer=OUTER,inner=INNER,trait=TRAIT,n_qtn=k,
      sample_id=yva$id,observed=yva$y,predicted=prediction)
  cat(sprintf('PASS extra_%s r=%.6f delta_r=%.6f RMSE=%.6f\n',k,r,r-base_r,rmse))
  rm(xt,xv,fit,prediction)
  invisible(gc())
}
ALL_METRICS <- do.call(rbind,metrics)
ALL_PRED <- do.call(rbind,preds)
if (nrow(ALL_METRICS)!=length(COUNT) || anyDuplicated(ALL_METRICS$n_qtn) ||
    any(ALL_METRICS$status != 'PASS') || nrow(ALL_PRED) != length(COUNT)*nrow(yva))
  stop('Final model/record checks failed')
write_atomic <- function(obj,p) {
  tmp <- paste0(p,'.tmp')
  write.table(obj,file=tmp,sep='\t',quote=FALSE,row.names=FALSE)
  if (!file.rename(tmp,p)) stop('Cannot publish ',p)
}
write_atomic(ALL_METRICS, file.path(OUT_DIR,'metrics.tsv'))
write_atomic(ALL_PRED, file.path(OUT_DIR,'validation_predictions.tsv'))
flag_tmp <- paste0(FLAG,'.tmp')
writeLines(c('version=1',paste0('outer=',OUTER),paste0('inner=',INNER),
             paste0('trait=',TRAIT),paste0('chip_strategy=',best),
             'models=20','n_qtn_max=700','baseline_reproduced=PASS',
             'training_only_imputation_and_centering=PASS',
             'validation_prediction_check=PASS','status=PASS'),flag_tmp)
if (!file.rename(flag_tmp,FLAG)) stop('Cannot publish completion flag')
cat('OVERALL PASS outer=',OUTER,' inner=',INNER,' trait=',TRAIT,' models=20\n',sep='')
