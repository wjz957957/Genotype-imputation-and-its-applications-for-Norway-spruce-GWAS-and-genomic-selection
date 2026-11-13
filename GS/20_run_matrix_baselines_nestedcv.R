#!/usr/bin/env Rscript
# Additional ABLUP / GBLUP (RR-BLUP-equivalent relationship kernel) baselines.
# These baselines use the EXACT same 10 held-out outer folds as the nested-CV analysis.
# Only training phenotypes enter model fitting; outer test phenotypes are read after fitting.
# No inner CV is needed for these fixed, all-marker baselines (no marker-count tuning).

suppressPackageStartupMessages(library(rrBLUP))

help_message <- paste(
  "Usage:",
  "  Rscript 20_run_matrix_baselines_nestedcv.R [options]",
  "Options:",
  "  --work DIRECTORY            Project root (default: current directory)",
  "  --cv-dir DIRECTORY          Folder containing outer_01/...outer_{train,test}.samples",
  "  --pheno FILE                phe-ag-blup.csv (default: <work>/phe-ag-blup.csv)",
  "  --ped FILE                  pedigree.csv (default: <work>/pedigree.csv)",
  "  --a FILE                    Optional existing pedigree A matrix (otherwise derive from pedigree.csv)",
  "  --g-chip FILE               Whole-array G (default: <work>/Snparray-plink-G-all.csv)",
  "  --g-imputed FILE            Imputed dataset G (default: <work>/Combine-plink-G-all.csv)",
  "  --out DIRECTORY             Output directory (default: <work>/matrix_baselines_outer10)",
  "  --models CSV                Model subset: ABLUP,RRBLUP-array,RRBLUP-imputed",
  "  --validate-only             Validate inputs and sample folds, without model fitting",
  "  --force                     Refit models even when outputs already exist",
  sep = "\n"
)

parse_args <- function(a) {
  out <- list()
  i <- 1L
  while (i <= length(a)) {
    token <- a[i]
    if (token == "--help") { cat(help_message, "\n"); quit(save = "no", status = 0L) }
    if (token %in% c("--validate-only", "--force")) {
      out[[sub("^--", "", token)]] <- TRUE
      i <- i + 1L
      next
    }
    if (!grepl("^--", token)) stop("Unknown argument: ", token)
    key <- sub("^--", "", token)
    if (grepl("=", key, fixed = TRUE)) {
      bits <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- bits[1]
      val <- paste(bits[-1], collapse = "=")
    } else {
      if (i == length(a) || grepl("^--", a[i+1L])) stop("Value missing for ", token)
      i <- i + 1L
      val <- a[i]
    }
    out[[key]] <- val
    i <- i + 1L
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
work <- normalizePath(if (is.null(opt$work)) getwd() else opt$work, mustWork = TRUE)
getpath <- function(name, default) {
  x <- opt[[name]]
  if (is.null(x)) x <- file.path(work, default)
  path.expand(x)
}
pheno_path <- getpath("pheno", "phe-ag-blup.csv")
ped_path   <- getpath("ped", "pedigree.csv")
a_path     <- if (is.null(opt$a)) NULL else path.expand(opt$a)
g_chip     <- getpath("g-chip", "Snparray-plink-G-all.csv")
g_imp      <- getpath("g-imputed", "Combine-plink-G-all.csv")
outdir     <- getpath("out", "matrix_baselines_outer10")

models <- if (is.null(opt$models)) c("ABLUP", "RRBLUP-array", "RRBLUP-imputed") else {
  trimws(strsplit(opt$models, ",", fixed=TRUE)[[1]])
}
allowed <- c("ABLUP", "RRBLUP-array", "RRBLUP-imputed")
if (!length(models) || anyDuplicated(models) || !all(models %in% allowed)) {
  stop("--models must contain unique choices from: ", paste(allowed, collapse = ","))
}
validate_only <- isTRUE(opt[["validate-only"]])
force <- isTRUE(opt$force)

cvdir <- opt[["cv-dir"]]
if (is.null(cvdir)) {
  choices <- c(file.path(work, "02_CV"), file.path(work, "CV"))
  matches <- choices[file.exists(file.path(choices, "outer_01", "outer_train.samples"))]
  if (!length(matches)) {
    stop("CV folds not found. Copy the original '02_CV/outer_XX/' folder into this project ",
         "or pass --cv-dir PATH. DO NOT regenerate fold partitions.")
  }
  cvdir <- matches[1]
}
cvdir <- normalizePath(cvdir, mustWork = TRUE)

read_ids <- function(file) {
  if (!file.exists(file)) stop("CV ID list not found: ", file)
  x <- trimws(readLines(file, warn = FALSE, encoding = "UTF-8"))
  x <- x[nzchar(x)]
  if (!length(x) || anyDuplicated(x)) stop("Empty or duplicated sample IDs: ", file)
  x
}

if (!file.exists(pheno_path)) stop("Phenotype file missing: ", pheno_path)
pheno <- read.csv(pheno_path, stringsAsFactors = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")
need <- c("Sample.Filename", "Adj_Budburst", "Adj_H6")
if (!all(need %in% names(pheno))) stop("Phenotype header should include: ", paste(need, collapse = ", "))
ids <- as.character(pheno[["Sample.Filename"]])
if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Phenotype IDs invalid or duplicated")
if (length(ids) != 904L) warning("Phenotype file contains ", length(ids), " individuals (expected 904)")
Y <- list(BB = as.numeric(pheno$Adj_Budburst), H6 = as.numeric(pheno$Adj_H6))
if (any(vapply(Y, function(x) any(!is.finite(x)), logical(1)))) stop("Phenotype contains missing / non-numeric values")
names(Y$BB) <- ids
names(Y$H6) <- ids

folds <- list()
for (i in 1:10) {
  fold <- sprintf("%02d", i)
  dir_i <- file.path(cvdir, paste0("outer_",fold))
  train <- read_ids(file.path(dir_i, "outer_train.samples"))
  test  <- read_ids(file.path(dir_i, "outer_test.samples"))
  if (length(intersect(train, test))) stop("Training/test overlap in fold ",fold)
  if (!setequal(c(train, test), ids)) stop("Fold ", fold, " does not partition the entire phenotype cohort")
  folds[[fold]] <- list(train=train, test=test)
  cat("Outer ", fold, ": train=", length(train), " test=", length(test), "\n", sep="")
}
test_pool <- unlist(lapply(folds, function(x) x$test), use.names = FALSE)
if (length(test_pool) != length(ids) || anyDuplicated(test_pool) || !setequal(test_pool, ids)) {
  stop("Outer folds are NOT standard disjoint 10-fold CV: test subjects must each appear exactly once")
}
cat("CV check PASS: all ", length(ids), " individuals tested once.\n", sep="")

read_relation <- function(file, label) {
  if (!file.exists(file)) stop(label, " matrix file missing: ", file)
  cat("Loading ", label, " matrix from ", file, "\n", sep="")
  # Expected file: one header row, first column = individual ID, 904 numeric columns.
  raw <- read.csv(file, row.names=1, header=TRUE, check.names=FALSE,
                  stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
  K <- as.matrix(raw)
  suppressWarnings(storage.mode(K) <- "double")
  if (nrow(K) != ncol(K) || !all(is.finite(K))) {
    stop(label, ": not a square finite numeric matrix. Check CSV header and leading ID column")
  }
  if (anyDuplicated(rownames(K)) || anyDuplicated(colnames(K))) stop(label, ": duplicated matrix IDs")
  if (!setequal(rownames(K),colnames(K))) stop(label, ": row/column IDs differ")
  if (!setequal(rownames(K),ids)) {
    missing <- head(setdiff(ids,rownames(K)), 3L)
    extra <- head(setdiff(rownames(K),ids), 3L)
    stop(label, ": matrix ID set differs from phenotype. Missing: ", paste(missing,collapse=","),
         "; extra: ",paste(extra,collapse=","))
  }
  K <- K[ids,ids,drop=FALSE]
  maxasym <- max(abs(K-t(K)))
  if (maxasym > 1e-6 * max(1, max(abs(K)))) stop(label, ": asymmetric matrix, max deviation=",maxasym)
  K <- (K+t(K))/2
  if (any(diag(K)<=0)) stop(label, ": nonpositive diagonal")
  # Check PSD once per input; small roundoff negatives are tolerated.
  ev <- eigen(K, symmetric=TRUE, only.values=TRUE)$values
  if (min(ev) < -1e-5 * max(1,max(ev))) {
    stop(label, ": appreciably negative eigenvalue, minimum=", min(ev),
         "; wrong matrix type / input format or matrix not PSD")
  }
  cat(label, ": ", nrow(K), " samples; eigenvalue min=",signif(min(ev),4),
      "; max asymmetry=",signif(maxasym,4),"\n",sep="")
  K
}

make_pedigree_A <- function(file, wanted_ids) {
  if (!file.exists(file)) stop("Pedigree missing: ",file)
  p <- read.csv(file, stringsAsFactors=FALSE, check.names=FALSE,
                colClasses="character",fileEncoding="UTF-8-BOM")
  if (!all(c("ID","Sire","Dam") %in% names(p))) stop("Pedigree needs ID,Sire,Dam columns")
  p <- p[,c("ID","Sire","Dam")]
  if (anyNA(p$ID) || anyDuplicated(p$ID)) stop("Invalid or duplicate pedigree IDs")
  pp <- setNames(seq_len(nrow(p)), p$ID)
  parent_index <- function(z) {
    zero <- z %in% c("0","",NA_character_)
    r <- match(z, p$ID)
    r[zero] <- 0L
    if (anyNA(r)) stop("Pedigree references missing parents: ", paste(head(z[is.na(r)], 5),collapse=","))
    r
  }
  s <- parent_index(p$Sire)
  d <- parent_index(p$Dam)
  n <- nrow(p)
  if (any(s >= seq_len(n) | d >= seq_len(n))) {
    stop("Pedigree parents must precede offspring. Sort parents before descendants first")
  }
  cat("Building pedigree A: ",n," rows (including ",sum(s==0 & d==0)," founders)\n",sep="")
  A <- matrix(0,n,n)
  for (i in seq_len(n)) {
    if (i > 1L) {
      j <- seq_len(i-1L)
      v <- numeric(length(j))
      if (s[i] > 0L) v <- v + A[s[i],j] / 2
      if (d[i] > 0L) v <- v + A[d[i],j] / 2
      A[i,j] <- v
      A[j,i] <- v
    }
    A[i,i] <- 1 + if (s[i] > 0L && d[i] > 0L) A[s[i],d[i]]/2 else 0
  }
  rownames(A) <- colnames(A) <- p$ID
  if (!all(wanted_ids %in% p$ID)) stop("Pedigree missing phenotype IDs")
  A[wanted_ids,wanted_ids,drop=FALSE]
}

Kmodels <- list()
for (m in models) {
  if (m == "ABLUP") {
    Kmodels[[m]] <- if (is.null(a_path)) make_pedigree_A(ped_path, ids) else read_relation(a_path,"Pedigree A")
  } else if (m == "RRBLUP-array") {
    Kmodels[[m]] <- read_relation(g_chip, "SNP-array G")
  } else {
    Kmodels[[m]] <- read_relation(g_imp, "Imputed G")
  }
}
if (validate_only) {
  cat("All requested input files and fold partitions validated. No models fitted.\n")
  quit(save="no",status=0L)
}

# FIT ONLY y[train], to ensure no outer-test phenotype enters model fitting.
# rrBLUP::mixed.solve estimates intercept and REML variance ratio on K[train,train].
# Predict held-out subjects with the standard GBLUP kernel cross-covariance formula.
fit_predict <- function(K, train, test, phenotype) {
  yy <- unname(phenotype[train])
  kt <- K[train,train,drop=FALSE]
  solution <- rrBLUP::mixed.solve(y=yy, K=kt, method="REML")
  beta <- as.numeric(solution$beta)[1L]
  ratio <- as.numeric(solution$Ve/solution$Vu)
  if (length(ratio)!=1L || !is.finite(ratio) || ratio<=0) stop("Invalid REML variance ratio")
  V <- kt + diag(ratio, length(train))
  alpha <- as.numeric(solve(V, yy-beta))
  pred <- beta + as.numeric(K[test,train,drop=FALSE] %*% alpha)
  names(pred) <- test
  list(pred=pred, Vg=as.numeric(solution$Vu), Ve=as.numeric(solution$Ve), lambda=ratio)
}

compute_metrics <- function(obs,pred,train_phenotype,info) {
  error <- pred-obs
  sd_train <- sd(train_phenotype)
  if (!is.finite(sd_train) || sd_train<=0) stop("Invalid training phenotype SD")
  data.frame(
    Trait=info$trait, Outer=info$outer, Model=info$model,
    N_Train=info$ntrain, N_Test=info$ntest,
    Pearson=cor(obs,pred,method="pearson"),
    Spearman=cor(obs,pred,method="spearman"),
    RMSE_raw=sqrt(mean(error^2)), MAE_raw=mean(abs(error)), Bias_raw=mean(error),
    RMSE_trainSD=sqrt(mean(error^2))/sd_train,
    MAE_trainSD=mean(abs(error))/sd_train,
    Bias_trainSD=mean(error)/sd_train,
    Vg=info$fit$Vg, Ve=info$fit$Ve, Lambda=info$fit$lambda,
    stringsAsFactors=FALSE
  )
}

cat("\n===== Starting relationship-matrix BLUP benchmarks =====\n")
for (m in models) {
  K <- Kmodels[[m]]
  for (trait in c("BB","H6")) {
    yy <- Y[[trait]]
    for (outer in names(folds)) {
      fold <- folds[[outer]]
      out_fold <- file.path(outdir,"per_fold",trait,paste0("outer_",outer))
      dir.create(out_fold,recursive=TRUE,showWarnings=FALSE)
      metric_file <- file.path(out_fold,paste0(m,"_metrics.tsv"))
      pred_file <- file.path(out_fold,paste0(m,"_predictions.tsv"))
      if (!force && file.exists(metric_file) && file.exists(pred_file)) {
        cat("SKIP ",trait," outer_",outer," ",m," (already complete)\n",sep="")
        next
      }
      res <- fit_predict(K,fold$train,fold$test,yy)
      ob <- unname(yy[fold$test])
      pr <- unname(res$pred[fold$test])
      met <- compute_metrics(ob,pr,yy[fold$train],
                             list(trait=trait,outer=outer,model=m,
                                  ntrain=length(fold$train),ntest=length(fold$test),fit=res))
      predicted <- data.frame(Trait=trait,Outer=outer,Model=m,
                              ID=fold$test,Observed_raw=ob,Predicted_raw=pr,
                              Residual_pred_minus_obs=pr-ob,stringsAsFactors=FALSE)
      write.table(predicted,pred_file,sep="\t",quote=FALSE,row.names=FALSE)
      write.table(met,metric_file,sep="\t",quote=FALSE,row.names=FALSE)
      cat(sprintf("%-15s  %-2s  outer_%s  Pearson=%.5f  RMSE(train-SD)=%.4f\n",
                  m,trait,outer,met$Pearson,met$RMSE_trainSD))
    }
  }
}

all_metrics <- list()
for (trait in c("BB","H6")) for (outer in names(folds)) for (m in models) {
  p <- file.path(outdir,"per_fold",trait,paste0("outer_",outer),paste0(m,"_metrics.tsv"))
  if (!file.exists(p)) stop("Incomplete result: ",p)
  all_metrics[[length(all_metrics)+1L]] <- read.delim(p,check.names=FALSE)
}
all <- do.call(rbind,all_metrics)
if (nrow(all) != 10L * 2L * length(models)) stop("Unexpected metric record count")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.table(all,file.path(outdir,"matrix_baselines_outer_metrics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

summary <- list()
for (trait in c("BB","H6")) for (m in models) {
  x <- all[all$Trait==trait & all$Model==m,,drop=FALSE]
  z <- data.frame(Trait=trait,Model=m,Folds=nrow(x))
  for (col in c("Pearson","Spearman","RMSE_raw","MAE_raw","Bias_raw",
                "RMSE_trainSD","MAE_trainSD","Bias_trainSD")) {
    z[[paste0(col,"_mean")]] <- mean(x[[col]])
    z[[paste0(col,"_sd")]] <- sd(x[[col]])
  }
  summary[[length(summary)+1L]] <- z
}
summary <- do.call(rbind, summary)
write.table(summary,file.path(outdir,"matrix_baselines_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("\n===== DONE =====\n")
cat("Results: ",outdir,"\n",sep="")
print(summary[,c("Trait","Model","Pearson_mean","Pearson_sd","RMSE_trainSD_mean")],row.names=FALSE)
cat("NOTE: RMSE_raw is in raw phenotype units. RMSE_trainSD uses training-fold SD,\n")
cat("      which is not identical to RMSE from an existing global-Z phenotype file.\n")
