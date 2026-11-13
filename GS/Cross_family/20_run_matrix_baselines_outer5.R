#!/usr/bin/env Rscript
# Fit the three matrix-based baselines for the current 5-fold, parent-connected
# outer CV: ABLUP, RRBLUP-array, and RRBLUP-imputed (Combine-plink-G-all.csv).
# Uses fold-specific train phenotype files. Evaluation-only test phenotypes are
# opened only after each model has generated its held-out predictions.

suppressPackageStartupMessages(library(rrBLUP))

help_message <- paste(
  "Usage: Rscript 20_run_matrix_baselines_outer5.R [options]",
  "Options:",
  "  --work DIRECTORY            Project root (default: current directory)",
  "  --fold-assignment FILE      Current progeny-level outer-fold assignment TSV",
  "  --model-data DIRECTORY      10_outer_model_data root",
  "  --ped FILE                  pedigree.csv (default: <work>/pedigree.csv)",
  "  --a FILE                    Optional existing pedigree A matrix",
  "  --g-chip FILE               Whole-array G (default: <work>/Snparray-plink-G-all.csv)",
  "  --g-imputed FILE            Full imputed/combined G (default: <work>/Combine-plink-G-all.csv)",
  "  --out DIRECTORY             Output root (default: Cross_family/results/20_outer_matrix_baselines)",
  "  --models CSV                Subset of ABLUP,RRBLUP-array,RRBLUP-imputed",
  "  --validate-only             Validate all inputs without fitting models",
  "  --force                     Refit completed model tasks",
  sep = "\n"
)

parse_args <- function(a) {
  out <- list(); i <- 1L
  while (i <= length(a)) {
    token <- a[i]
    if (token == "--help") { cat(help_message, "\n"); quit(save="no", status=0L) }
    if (token %in% c("--validate-only", "--force")) {
      out[[sub("^--", "", token)]] <- TRUE; i <- i + 1L; next
    }
    if (!grepl("^--", token)) stop("Unknown argument: ", token)
    key <- sub("^--", "", token)
    if (grepl("=", key, fixed=TRUE)) {
      bits <- strsplit(key, "=", fixed=TRUE)[[1]]
      key <- bits[1]; val <- paste(bits[-1], collapse="=")
    } else {
      if (i == length(a) || grepl("^--", a[i+1L])) stop("Value missing for ", token)
      i <- i + 1L; val <- a[i]
    }
    out[[key]] <- val; i <- i + 1L
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly=TRUE))
work <- normalizePath(if (is.null(opt$work)) getwd() else opt$work, mustWork=TRUE)
getpath <- function(name, default) {
  x <- opt[[name]]
  if (is.null(x)) x <- file.path(work, default)
  path.expand(x)
}
assignment_path <- getpath("fold-assignment", "Cross_family/config/progeny_parent_connected_outer_fold_assignment.tsv")
model_data_root <- getpath("model-data", "Cross_family/results/10_outer_model_data")
ped_path <- getpath("ped", "pedigree.csv")
a_path <- if (is.null(opt$a)) NULL else path.expand(opt$a)
g_chip_path <- getpath("g-chip", "Snparray-plink-G-all.csv")
g_imp_path <- getpath("g-imputed", "Combine-plink-G-all.csv")
outdir <- getpath("out", "Cross_family/results/20_outer_matrix_baselines")
models <- if (is.null(opt$models)) c("ABLUP", "RRBLUP-array", "RRBLUP-imputed") else
  trimws(strsplit(opt$models, ",", fixed=TRUE)[[1]])
allowed <- c("ABLUP", "RRBLUP-array", "RRBLUP-imputed")
if (!length(models) || anyDuplicated(models) || !all(models %in% allowed)) {
  stop("--models must contain unique choices from: ", paste(allowed, collapse=","))
}
validate_only <- isTRUE(opt[["validate-only"]])
force <- isTRUE(opt$force)

need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)
require_file <- function(path, label=path) {
  need(file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0,
       paste("Missing or empty", label, ":", path))
  invisible(path)
}
read_ids <- function(path) {
  require_file(path, "ID file")
  x <- trimws(readLines(path, warn=FALSE, encoding="UTF-8")); x <- x[nzchar(x)]
  need(length(x)>0L && !anyDuplicated(x), paste("Empty/duplicate IDs:", path))
  x
}
read_flag <- function(path) {
  require_file(path, "completion flag")
  lines <- readLines(path, warn=FALSE); out <- list()
  for (line in lines) {
    p <- strsplit(line, "=", fixed=TRUE)[[1]]
    need(length(p)>=2L, paste("Malformed flag line:", path, line))
    k <- p[1]; v <- paste(p[-1], collapse="=")
    need(is.null(out[[k]]), paste("Duplicate flag key:", k, path))
    out[[k]] <- v
  }
  out
}
read_tsv <- function(path) {
  require_file(path, "TSV")
  read.delim(path, header=TRUE, sep="\t", check.names=FALSE,
             stringsAsFactors=FALSE, quote="", comment.char="")
}
write_tsv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  write.table(x, file=tmp, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  if (file.exists(path)) unlink(path)
  need(file.rename(tmp, path), paste("Could not publish", path))
}

# Validate the current 5-fold family-connected assignment, not the old 10-fold CV files.
require_file(assignment_path, "outer-fold assignment")
assignment <- read_tsv(assignment_path)
need(all(c("Sample.Filename", "Outer_test_fold", "Parent_connected_group") %in% names(assignment)),
     "Fold assignment must contain Sample.Filename, Outer_test_fold, Parent_connected_group")
all_ids <- as.character(assignment$Sample.Filename)
fold_num <- suppressWarnings(as.integer(as.character(assignment$Outer_test_fold)))
groups <- as.character(assignment$Parent_connected_group)
need(nrow(assignment)==904L, paste("Expected 904 assignment rows; found", nrow(assignment)))
need(!anyNA(all_ids) && all(nzchar(all_ids)) && !anyDuplicated(all_ids),
     "Invalid/duplicate Sample.Filename in assignment")
need(all(fold_num %in% 1:5), "Outer_test_fold must be 01..05")
need(!anyNA(groups) && all(nzchar(groups)), "Missing parent-connected group labels")
group_fold_counts <- tapply(fold_num, groups, function(x) length(unique(x)))
need(all(group_fold_counts==1L), "A parent-connected group is split across outer folds")
expected_test_n <- c(`01`=175L, `02`=178L, `03`=176L, `04`=185L, `05`=190L)
observed_test_n <- table(factor(sprintf("%02d", fold_num), levels=names(expected_test_n)))
need(all(as.integer(observed_test_n)==as.integer(expected_test_n)),
     paste("Unexpected outer-test counts:", paste(names(observed_test_n), observed_test_n, sep="=", collapse=", ")))
cat("Outer-fold assignment PASS: 904 individuals; test sizes 175/178/176/185/190; parent-connected groups intact.\n")

# Each fold/trait must have its own training phenotype and evaluation-only test phenotype.
fold_inputs <- list()
for (outer in sprintf("%02d", 1:5)) {
  expected_test <- all_ids[sprintf("%02d", fold_num)==outer]
  expected_train <- all_ids[sprintf("%02d", fold_num)!=outer]
  for (trait in c("BB", "H6")) {
    base <- file.path(model_data_root, paste0("outer_", outer), trait)
    flag <- read_flag(file.path(base, "complete.flag"))
    need(identical(flag$status, "PASS") && identical(flag$outer, outer) && identical(flag$trait, trait),
         paste("Input complete.flag does not match", outer, trait))
    train_ids <- read_ids(file.path(base, "train.ids"))
    test_ids <- read_ids(file.path(base, "test.ids"))
    need(setequal(train_ids, expected_train), paste("train.ids differ from fold assignment:", outer, trait))
    need(setequal(test_ids, expected_test), paste("test.ids differ from fold assignment:", outer, trait))
    need(length(intersect(train_ids, test_ids))==0L, paste("Train/test overlap:", outer, trait))
    train_pheno_path <- file.path(base, "train.pheno.tsv")
    test_pheno_path <- file.path(base, "evaluation_only", "test.pheno.tsv")
    require_file(train_pheno_path, "training phenotype")
    require_file(test_pheno_path, "evaluation-only test phenotype")
    train_pheno <- read_tsv(train_pheno_path)
    need(identical(names(train_pheno), c("FID", "IID", "phenotype")),
         paste("Unexpected training phenotype header:", train_pheno_path))
    need(nrow(train_pheno)==length(train_ids) && identical(as.character(train_pheno$FID), train_ids) &&
         identical(as.character(train_pheno$IID), train_ids),
         paste("Training phenotype ID/order mismatch:", outer, trait))
    y <- suppressWarnings(as.numeric(train_pheno$phenotype))
    need(all(is.finite(y)) && all(y != -9) && is.finite(sd(y)) && sd(y)>0,
         paste("Invalid training phenotype values:", outer, trait))
    # Only verify the evaluation-only file exists at this stage. Its values are not read until prediction.
    fold_inputs[[paste(outer, trait, sep="|")]] <- list(
      base=base, train_ids=train_ids, test_ids=test_ids,
      train_pheno_path=train_pheno_path, test_pheno_path=test_pheno_path,
      y_train=y
    )
  }
}
cat("Fold-specific phenotype inputs PASS: training values validated; evaluation-only values not read yet.\n")

read_relation <- function(file, label, wanted_ids) {
  require_file(file, paste(label, "matrix"))
  cat("Loading ", label, " from ", file, "\n", sep="")
  raw <- read.csv(file, row.names=1, header=TRUE, check.names=FALSE,
                  stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
  K <- as.matrix(raw); suppressWarnings(storage.mode(K) <- "double")
  need(nrow(K)==ncol(K) && all(is.finite(K)), paste(label, "is not a square finite numeric matrix"))
  need(!anyDuplicated(rownames(K)) && !anyDuplicated(colnames(K)), paste(label, "has duplicate matrix IDs"))
  need(setequal(rownames(K), colnames(K)), paste(label, "row/column IDs differ"))
  need(setequal(rownames(K), wanted_ids), paste(label, "ID set differs from 5-fold cohort"))
  K <- K[wanted_ids, wanted_ids, drop=FALSE]
  maxasym <- max(abs(K-t(K)))
  need(maxasym <= 1e-6*max(1, max(abs(K))), paste(label, "is asymmetric; max deviation=", maxasym))
  K <- (K+t(K))/2
  need(all(diag(K)>0), paste(label, "has nonpositive diagonal"))
  ev <- eigen(K, symmetric=TRUE, only.values=TRUE)$values
  need(min(ev) >= -1e-5*max(1, max(ev)), paste(label, "is not positive semidefinite; minimum eigenvalue=", min(ev)))
  cat(label, ": n=", nrow(K), "; min_eigenvalue=", signif(min(ev), 5),
      "; max_asymmetry=", signif(maxasym, 5), "\n", sep="")
  K
}

make_pedigree_A <- function(file, wanted_ids) {
  require_file(file, "pedigree")
  p <- read.csv(file, stringsAsFactors=FALSE, check.names=FALSE,
                colClasses="character", fileEncoding="UTF-8-BOM")
  need(all(c("ID", "Sire", "Dam") %in% names(p)), "Pedigree must contain ID,Sire,Dam columns")
  p <- p[, c("ID", "Sire", "Dam")]
  need(!anyNA(p$ID) && !anyDuplicated(p$ID), "Invalid/duplicate pedigree IDs")
  parent_index <- function(z) {
    zero <- is.na(z) | z %in% c("0", "")
    r <- match(z, p$ID); r[zero] <- 0L
    need(!anyNA(r), paste("Pedigree references missing parents:", paste(head(z[is.na(r)], 5), collapse=",")))
    r
  }
  s <- parent_index(p$Sire); d <- parent_index(p$Dam); n <- nrow(p)
  need(all(!(s >= seq_len(n) | d >= seq_len(n))),
       "Pedigree parents must precede offspring; sort parents before descendants first")
  cat("Building pedigree A from ", n, " pedigree records.\n", sep="")
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (i > 1L) {
      j <- seq_len(i-1L); v <- numeric(length(j))
      if (s[i] > 0L) v <- v + A[s[i], j]/2
      if (d[i] > 0L) v <- v + A[d[i], j]/2
      A[i, j] <- v; A[j, i] <- v
    }
    A[i, i] <- 1 + if (s[i] > 0L && d[i] > 0L) A[s[i], d[i]]/2 else 0
  }
  rownames(A) <- colnames(A) <- p$ID
  need(all(wanted_ids %in% p$ID), "Pedigree is missing one or more cohort IDs")
  A[wanted_ids, wanted_ids, drop=FALSE]
}

# Model phenotypes are deliberately sourced only from the fold-specific files above.
# Do not load the all-individual phenotype table here: the held-out phenotype remains
# in evaluation_only/test.pheno.tsv until each model has generated predictions.

Kmodels <- list()
if ("ABLUP" %in% models) {
  Kmodels$ABLUP <- if (is.null(a_path)) make_pedigree_A(ped_path, all_ids) else
    read_relation(a_path, "Pedigree A", all_ids)
  # Validate symmetry/PSD after construction, just as for supplied matrices.
  A <- Kmodels$ABLUP
  need(all(is.finite(A)) && max(abs(A-t(A))) < 1e-8, "Constructed pedigree A is invalid/asymmetric")
  ev <- eigen(A, symmetric=TRUE, only.values=TRUE)$values
  need(min(ev) >= -1e-5*max(1, max(ev)), "Constructed pedigree A is not PSD")
  cat("Pedigree A PASS: n=", nrow(A), "; min_eigenvalue=", signif(min(ev), 5), "\n", sep="")
}
if ("RRBLUP-array" %in% models) Kmodels[["RRBLUP-array"]] <- read_relation(g_chip_path, "SNP-array G", all_ids)
if ("RRBLUP-imputed" %in% models) Kmodels[["RRBLUP-imputed"]] <- read_relation(g_imp_path, "Full imputed/combined G", all_ids)

if (validate_only) {
  cat("VALIDATION PASS: all requested matrix inputs and current 5-fold inputs validated; no models fitted.\n")
  quit(save="no", status=0L)
}

fit_predict <- function(K, train, test, y_train) {
  K_train <- K[train, train, drop=FALSE]
  fit <- rrBLUP::mixed.solve(y=unname(y_train), K=K_train, method="REML")
  beta <- as.numeric(fit$beta)[1L]
  vu <- as.numeric(fit$Vu); ve <- as.numeric(fit$Ve)
  need(length(vu)==1L && length(ve)==1L && is.finite(vu) && is.finite(ve) && vu>0 && ve>0,
       "Invalid variance estimates from rrBLUP::mixed.solve")
  lambda <- ve/vu
  need(is.finite(lambda) && lambda>0, "Invalid Ve/Vu ratio")
  V <- K_train + diag(lambda, length(train))
  alpha <- as.numeric(solve(V, unname(y_train-beta)))
  pred <- beta + as.numeric(K[test, train, drop=FALSE] %*% alpha)
  need(length(pred)==length(test) && all(is.finite(pred)), "Non-finite held-out predictions")
  list(pred=pred, Vu=vu, Ve=ve, Lambda=lambda)
}
safe_cor <- function(obs, pred, method="pearson") {
  # A predictor with zero (or purely round-off) variation has undefined correlation.
  # Never replace its correlation with zero or fabricate a small jitter.
  if (length(obs)<2L || !all(is.finite(obs)) || !all(is.finite(pred)) ||
      !is.finite(sd(obs)) || sd(obs)<=1e-12 ||
      (max(pred)-min(pred)) <= 1e-10) return(NA_real_)
  as.numeric(cor(obs, pred, method=method))
}
compute_metrics <- function(obs, pred, y_train, info) {
  err <- pred-obs; train_sd <- sd(y_train)
  need(is.finite(train_sd) && train_sd>0, "Invalid training phenotype SD")
  pearson <- safe_cor(obs, pred, "pearson")
  spearman <- safe_cor(obs, pred, "spearman")
  pred_range <- max(pred)-min(pred)
  no_pedigree_link <- info$model == "ABLUP" && info$max_cross <= 1e-12
  degenerate <- !is.finite(pearson) || !is.finite(spearman)
  if (degenerate) {
    need(no_pedigree_link && pred_range <= 1e-10,
      paste("Unexpected constant/invalid prediction (not explained by zero pedigree cross-relationships):",
            info$model, info$trait, info$outer,
            "max_abs_K_test_train=", format(info$max_cross, digits=8),
            "prediction_range=", format(pred_range, digits=8)))
  }
  status <- if (degenerate) "UNDEFINED_NO_PEDIGREE_LINK" else "DEFINED"
  data.frame(
    Trait=info$trait, Outer=as.integer(info$outer), Model=info$model,
    N_Train=length(y_train), N_Test=length(obs),
    Pearson=pearson, Spearman=spearman,
    Prediction_status=status,
    Max_abs_K_test_train=info$max_cross,
    Prediction_range=pred_range,
    Prediction_SD=sd(pred),
    RMSE_raw=sqrt(mean(err^2)), MAE_raw=mean(abs(err)), Bias_raw=mean(err),
    RMSE_trainSD=sqrt(mean(err^2))/train_sd,
    MAE_trainSD=mean(abs(err))/train_sd, Bias_trainSD=mean(err)/train_sd,
    Vg=info$fit$Vu, Ve=info$fit$Ve, Lambda=info$fit$Lambda,
    stringsAsFactors=FALSE
  )
}

cat("\n===== Starting 5-fold matrix-baseline models =====\n")
for (model in models) {
  K <- Kmodels[[model]]
  for (trait in c("BB", "H6")) {
    for (outer in sprintf("%02d", 1:5)) {
      item <- fold_inputs[[paste(outer, trait, sep="|")]]
      out_fold <- file.path(outdir, "per_fold", trait, paste0("outer_", outer))
      dir.create(out_fold, recursive=TRUE, showWarnings=FALSE)
      metric_file <- file.path(out_fold, paste0(model, "_metrics.tsv"))
      pred_file <- file.path(out_fold, paste0(model, "_predictions.tsv"))
      flag_file <- file.path(out_fold, paste0(model, ".complete.flag"))
      if (!force && file.exists(flag_file) && file.exists(metric_file) && file.exists(pred_file)) {
        old <- read_flag(flag_file)
        same <- identical(old$status, "PASS") && identical(old$outer, outer) &&
          identical(old$trait, trait) && identical(old$model, model)
        if (same) {
          cat("[SKIP] ", model, " ", trait, " outer_", outer, " already complete\n", sep="")
          next
        }
        stop("Existing completion flag mismatches current task: ", flag_file,
             " (use --force only after inspecting old outputs)")
      }
      if (file.exists(flag_file)) unlink(flag_file)
      train_ids <- item$train_ids; test_ids <- item$test_ids; y_train <- item$y_train
      cat(sprintf("[RUN] %-15s trait=%s outer_%s train=%d test=%d\n",
                  model, trait, outer, length(train_ids), length(test_ids)))
      flush.console()
      # Cross-relationships determine whether an ABLUP model can distinguish
      # held-out individuals. Parent-connected CV may make every A[test, train]=0.
      max_cross <- max(abs(K[test_ids, train_ids, drop=FALSE]))
      if (model == "ABLUP") cat(sprintf(
        "[DIAG] ABLUP trait=%s outer_%s max_abs_A_test_train=%.12g\n",
        trait, outer, max_cross))
      fit <- fit_predict(K, train_ids, test_ids, y_train)
      pred <- fit$pred
      # Prediction vector is now complete. Only now open the evaluation-only test phenotype.
      test_pheno <- read_tsv(item$test_pheno_path)
      need(identical(names(test_pheno), c("ID", "phenotype")),
           paste("Unexpected evaluation-only test phenotype header:", item$test_pheno_path))
      need(nrow(test_pheno)==length(test_ids) && identical(as.character(test_pheno$ID), test_ids),
           paste("Evaluation-only test phenotype ID/order mismatch:", outer, trait))
      y_test <- suppressWarnings(as.numeric(test_pheno$phenotype))
      need(all(is.finite(y_test)) && all(y_test != -9),
           paste("Invalid evaluation-only test phenotype values:", outer, trait))
      met <- compute_metrics(y_test, pred, y_train,
        list(trait=trait, outer=outer, model=model, fit=fit, max_cross=max_cross))
      predictions <- data.frame(Trait=trait, Outer=as.integer(outer), Model=model,
        ID=test_ids, Observed_raw=y_test, Predicted_raw=pred,
        Residual_pred_minus_obs=pred-y_test, stringsAsFactors=FALSE)
      write_tsv_atomic(predictions, pred_file)
      write_tsv_atomic(met, metric_file)
      flag_tmp <- paste0(flag_file, ".tmp")
      writeLines(c("version=1", "status=PASS", paste0("outer=", outer),
        paste0("trait=", trait), paste0("model=", model),
        paste0("n_train=", length(train_ids)), paste0("n_test=", length(test_ids)),
        paste0("matrix_source=", if (model=="ABLUP") "pedigree_A" else if (model=="RRBLUP-array") basename(g_chip_path) else basename(g_imp_path)),
        "parent_connected_group_check=PASS", "test_phenotype_read_after_prediction=YES",
        paste0("prediction_status=", met$Prediction_status),
        paste0("max_abs_K_test_train=", format(max_cross, digits=16)),
        paste0("prediction_range=", format(met$Prediction_range, digits=16))), flag_tmp)
      if (file.exists(flag_file)) unlink(flag_file)
      need(file.rename(flag_tmp, flag_file), paste("Could not publish completion flag:", flag_file))
      cat(sprintf("[PASS] %-15s %s outer_%s Pearson=%s Spearman=%s RMSE=%.5f [%s]\n",
                  model, trait, outer,
                  if (is.na(met$Pearson)) "NA" else sprintf("%.5f", met$Pearson),
                  if (is.na(met$Spearman)) "NA" else sprintf("%.5f", met$Spearman),
                  met$RMSE_raw, met$Prediction_status))
      flush.console()
      rm(test_pheno, y_test, predictions, met, fit, pred); invisible(gc())
    }
  }
}

# Rebuild the aggregate outputs from all completed per-fold files.
all_metrics <- list(); joined_rows <- list(); k <- 0L
for (trait in c("BB", "H6")) for (outer in sprintf("%02d", 1:5)) for (model in models) {
  d <- file.path(outdir, "per_fold", trait, paste0("outer_", outer))
  mf <- file.path(d, paste0(model, "_metrics.tsv")); pf <- file.path(d, paste0(model, "_predictions.tsv"))
  ff <- file.path(d, paste0(model, ".complete.flag"))
  flag <- read_flag(ff)
  need(identical(flag$status, "PASS") && identical(flag$outer, outer) &&
       identical(flag$trait, trait) && identical(flag$model, model),
       paste("Missing/mismatched completion flag:", ff))
  met <- read_tsv(mf); pred <- read_tsv(pf)
  need(nrow(met)==1L && nrow(pred)==as.integer(flag$n_test), paste("Output row count mismatch:", d))
  need(all(pred$Trait==trait) && all(as.integer(pred$Outer)==as.integer(outer)) && all(pred$Model==model),
       paste("Prediction metadata mismatch:", pf))
  k <- k+1L; all_metrics[[k]] <- met
  joined_rows[[k]] <- data.frame(Trait=trait, Outer=as.integer(outer), Model=model,
    Pearson=as.numeric(met$Pearson), Spearman=as.numeric(met$Spearman),
    RMSE_original=as.numeric(met$RMSE_raw), MAE_original=as.numeric(met$MAE_raw),
    RMSE_trainSD=as.numeric(met$RMSE_trainSD), MAE_trainSD=as.numeric(met$MAE_trainSD),
    Prediction_status=as.character(met$Prediction_status),
    Source="Matrix-baselines", stringsAsFactors=FALSE)
}
all <- do.call(rbind, all_metrics)
need(nrow(all)==length(models)*2L*5L, "Unexpected aggregate metric count")
need(!anyDuplicated(paste(all$Trait, all$Outer, all$Model, sep="|")), "Duplicate aggregate model keys")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
write_tsv_atomic(all, file.path(outdir, "matrix_baselines_outer5_metrics.tsv"))
joined <- do.call(rbind, joined_rows)
write_tsv_atomic(joined, file.path(outdir, "seven_models_baseline_outer5.tsv"))
summary_rows <- list(); s <- 0L
for (trait in c("BB", "H6")) for (model in models) {
  x <- all[all$Trait==trait & all$Model==model, , drop=FALSE]
  need(nrow(x)==5L, paste("Expected five folds for", trait, model))
  z <- data.frame(Trait=trait, Model=model, N_outer_folds=nrow(x),
                  N_defined_Pearson=sum(is.finite(x$Pearson)), stringsAsFactors=FALSE)
  for (col in c("Pearson", "Spearman", "RMSE_raw", "MAE_raw", "Bias_raw",
                "RMSE_trainSD", "MAE_trainSD", "Bias_trainSD")) {
    valid <- x[[col]][is.finite(x[[col]])]
    z[[paste0(col, "_mean")]] <- if (length(valid)) mean(valid) else NA_real_
    z[[paste0(col, "_sd")]] <- if (length(valid)>1L) sd(valid) else NA_real_
  }
  s <- s+1L; summary_rows[[s]] <- z
}
summary <- do.call(rbind, summary_rows)
write_tsv_atomic(summary, file.path(outdir, "matrix_baselines_outer5_summary.tsv"))
writeLines(c("version=1", "status=PASS", paste0("n_models=", length(models)),
             paste0("n_metric_records=", nrow(all)), "n_outer_folds=5",
             "fold_design=parent_connected_group", "test_phenotype_read_after_prediction=YES",
             "baseline_labels=ABLUP,RRBLUP-array,RRBLUP-imputed",
             "undefined_ABLUP_correlation_allowed_only_for_zero_A_test_train=YES"),
           file.path(outdir, "complete.flag"))
cat("\n===== MATRIX BASELINES COMPLETE =====\n")
cat("Metric records: ", nrow(all), "\n", sep="")
cat("Seven-model join input: ", file.path(outdir, "seven_models_baseline_outer5.tsv"), "\n", sep="")
cat("Full metrics: ", file.path(outdir, "matrix_baselines_outer5_metrics.tsv"), "\n", sep="")
cat("Summary: ", file.path(outdir, "matrix_baselines_outer5_summary.tsv"), "\n", sep="")
print(summary[, c("Trait", "Model", "Pearson_mean", "Pearson_sd", "RMSE_trainSD_mean")], row.names=FALSE)
