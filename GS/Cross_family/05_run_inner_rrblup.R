#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 05_run_inner_rrblup.R WORK OUTER INNER")
}

WORK <- normalizePath(args[1], mustWork = TRUE)
OUTER_NUM <- suppressWarnings(as.integer(args[2]))
INNER_NUM <- suppressWarnings(as.integer(args[3]))
if (is.na(OUTER_NUM) || OUTER_NUM < 1L || OUTER_NUM > 5L ||
    is.na(INNER_NUM) || INNER_NUM < 1L || INNER_NUM > 3L) {
  stop("OUTER must be 1..5 and INNER must be 1..3")
}
OUTER <- sprintf("%02d", OUTER_NUM)
INNER <- sprintf("%02d", INNER_NUM)

N_EXPECTED_SNP <- 43267L
SCRIPT_VERSION <- 1L
TRAITS <- c("BB", "H6")
SIZES <- c(50L, 100L, 200L, 300L, 400L, 500L,
           1000L, 2000L, 3000L, 4000L, 5000L, 8000L,
           10000L, 20000L, 30000L, 40000L)
MARKER_SETS <- c(paste0("chip_", SIZES), "chip_LD_ALL", "chip_ALL")

if (!requireNamespace("rrBLUP", quietly = TRUE)) {
  stop("Package rrBLUP is unavailable in this R environment")
}

GENO_DIR <- file.path(WORK, "Cross_family", "results",
                      "04_inner_rrblup_genotypes",
                      paste0("outer_", OUTER), paste0("inner_", INNER))
PHENO_DIR <- file.path(WORK, "Cross_family", "results",
                       "02_inner_phenotypes",
                       paste0("outer_", OUTER), paste0("inner_", INNER))
MARKER_DIR <- file.path(WORK, "Cross_family", "results",
                        "03_inner_chip_gwas_clumping",
                        paste0("outer_", OUTER), paste0("inner_", INNER),
                        "marker_sets")
OUT_DIR <- file.path(WORK, "Cross_family", "results",
                     "05_inner_rrblup",
                     paste0("outer_", OUTER), paste0("inner_", INNER))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FLAG <- file.path(OUT_DIR, "complete.flag")
METRICS_FILE <- file.path(OUT_DIR, "metrics.tsv")
PRED_FILE <- file.path(OUT_DIR, "validation_predictions.tsv")

cat("========================================\n")
cat("Inner-fold RRBLUP\n")
cat("Outer:", OUTER, " Inner:", INNER, "\n")
cat("R version:", as.character(getRversion()), "\n")
cat("rrBLUP version:", as.character(packageVersion("rrBLUP")), "\n")
cat("========================================\n")

require_nonempty <- function(path, label = path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Missing or empty ", label, ": ", path)
  }
  invisible(path)
}

# Confirm that the upstream genotype preparation completed and passed its key checks.
GENO_FLAG <- file.path(GENO_DIR, "complete.flag")
require_nonempty(GENO_FLAG, "genotype completion flag")
flag_lines <- readLines(GENO_FLAG, warn = FALSE)
required_flags <- c("version=2",
                    "train_valid_snp_header_match=PASS",
                    "train_valid_snp_order_match=PASS",
                    "outer_test_excluded=PASS")
if (!all(required_flags %in% flag_lines)) {
  stop("Genotype complete.flag is missing one or more required PASS entries: ", GENO_FLAG)
}

# Read PLINK --recode A output while preserving the original SNP column names.
read_raw <- function(path, label) {
  require_nonempty(path, label)
  d <- read.table(path,
                  header = TRUE,
                  sep = "",
                  check.names = FALSE,
                  stringsAsFactors = FALSE,
                  na.strings = c("NA", "NaN"),
                  comment.char = "",
                  quote = "",
                  fill = FALSE)

  if (ncol(d) != 6L + N_EXPECTED_SNP) {
    stop(label, ": expected ", 6L + N_EXPECTED_SNP,
         " columns; found ", ncol(d))
  }
  if (!identical(toupper(names(d)[1:2]), c("FID", "IID"))) {
    stop(label, ": first two columns must be FID and IID")
  }
  d$FID <- as.character(d$FID)
  d$IID <- as.character(d$IID)
  if (anyNA(d$FID) || anyNA(d$IID) || anyDuplicated(d$IID)) {
    stop(label, ": missing or duplicated sample IDs")
  }
  if (any(d$FID != d$IID)) {
    stop(label, ": FID/IID mismatch")
  }
  d
}

RAW_TRAIN <- read_raw(file.path(GENO_DIR, "chip_all_train.raw"), "Training raw")
RAW_VALID <- read_raw(file.path(GENO_DIR, "chip_all_valid.raw"), "Validation raw")

TRAIN_HEADERS_FILE <- file.path(GENO_DIR, "train_snp_headers.txt")
VALID_HEADERS_FILE <- file.path(GENO_DIR, "valid_snp_headers.txt")
require_nonempty(TRAIN_HEADERS_FILE, "training SNP header list")
require_nonempty(VALID_HEADERS_FILE, "validation SNP header list")
TRAIN_HEADERS <- readLines(TRAIN_HEADERS_FILE, warn = FALSE)
VALID_HEADERS <- readLines(VALID_HEADERS_FILE, warn = FALSE)
RAW_HEADERS_TRAIN <- names(RAW_TRAIN)[7:ncol(RAW_TRAIN)]
RAW_HEADERS_VALID <- names(RAW_VALID)[7:ncol(RAW_VALID)]
if (!identical(RAW_HEADERS_TRAIN, RAW_HEADERS_VALID)) {
  stop("Training and validation SNP headers differ")
}
if (!identical(TRAIN_HEADERS, RAW_HEADERS_TRAIN) ||
    !identical(VALID_HEADERS, RAW_HEADERS_VALID)) {
  stop("Saved SNP header lists do not match the .raw headers")
}
RAW_SNP_IDS <- sub("_[^_]+$", "", RAW_HEADERS_TRAIN)
if (length(RAW_SNP_IDS) != N_EXPECTED_SNP || anyDuplicated(RAW_SNP_IDS)) {
  stop("Raw genotype SNP IDs are invalid or duplicated")
}

# Phenotype files are headerless three-column files: FID, IID, phenotype.
read_pheno <- function(path) {
  require_nonempty(path, "phenotype file")
  d <- read.table(path,
                  header = FALSE,
                  sep = "",
                  colClasses = "character",
                  stringsAsFactors = FALSE,
                  na.strings = c("NA", "NaN"),
                  comment.char = "",
                  quote = "",
                  fill = FALSE)
  if (ncol(d) != 3L) {
    stop("Phenotype file must have exactly 3 columns (FID IID y): ", path)
  }
  out <- data.frame(FID = as.character(d[[1]]),
                    IID = as.character(d[[2]]),
                    y = suppressWarnings(as.numeric(d[[3]])),
                    stringsAsFactors = FALSE)
  if (nrow(out) == 0L || anyNA(out$FID) || anyNA(out$IID) ||
      anyDuplicated(out$IID)) {
    stop("Invalid or duplicated phenotype IDs: ", path)
  }
  if (any(out$FID != out$IID)) {
    stop("Phenotype FID/IID mismatch: ", path)
  }
  if (any(!is.finite(out$y))) {
    stop("Missing or non-finite phenotype values: ", path)
  }
  out
}

# Return genotype row indices in the same order as the phenotype records.
match_samples <- function(raw, pheno, label) {
  if (nrow(raw) != nrow(pheno) || !setequal(raw$IID, pheno$IID)) {
    stop(label, ": phenotype and genotype sample IDs differ")
  }
  idx <- match(pheno$IID, raw$IID)
  if (anyNA(idx) || any(raw$FID[idx] != pheno$FID)) {
    stop(label, ": unable to match phenotype IDs and FIDs to genotypes")
  }
  idx
}

# Remove any stale completion marker before beginning a new run.
if (file.exists(FLAG)) unlink(FLAG)
METRIC_ROWS <- list()
PRED_ROWS <- list()
K <- 0L

for (TRAIT in TRAITS) {
  PHENO_TRAIN <- read_pheno(file.path(PHENO_DIR, paste0(TRAIT, ".train.pheno")))
  PHENO_VALID <- read_pheno(file.path(PHENO_DIR, paste0(TRAIT, ".valid.pheno")))
  TRAIN_IDX <- match_samples(RAW_TRAIN, PHENO_TRAIN, paste(TRAIT, "training"))
  VALID_IDX <- match_samples(RAW_VALID, PHENO_VALID, paste(TRAIT, "validation"))
  Y_TRAIN <- PHENO_TRAIN$y
  Y_VALID <- PHENO_VALID$y

  for (SET in MARKER_SETS) {
    MARKER_FILE <- file.path(MARKER_DIR, paste0(TRAIT, "_", SET, ".snps"))
    require_nonempty(MARKER_FILE, paste(TRAIT, SET, "marker set"))
    MARKERS <- trimws(readLines(MARKER_FILE, warn = FALSE))
    MARKERS <- MARKERS[nzchar(MARKERS)]

    if (length(MARKERS) == 0L || anyDuplicated(MARKERS)) {
      stop("Empty or duplicated marker set: ", MARKER_FILE)
    }
    if (grepl("^chip_[0-9]+$", SET)) {
      expected_n <- as.integer(sub("^chip_", "", SET))
      if (length(MARKERS) != expected_n) {
        stop(SET, ": expected ", expected_n, " markers; found ", length(MARKERS))
      }
    }
    if (SET == "chip_ALL" && length(MARKERS) != N_EXPECTED_SNP) {
      stop("chip_ALL does not contain all expected chip SNPs")
    }
    COL_IDX <- match(MARKERS, RAW_SNP_IDS)
    if (anyNA(COL_IDX)) {
      absent <- MARKERS[is.na(COL_IDX)]
      stop(TRAIT, "/", SET, ": ", length(absent),
           " markers absent from raw genotypes; first: ", absent[1])
    }
    COLS <- COL_IDX + 6L

    X_TRAIN <- as.matrix(RAW_TRAIN[TRAIN_IDX, COLS, drop = FALSE])
    X_VALID <- as.matrix(RAW_VALID[VALID_IDX, COLS, drop = FALSE])
    storage.mode(X_TRAIN) <- "double"
    storage.mode(X_VALID) <- "double"
    X_TRAIN[!is.finite(X_TRAIN)] <- NA_real_
    X_VALID[!is.finite(X_VALID)] <- NA_real_

    # Impute using training-set means only.
    MEANS <- colMeans(X_TRAIN, na.rm = TRUE)
    KEEP_MEAN <- is.finite(MEANS)
    if (!all(KEEP_MEAN)) {
      X_TRAIN <- X_TRAIN[, KEEP_MEAN, drop = FALSE]
      X_VALID <- X_VALID[, KEEP_MEAN, drop = FALSE]
      MEANS <- MEANS[KEEP_MEAN]
    }
    if (ncol(X_TRAIN) == 0L) {
      stop(TRAIT, "/", SET, ": no usable markers after missingness filtering")
    }
    for (j in seq_len(ncol(X_TRAIN))) {
      miss_train <- !is.finite(X_TRAIN[, j])
      miss_valid <- !is.finite(X_VALID[, j])
      if (any(miss_train)) X_TRAIN[miss_train, j] <- MEANS[j]
      if (any(miss_valid)) X_VALID[miss_valid, j] <- MEANS[j]
    }

    # Remove markers that are monomorphic in the training set only.
    SD_TRAIN <- apply(X_TRAIN, 2L, sd)
    KEEP_VAR <- is.finite(SD_TRAIN) & SD_TRAIN > 0
    X_TRAIN <- X_TRAIN[, KEEP_VAR, drop = FALSE]
    X_VALID <- X_VALID[, KEEP_VAR, drop = FALSE]
    if (ncol(X_TRAIN) == 0L) {
      stop(TRAIT, "/", SET, ": all markers are monomorphic in training")
    }

    # Center both partitions using training-set allele-dosage means.
    CENTER <- colMeans(X_TRAIN)
    X_TRAIN <- sweep(X_TRAIN, 2L, CENTER, FUN = "-")
    X_VALID <- sweep(X_VALID, 2L, CENTER, FUN = "-")

    cat("Fitting:", TRAIT, SET,
        "| requested SNPs:", length(MARKERS),
        "| usable SNPs:", ncol(X_TRAIN), "\n")

    X_FIXED <- matrix(1, nrow = length(Y_TRAIN), ncol = 1L)
    FIT <- tryCatch(
      rrBLUP::mixed.solve(y = Y_TRAIN, Z = X_TRAIN, X = X_FIXED),
      error = function(e) {
        stop("RRBLUP failed for ", TRAIT, "/", SET, ": ", conditionMessage(e))
      }
    )
    BETA <- as.numeric(FIT$beta[1])
    U <- as.numeric(FIT$u)
    PRED <- as.numeric(BETA + X_VALID %*% U)
    if (length(PRED) != length(Y_VALID) || any(!is.finite(PRED))) {
      stop(TRAIT, "/", SET, ": invalid predictions")
    }

    PEARSON_R <- if (length(Y_VALID) > 1L && sd(Y_VALID) > 0 && sd(PRED) > 0) {
      as.numeric(cor(Y_VALID, PRED, method = "pearson"))
    } else {
      NA_real_
    }
    RMSE <- sqrt(mean((PRED - Y_VALID)^2))
    MAE <- mean(abs(PRED - Y_VALID))
    BIAS <- mean(PRED - Y_VALID)
    VU <- as.numeric(FIT$Vu)
    VE <- as.numeric(FIT$Ve)
    LAMBDA <- if (is.finite(VU) && VU > 0) VE / VU else NA_real_

    K <- K + 1L
    METRIC_ROWS[[K]] <- data.frame(
      outer = OUTER, inner = INNER, trait = TRAIT, marker_set = SET,
      n_train = length(Y_TRAIN), n_valid = length(Y_VALID),
      n_markers_requested = length(MARKERS), n_markers_used = ncol(X_TRAIN),
      pearson_r = PEARSON_R, RMSE = RMSE, MAE = MAE, bias = BIAS,
      Vu = VU, Ve = VE, lambda = LAMBDA, status = "PASS",
      stringsAsFactors = FALSE
    )
    PRED_ROWS[[K]] <- data.frame(
      outer = OUTER, inner = INNER, trait = TRAIT, marker_set = SET,
      sample_id = PHENO_VALID$IID, observed = Y_VALID, predicted = PRED,
      stringsAsFactors = FALSE
    )
    cat(sprintf("PASS %s/%s | r=%s | RMSE=%.4f | MAE=%.4f\n",
                TRAIT, SET,
                if (is.finite(PEARSON_R)) sprintf("%.4f", PEARSON_R) else "NA",
                RMSE, MAE))
    rm(X_TRAIN, X_VALID, FIT, U, PRED)
    invisible(gc())
  }
}

METRICS <- do.call(rbind, METRIC_ROWS)
PREDICTIONS <- do.call(rbind, PRED_ROWS)
if (nrow(METRICS) != 36L || length(unique(METRICS$marker_set)) != 18L ||
    any(METRICS$status != "PASS") || nrow(PREDICTIONS) == 0L) {
  stop("Final model count/status validation failed")
}
if (anyDuplicated(paste(METRICS$trait, METRICS$marker_set, sep = "::"))) {
  stop("Duplicate trait/marker-set metric rows")
}

METRICS_TMP <- paste0(METRICS_FILE, ".tmp")
PRED_TMP <- paste0(PRED_FILE, ".tmp")
write.table(METRICS, file = METRICS_TMP, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(PREDICTIONS, file = PRED_TMP, sep = "\t", quote = FALSE, row.names = FALSE)
if (file.exists(METRICS_FILE)) unlink(METRICS_FILE)
if (file.exists(PRED_FILE)) unlink(PRED_FILE)
if (!file.rename(METRICS_TMP, METRICS_FILE)) stop("Could not publish metrics.tsv")
if (!file.rename(PRED_TMP, PRED_FILE)) stop("Could not publish validation_predictions.tsv")

FLAG_TMP <- paste0(FLAG, ".tmp")
writeLines(c(
  paste0("version=", SCRIPT_VERSION),
  paste0("outer=", OUTER),
  paste0("inner=", INNER),
  "traits=BB,H6",
  "marker_sets_per_trait=18",
  "models=36",
  "model_status=PASS",
  "sample_ID_check=PASS",
  "marker_ID_check=PASS",
  "training_only_imputation_and_centering=PASS",
  "validation_prediction_check=PASS"
), con = FLAG_TMP)
if (!file.rename(FLAG_TMP, FLAG)) stop("Could not publish completion flag")

cat("========================================\n")
cat("OVERALL PASS\n")
cat("Outer:", OUTER, " Inner:", INNER, "\n")
cat("Models fitted:", nrow(METRICS), "\n")
cat("Metrics:", METRICS_FILE, "\n")
cat("Predictions:", PRED_FILE, "\n")
cat("========================================\n")
