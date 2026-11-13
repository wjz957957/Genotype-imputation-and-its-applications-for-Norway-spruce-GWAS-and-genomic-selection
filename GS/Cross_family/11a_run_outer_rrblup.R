#!/usr/bin/env Rscript
# Outer-fold RR-BLUP for the current 5-fold Cross_family workflow.
# Preserves the reference model: rrBLUP::mixed.solve(y=y_train, Z=G_train).
# The outer-test phenotype file is opened only after both model predictions
# have been generated and written to a prediction-only file.

suppressPackageStartupMessages(library(rrBLUP))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 11a_run_outer_rrblup.R WORK OUTER TRAIT")
}
WORK <- normalizePath(args[1], mustWork=TRUE)
OUTER <- sprintf("%02d", as.integer(args[2]))
TRAIT <- args[3]
if (!(OUTER %in% sprintf("%02d", 1:5)) || !(TRAIT %in% c("BB", "H6"))) {
  stop("Invalid OUTER/TRAIT; expected outer 01..05 and BB/H6")
}

if (requireNamespace("data.table", quietly=TRUE)) {
  nthreads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
  if (!is.finite(nthreads) || nthreads < 1L) nthreads <- 1L
  data.table::setDTthreads(nthreads)
}

BASE <- file.path(WORK, "Cross_family", "results", "10_outer_model_data",
                  paste0("outer_", OUTER), TRAIT)
OUT <- file.path(WORK, "Cross_family", "results", "11_outer_rrblup",
                 paste0("outer_", OUTER), TRAIT)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

need <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call.=FALSE)
}
read_flag <- function(path) {
  need(file.exists(path) && file.info(path)$size > 0, paste("Missing flag:", path))
  x <- readLines(path, warn=FALSE)
  out <- list()
  for (line in x) {
    p <- strsplit(line, "=", fixed=TRUE)[[1]]
    need(length(p) >= 2L, paste("Malformed flag line:", path, line))
    key <- p[1]
    value <- paste(p[-1], collapse="=")
    need(is.null(out[[key]]), paste("Duplicate flag key:", key, path))
    out[[key]] <- value
  }
  out
}
read_ids <- function(path) {
  need(file.exists(path) && file.info(path)$size > 0, paste("Missing ID list:", path))
  ids <- readLines(path, warn=FALSE)
  need(length(ids) > 0L && all(nzchar(ids)) && !anyDuplicated(ids),
       paste("Invalid/duplicate IDs:", path))
  ids
}
read_raw_header <- function(path) {
  need(file.exists(path) && file.info(path)$size > 0, paste("Missing RAW:", path))
  line <- readLines(path, n=1L, warn=FALSE)
  h <- strsplit(trimws(line), "[[:space:]]+")[[1]]
  need(length(h) >= 7L &&
       identical(h[1:6], c("FID","IID","PAT","MAT","SEX","PHENOTYPE")),
       paste("Invalid PLINK RAW header:", path))
  snps <- sub("_[^_]+$", "", h[7:length(h)])
  need(all(nzchar(snps)) && !anyDuplicated(snps),
       paste("Empty/duplicate SNP IDs in RAW header:", path))
  list(header=h, snps=snps)
}
read_raw_matrix <- function(path, expected_ids, expected_header) {
  if (requireNamespace("data.table", quietly=TRUE)) {
    d <- data.table::fread(path, header=TRUE, data.table=FALSE,
                           check.names=FALSE, showProgress=FALSE,
                           na.strings=c("NA","NaN","."))
  } else {
    d <- read.table(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE,
                    comment.char="", quote="", na.strings=c("NA","NaN","."))
  }
  need(identical(names(d), expected_header), paste("RAW header changed while reading:", path))
  need(nrow(d) == length(expected_ids), paste("RAW row count mismatch:", path))
  ids <- as.character(d[["IID"]])
  fids <- as.character(d[["FID"]])
  need(identical(ids, expected_ids) && identical(fids, expected_ids),
       paste("RAW sample IDs/order mismatch:", path))
  pheno_meta <- suppressWarnings(as.numeric(d[["PHENOTYPE"]]))
  need(all(is.finite(pheno_meta)) && all(pheno_meta == -9),
       paste("RAW PHENOTYPE metadata should be -9:", path))
  G <- as.matrix(d[, 7:ncol(d), drop=FALSE])
  storage.mode(G) <- "numeric"
  need(all(is.finite(G)), paste("Missing/non-finite genotype values:", path))
  need(all(G %in% c(0,1,2)), paste("Genotype values outside 0/1/2:", path))
  G
}
read_train_pheno <- function(path, expected_ids) {
  need(file.exists(path) && file.info(path)$size > 0, paste("Missing training phenotype:", path))
  p <- read.delim(path, header=TRUE, sep="\t", check.names=FALSE,
                  stringsAsFactors=FALSE)
  need(identical(names(p), c("FID","IID","phenotype")), paste("Unexpected train phenotype header:", path))
  need(nrow(p) == length(expected_ids), paste("Training phenotype count mismatch:", path))
  ids <- as.character(p$IID)
  need(identical(as.character(p$FID), expected_ids) && identical(ids, expected_ids),
       paste("Training phenotype IDs/order mismatch:", path))
  y <- suppressWarnings(as.numeric(p$phenotype))
  need(all(is.finite(y)) && all(y != -9), paste("Missing/non-finite training phenotype:", path))
  y
}
read_test_pheno <- function(path, expected_ids) {
  need(file.exists(path) && file.info(path)$size > 0, paste("Missing evaluation-only test phenotype:", path))
  p <- read.delim(path, header=TRUE, sep="\t", check.names=FALSE,
                  stringsAsFactors=FALSE)
  need(identical(names(p), c("ID","phenotype")), paste("Unexpected test phenotype header:", path))
  need(nrow(p) == length(expected_ids), paste("Test phenotype count mismatch:", path))
  ids <- as.character(p$ID)
  need(identical(ids, expected_ids), paste("Test phenotype IDs/order mismatch:", path))
  y <- suppressWarnings(as.numeric(p$phenotype))
  need(all(is.finite(y)) && all(y != -9), paste("Missing/non-finite test phenotype:", path))
  y
}
read_tsv <- function(path) {
  need(file.exists(path) && file.info(path)$size > 0, paste("Missing TSV:", path))
  read.delim(path, header=TRUE, sep="\t", check.names=FALSE,
             stringsAsFactors=FALSE)
}
write_tsv <- function(x, path) {
  write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
safe_cor <- function(x, y, method="pearson") {
  if (length(x) < 2L || !is.finite(sd(x)) || !is.finite(sd(y)) ||
      sd(x) == 0 || sd(y) == 0) return(NA_real_)
  as.numeric(cor(x, y, method=method))
}
fit_predict <- function(G_train, G_test, y_train) {
  fit <- rrBLUP::mixed.solve(y=y_train, Z=G_train, SE=FALSE, return.Hinv=FALSE)
  pred <- as.numeric(G_test %*% fit$u + as.numeric(fit$beta))
  need(length(pred) == nrow(G_test) && all(is.finite(pred)),
       "Invalid RRBLUP predictions")
  list(pred=pred, Vu=as.numeric(fit$Vu), Ve=as.numeric(fit$Ve))
}
metric_one <- function(y, pred, train_sd) {
  err <- pred - y
  c(Pearson=safe_cor(y, pred, "pearson"),
    Spearman=safe_cor(y, pred, "spearman"),
    RMSE_raw=sqrt(mean(err^2)),
    RMSE_trainSD=sqrt(mean(err^2))/train_sd,
    MAE_raw=mean(abs(err)),
    Bias_raw=mean(err))
}

flag_path <- file.path(OUT, "complete.flag")
manifest_path <- file.path(BASE, "final_genotype_manifest.tsv")
need(file.exists(file.path(BASE, "complete.flag")) && file.exists(manifest_path),
     paste("Incomplete model-data input directory:", BASE))
input_flag <- read_flag(file.path(BASE, "complete.flag"))
need(identical(input_flag[["status"]], "PASS") &&
     identical(input_flag[["outer"]], OUTER) &&
     identical(input_flag[["trait"]], TRAIT),
     "Input complete.flag does not match this task")

manifest <- read_tsv(manifest_path)
need(nrow(manifest) == 1L, "Expected one row in final_genotype_manifest.tsv")
m <- manifest[1, , drop=FALSE]
# read.delim() parses a manifest value such as "01" as integer 1.
# Compare normalized numeric fold IDs so leading-zero formatting cannot cause
# a false mismatch; keep the separate complete.flag check for exact task identity.
manifest_outer <- suppressWarnings(as.integer(as.character(m$Outer)))
task_outer <- suppressWarnings(as.integer(OUTER))
need(length(manifest_outer) == 1L && !is.na(manifest_outer) &&
     length(task_outer) == 1L && !is.na(task_outer) &&
     manifest_outer == task_outer && as.character(m$Trait) == TRAIT,
     "Input manifest fold/trait mismatch")
N_CHIP <- as.integer(m$N_chip)
N_QTN <- as.integer(m$N_QTN)
N_TRAIN <- as.integer(m$N_train)
N_TEST <- as.integer(m$N_test)
need(all(is.finite(c(N_CHIP,N_QTN,N_TRAIN,N_TEST))) &&
     N_CHIP > 0L && N_QTN >= 0L && N_TRAIN > 0L && N_TEST > 0L,
     "Invalid counts in model-data manifest")

if (file.exists(flag_path) && file.info(flag_path)$size > 0) {
  old <- read_flag(flag_path)
  same <- identical(old[["status"]], "PASS") &&
    identical(old[["outer"]], OUTER) && identical(old[["trait"]], TRAIT) &&
    identical(old[["N_chip"]], as.character(N_CHIP)) &&
    identical(old[["N_qtn"]], as.character(N_QTN)) &&
    file.exists(file.path(OUT,"outer_metrics.tsv")) &&
    file.exists(file.path(OUT,"outer_test_predictions.tsv"))
  if (same) {
    cat(sprintf("[SKIP] already completed outer_%s/%s\n", OUTER, TRAIT))
    quit(save="no", status=0L)
  }
  stop("Existing output complete.flag mismatches current inputs; inspect output before rerunning")
}
unlink(flag_path)

train_ids <- read_ids(file.path(BASE, "train.ids"))
test_ids <- read_ids(file.path(BASE, "test.ids"))
need(length(train_ids) == N_TRAIN && length(test_ids) == N_TEST,
     "ID counts differ from model-data manifest")
need(length(intersect(train_ids,test_ids)) == 0L, "Train/test IDs overlap")

y_train <- read_train_pheno(file.path(BASE, "train.pheno.tsv"), train_ids)
need(is.finite(sd(y_train)) && sd(y_train) > 0, "Invalid training phenotype SD")

chip_train_path <- file.path(BASE, "chip_train.raw")
chip_test_path <- file.path(BASE, "chip_test.raw")
combined_train_path <- file.path(BASE, "combined_train.raw")
combined_test_path <- file.path(BASE, "combined_test.raw")
for (p in c(chip_train_path,chip_test_path,combined_train_path,combined_test_path)) {
  need(file.exists(p) && file.info(p)$size > 0, paste("Missing model matrix:", p))
}
ch_tr <- read_raw_header(chip_train_path)
ch_te <- read_raw_header(chip_test_path)
co_tr <- read_raw_header(combined_train_path)
co_te <- read_raw_header(combined_test_path)
need(identical(ch_tr$header,ch_te$header), "CHIP train/test headers differ")
need(identical(co_tr$header,co_te$header), "Combined train/test headers differ")
need(length(ch_tr$snps) == N_CHIP, "CHIP SNP count differs from manifest")
need(length(co_tr$snps) == N_CHIP + N_QTN, "Combined SNP count differs from manifest")
need(identical(co_tr$snps[seq_len(N_CHIP)], ch_tr$snps),
     "Combined matrix does not start with CHIP SNPs in the same order")
need(!anyDuplicated(co_tr$snps), "Duplicate SNPs in combined matrix")
chip_selected <- read_ids(file.path(BASE, "selected_chip.snps"))
qtn_selected <- if (N_QTN > 0L) read_ids(file.path(BASE, "selected_qtn.snps")) else character(0)
need(length(chip_selected) == N_CHIP && setequal(chip_selected,ch_tr$snps),
     "Selected CHIP list differs from CHIP matrix SNP set")
need(length(qtn_selected) == N_QTN &&
     setequal(qtn_selected,co_tr$snps[-seq_len(N_CHIP)]),
     "Selected QTN list differs from Combined matrix QTN set")
need(length(intersect(chip_selected,qtn_selected)) == 0L, "CHIP/QTN SNP overlap")

cat(sprintf("[RUN] outer_%s/%s: CHIP=%d QTN=%d train=%d test=%d\n",
            OUTER, TRAIT, N_CHIP, N_QTN, N_TRAIN, N_TEST), flush=TRUE)

# Fit the CHIP baseline first, then release its matrices before loading the larger combined matrix.
G_chip_train <- read_raw_matrix(chip_train_path, train_ids, ch_tr$header)
G_chip_test <- read_raw_matrix(chip_test_path, test_ids, ch_te$header)
chip_fit <- fit_predict(G_chip_train, G_chip_test, y_train)
pred_chip <- chip_fit$pred
rm(G_chip_train, G_chip_test)
gc()

if (N_QTN == 0L) {
  need(identical(co_tr$header,ch_tr$header), "Zero-QTN combined header differs from CHIP")
  pred_combined <- pred_chip
  combined_fit <- chip_fit
} else {
  G_combined_train <- read_raw_matrix(combined_train_path, train_ids, co_tr$header)
  G_combined_test <- read_raw_matrix(combined_test_path, test_ids, co_te$header)
  combined_fit <- fit_predict(G_combined_train, G_combined_test, y_train)
  pred_combined <- combined_fit$pred
  rm(G_combined_train, G_combined_test)
  gc()
}

# Freeze model predictions before reading the evaluation-only phenotype file.
pred_only <- data.frame(Outer=OUTER, Trait=TRAIT, IID=test_ids,
                        RRBLUP_CHIP_raw=pred_chip,
                        RRBLUP_COMBINED_raw=pred_combined,
                        stringsAsFactors=FALSE)
write_tsv(pred_only, file.path(OUT, "outer_test_predictions.predictions_only.tsv"))

test_pheno_path <- file.path(BASE, "evaluation_only", "test.pheno.tsv")
y_test <- read_test_pheno(test_pheno_path, test_ids)
m_chip <- metric_one(y_test, pred_chip, sd(y_train))
m_combined <- metric_one(y_test, pred_combined, sd(y_train))

pred_final <- data.frame(Outer=OUTER, Trait=TRAIT, IID=test_ids,
                         Observed_raw=y_test,
                         RRBLUP_CHIP_raw=pred_chip,
                         RRBLUP_COMBINED_raw=pred_combined,
                         stringsAsFactors=FALSE)
write_tsv(pred_final, file.path(OUT, "outer_test_predictions.tsv"))

met <- data.frame(
  Outer=OUTER, Trait=TRAIT,
  Best_Chip=as.character(m$Best_Chip),
  N_chip_selected=N_CHIP, N_QTN_selected=N_QTN,
  N_total_selected=N_CHIP+N_QTN,
  N_train=N_TRAIN, N_test=N_TEST,
  Pearson_CHIP=m_chip[["Pearson"]],
  Pearson_COMBINED=m_combined[["Pearson"]],
  Spearman_CHIP=m_chip[["Spearman"]],
  Spearman_COMBINED=m_combined[["Spearman"]],
  RMSE_raw_CHIP=m_chip[["RMSE_raw"]],
  RMSE_raw_COMBINED=m_combined[["RMSE_raw"]],
  RMSE_trainSD_CHIP=m_chip[["RMSE_trainSD"]],
  RMSE_trainSD_COMBINED=m_combined[["RMSE_trainSD"]],
  MAE_raw_CHIP=m_chip[["MAE_raw"]],
  MAE_raw_COMBINED=m_combined[["MAE_raw"]],
  Bias_raw_CHIP=m_chip[["Bias_raw"]],
  Bias_raw_COMBINED=m_combined[["Bias_raw"]],
  Delta_Pearson=m_combined[["Pearson"]]-m_chip[["Pearson"]],
  Delta_RMSE_raw=m_combined[["RMSE_raw"]]-m_chip[["RMSE_raw"]],
  Vu_CHIP=chip_fit$Vu, Ve_CHIP=chip_fit$Ve,
  Vu_COMBINED=combined_fit$Vu, Ve_COMBINED=combined_fit$Ve,
  stringsAsFactors=FALSE
)
write_tsv(met, file.path(OUT, "outer_metrics.tsv"))

model_manifest <- data.frame(
  Outer=OUTER, Trait=TRAIT,
  Model=c("CHIP","CHIP_PLUS_QTN"),
  Method="rrBLUP::mixed.solve",
  N_chip=N_CHIP, N_QTN=c(0L,N_QTN),
  N_total=c(N_CHIP,N_CHIP+N_QTN),
  N_train=N_TRAIN, N_test=N_TEST,
  Genotype_encoding="raw 0/1/2 dosage; no marker standardization",
  Test_phenotype_used_for_fitting="NO",
  stringsAsFactors=FALSE
)
write_tsv(model_manifest, file.path(OUT, "model_manifest.tsv"))
writeLines(c(
  "status=PASS", "version=1", paste0("outer=",OUTER), paste0("trait=",TRAIT),
  paste0("N_chip=",N_CHIP), paste0("N_qtn=",N_QTN),
  paste0("N_train=",N_TRAIN), paste0("N_test=",N_TEST),
  "method=rrBLUP::mixed.solve", "genotype_encoding=0/1/2",
  "test_phenotype_used_for_fitting=NO",
  "test_phenotype_loaded_after_prediction_generation=YES"
), flag_path)

cat(sprintf("[PASS] outer_%s/%s: CHIP Pearson=%s; combined Pearson=%s; delta=%s\n",
            OUTER, TRAIT,
            format(met$Pearson_CHIP, digits=8),
            format(met$Pearson_COMBINED, digits=8),
            format(met$Delta_Pearson, digits=8)))
