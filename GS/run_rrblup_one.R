args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript run_rrblup_one.R <BB|H6> <chip|extra|combined> <01-10>")
}

trait <- args[1]
geno_type <- args[2]
outer <- sprintf("%02d", as.integer(args[3]))

if (!trait %in% c("BB", "H6"))
    stop("trait must be BB or H6")

if (!geno_type %in% c("chip", "extra", "combined"))
    stop("geno_type must be chip, extra, or combined")

cat("===== RR-BLUP =====\n")
cat("trait:", trait, "\n")
cat("genotype:", geno_type, "\n")
cat("outer:", outer, "\n")


# ============================================================
# LOAD PACKAGE
# ============================================================

suppressPackageStartupMessages(
    library(rrBLUP)
)


# ============================================================
# FILE PATHS
# ============================================================

train_file <- file.path(
    "final_genotypes",
    paste0("outer_", outer),
    trait,
    paste0(geno_type, "_train.raw")
)

test_file <- file.path(
    "final_genotypes",
    paste0("outer_", outer),
    trait,
    paste0(geno_type, "_test.raw")
)

cv_train_file <- file.path(
    "CV",
    paste0("outer_", outer),
    "outer_train.samples"
)

cv_test_file <- file.path(
    "CV",
    paste0("outer_", outer),
    "outer_test.samples"
)

pheno_file <- file.path(
    "phenotype",
    paste0("Trait-Adj_", trait, "-sort-Z-stander.txt")
)


# ============================================================
# CHECK FILES
# ============================================================

if (!file.exists(train_file))
    stop("Missing: ", train_file)

if (!file.exists(test_file))
    stop("Missing: ", test_file)

if (!file.exists(cv_train_file))
    stop("Missing: ", cv_train_file)

if (!file.exists(cv_test_file))
    stop("Missing: ", cv_test_file)

if (!file.exists(pheno_file))
    stop("Missing: ", pheno_file)


# ============================================================
# READ GENOTYPES
# ============================================================

cat("\n===== READ GENOTYPES =====\n")

read_genotype <- function(path) {

    d <- read.table(
        path,
        header = TRUE,
        sep = "",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (ncol(d) < 7)
        stop("Expected at least 7 columns in: ", path)

    # IID is column 2
    ids <- as.character(d[[2]])

    # Marker columns start from column 7
    markers <- colnames(d)[7:ncol(d)]

    # Check duplicated sample IDs
    if (anyDuplicated(ids))
        stop("Duplicated sample IDs in: ", path)

    # Check duplicated marker names
    if (anyDuplicated(markers))
        stop("Duplicated marker names in: ", path)

    # Convert genotype matrix
    G <- data.matrix(
        d[, 7:ncol(d), drop = FALSE]
    )

    # Check missing values
    if (anyNA(G))
        stop("NA found in genotype matrix: ", path)

    # Check genotype values
    if (any(!(G %in% c(0, 1, 2, -9))))
        stop("Unexpected genotype value found in: ", path)

    list(
        ids = ids,
        markers = markers,
        G = G
    )
}


train <- read_genotype(train_file)
test  <- read_genotype(test_file)

cat("train samples :", nrow(train$G), "\n")
cat("test samples  :", nrow(test$G), "\n")
cat("markers       :", ncol(train$G), "\n")


# ============================================================
# CHECK TRAIN / TEST MARKER CONSISTENCY
# ============================================================

if (!identical(train$markers, test$markers))
    stop("Train/test marker columns are not identical")

if (ncol(train$G) != ncol(test$G))
    stop("Train/test marker number differs")


# ============================================================
# SAMPLE ORDER CHECK
# ============================================================

cat("\n===== SAMPLE ORDER CHECK =====\n")

cv_train_ids <- scan(
    cv_train_file,
    what = "",
    quiet = TRUE
)

cv_test_ids <- scan(
    cv_test_file,
    what = "",
    quiet = TRUE
)


if (!identical(train$ids, cv_train_ids))
    stop("Training sample order does not match CV/outer_train.samples")

if (!identical(test$ids, cv_test_ids))
    stop("Testing sample order does not match CV/outer_test.samples")


cat("train order: PASS\n")
cat("test order : PASS\n")


# ============================================================
# READ PHENOTYPE
# ============================================================

cat("\n===== READ PHENOTYPE =====\n")

pheno <- read.table(
    pheno_file,
    header = TRUE,
    sep = "",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (ncol(pheno) < 2)
    stop("Phenotype file must contain at least 2 columns")

pheno_ids <- as.character(pheno[[1]])
pheno_values <- as.numeric(pheno[[2]])

if (anyNA(pheno_values))
    stop("NA/non-numeric phenotype values found")

if (anyDuplicated(pheno_ids))
    stop("Duplicated phenotype IDs")


# Map phenotype by sample ID
pheno_map <- setNames(
    pheno_values,
    pheno_ids
)

y_train <- unname(
    pheno_map[train$ids]
)

y_test <- unname(
    pheno_map[test$ids]
)


if (anyNA(y_train))
    stop("Missing phenotype for training sample")

if (anyNA(y_test))
    stop("Missing phenotype for test sample")


cat("phenotype train:", length(y_train), "\n")
cat("phenotype test :", length(y_test), "\n")


# ============================================================
# FIT RR-BLUP
# ============================================================

cat("\n===== FIT RR-BLUP =====\n")

fit <- rrBLUP::mixed.solve(
    y = y_train,
    Z = train$G,
    SE = FALSE,
    return.Hinv = FALSE
)

cat("Vu   =", fit$Vu, "\n")
cat("Ve   =", fit$Ve, "\n")
cat("beta =", as.numeric(fit$beta), "\n")

if (fit$Vu != 0)
    cat("lambda =", fit$Ve / fit$Vu, "\n")


# ============================================================
# PREDICTION
# ============================================================

cat("\n===== PREDICTION =====\n")

pred <- as.numeric(
    test$G %*% fit$u + as.numeric(fit$beta)
)


if (length(pred) != length(y_test))
    stop("Prediction length does not match test phenotype length")


# ============================================================
# METRICS
# ============================================================

pearson <- cor(
    y_test,
    pred,
    method = "pearson"
)

spearman <- cor(
    y_test,
    pred,
    method = "spearman"
)

rmse <- sqrt(
    mean((pred - y_test)^2)
)

mae <- mean(
    abs(pred - y_test)
)

bias <- mean(
    pred - y_test
)


cat("Pearson  =", pearson, "\n")
cat("Spearman =", spearman, "\n")
cat("RMSE     =", rmse, "\n")
cat("MAE      =", mae, "\n")
cat("Bias     =", bias, "\n")


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

outdir <- file.path(
    "results",
    "rrBLUP",
    trait,
    geno_type
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# WRITE PREDICTIONS
# ============================================================

pred_file <- file.path(
    outdir,
    paste0(
        "outer_",
        outer,
        "_predictions.tsv"
    )
)

pred_out <- data.frame(
    ID = test$ids,
    observed = y_test,
    predicted = pred
)

write.table(
    pred_out,
    pred_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# WRITE METRICS
# ============================================================

metrics_file <- file.path(
    outdir,
    paste0(
        "outer_",
        outer,
        "_metrics.tsv"
    )
)

metrics_out <- data.frame(
    trait = trait,
    genotype = geno_type,
    outer = outer,
    n_train = nrow(train$G),
    n_test = nrow(test$G),
    n_markers = ncol(train$G),
    Pearson = pearson,
    Spearman = spearman,
    RMSE = rmse,
    MAE = mae,
    Bias = bias,
    Vu = fit$Vu,
    Ve = fit$Ve,
    beta = as.numeric(fit$beta)
)

write.table(
    metrics_out,
    metrics_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# DONE
# ============================================================

cat("\n===== DONE =====\n")
cat("predictions:", pred_file, "\n")
cat("metrics    :", metrics_file, "\n")
