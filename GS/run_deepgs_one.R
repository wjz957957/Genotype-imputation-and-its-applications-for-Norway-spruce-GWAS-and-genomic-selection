args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript run_deepgs_one.R <BB|H6> <chip|combined> <01-10>")
}

trait <- args[1]
geno_type <- args[2]
outer <- sprintf("%02d", as.integer(args[3]))

if (!trait %in% c("BB", "H6"))
    stop("trait must be BB or H6")

if (!geno_type %in% c("chip", "combined"))
    stop("geno_type must be chip or combined")

cat("========================================\n")
cat("              DeepGS\n")
cat("========================================\n")
cat("trait     :", trait, "\n")
cat("genotype  :", geno_type, "\n")
cat("outer     :", outer, "\n")

suppressPackageStartupMessages(
    library(DeepGS)
)

cat("DeepGS version:", as.character(packageVersion("DeepGS")), "\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base_dir <- file.path(
    "final_genotypes",
    paste0("outer_", outer),
    trait
)

train_file <- file.path(
    base_dir,
    paste0(geno_type, "_train.raw")
)

test_file <- file.path(
    base_dir,
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

outdir <- file.path(
    "deepgs_results",
    trait,
    geno_type
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pred_file <- file.path(
    outdir,
    paste0("outer_", outer, "_predictions.tsv")
)

metrics_file <- file.path(
    outdir,
    paste0("outer_", outer, "_metrics.tsv")
)

complete_file <- file.path(
    outdir,
    paste0("outer_", outer, ".complete")
)

# ------------------------------------------------------------
# Check existing result
# ------------------------------------------------------------

if (file.exists(pred_file) &&
    file.exists(metrics_file) &&
    file.exists(complete_file)) {

    cat("\nResult already exists:\n")
    cat(metrics_file, "\n")
    cat("Skipping.\n")

    quit(save = "no", status = 0)
}

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

required_files <- c(
    train_file,
    test_file,
    cv_train_file,
    cv_test_file,
    pheno_file
)

for (f in required_files) {
    if (!file.exists(f))
        stop("Missing file: ", f)
}

# ------------------------------------------------------------
# Special case:
# H6 outer folds without QTN
#
# If combined is identical to CHIP because no QTN exists,
# do NOT retrain the same genotype twice.
# The combined result will inherit the CHIP prediction.
# ------------------------------------------------------------

if (geno_type == "combined") {

    extra_train_file <- file.path(
        base_dir,
        "extra_train.raw"
    )

    extra_test_file <- file.path(
        base_dir,
        "extra_test.raw"
    )

    if (!file.exists(extra_train_file) ||
        !file.exists(extra_test_file)) {

        cat("\n===== NO QTN DETECTED =====\n")
        cat("No extra_train.raw / extra_test.raw found.\n")

        chip_outdir <- file.path(
            "deepgs_results",
            trait,
            "chip"
        )

        chip_pred_file <- file.path(
            chip_outdir,
            paste0("outer_", outer, "_predictions.tsv")
        )

        chip_metrics_file <- file.path(
            chip_outdir,
            paste0("outer_", outer, "_metrics.tsv")
        )

        chip_complete_file <- file.path(
            chip_outdir,
            paste0("outer_", outer, ".complete")
        )

        if (!file.exists(chip_pred_file) ||
            !file.exists(chip_metrics_file) ||
            !file.exists(chip_complete_file)) {

            stop(
                "No-QTN fold detected, but CHIP DeepGS result is missing. ",
                "Run CHIP first: ", trait, " outer_", outer
            )
        }

        cat("Using CHIP result because CHIP+QTN genotype is identical to CHIP.\n")

        pred_tmp <- read.table(
            chip_pred_file,
            header = TRUE,
            sep = "\t",
            quote = "",
            comment.char = "",
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        metrics_tmp <- read.table(
            chip_metrics_file,
            header = TRUE,
            sep = "\t",
            quote = "",
            comment.char = "",
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        metrics_tmp$genotype <- "combined"
        metrics_tmp$source <- "identical_to_chip_no_QTN"

        write.table(
            pred_tmp,
            pred_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        write.table(
            metrics_tmp,
            metrics_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        writeLines(
            "combined genotype identical to CHIP; no QTN available",
            complete_file
        )

        cat("Created:", metrics_file, "\n")
        cat("DONE.\n")

        quit(save = "no", status = 0)
    }
}

# ------------------------------------------------------------
# Read genotype
# ------------------------------------------------------------

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
        stop("Expected at least 7 columns: ", path)

    ids <- as.character(d[[2]])

    markers <- colnames(d)[7:ncol(d)]

    if (anyDuplicated(ids))
        stop("Duplicated sample IDs: ", path)

    if (anyDuplicated(markers))
        stop("Duplicated marker names: ", path)

    G <- data.matrix(
        d[, 7:ncol(d), drop = FALSE]
    )

    if (anyNA(G))
        stop("NA found in genotype matrix: ", path)

    if (any(!(G %in% c(0, 1, 2))))
        stop(
            "Unexpected genotype value found in: ",
            path,
            ". DeepGS input must contain 0/1/2."
        )

    list(
        ids = ids,
        markers = markers,
        G = G
    )
}

cat("\n===== READ GENOTYPES =====\n")

train <- read_genotype(train_file)
test <- read_genotype(test_file)

cat("outer train samples :", nrow(train$G), "\n")
cat("outer test samples  :", nrow(test$G), "\n")
cat("markers             :", ncol(train$G), "\n")

if (!identical(train$markers, test$markers))
    stop("Train/test marker columns are not identical")

# ------------------------------------------------------------
# Read outer CV sample IDs
# ------------------------------------------------------------

cat("\n===== OUTER SAMPLE ORDER CHECK =====\n")

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
    stop("Training sample order does not match outer_train.samples")

if (!identical(test$ids, cv_test_ids))
    stop("Testing sample order does not match outer_test.samples")

cat("train order: PASS\n")
cat("test order : PASS\n")

# ------------------------------------------------------------
# Read phenotype
# ------------------------------------------------------------

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
    stop("NA/non-numeric phenotype values")

if (anyDuplicated(pheno_ids))
    stop("Duplicated phenotype IDs")

pheno_map <- setNames(
    pheno_values,
    pheno_ids
)

y_outer_train <- unname(
    pheno_map[train$ids]
)

y_test <- unname(
    pheno_map[test$ids]
)

if (anyNA(y_outer_train))
    stop("Missing phenotype for outer training sample")

if (anyNA(y_test))
    stop("Missing phenotype for outer test sample")

cat("phenotype outer train:", length(y_outer_train), "\n")
cat("phenotype test       :", length(y_test), "\n")

# ------------------------------------------------------------
# Fixed inner validation set
# ------------------------------------------------------------

cat("\n===== INNER VALIDATION =====\n")

validation_dir <- file.path(
    "deepgs_results",
    "validation"
)

dir.create(
    validation_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

validation_file <- file.path(
    validation_dir,
    paste0("outer_", outer, "_valid.samples")
)

if (file.exists(validation_file)) {

    valid_ids <- scan(
        validation_file,
        what = "",
        quiet = TRUE
    )

    valid_idx <- match(
        valid_ids,
        train$ids
    )

    if (anyNA(valid_idx))
        stop(
            "Validation IDs contain samples not found in outer training set"
        )

} else {

    # Same seed for CHIP and CHIP+QTN within the same outer fold.
    validation_seed <- 20261000 + as.integer(outer)

    set.seed(validation_seed)

    n_outer_train <- nrow(train$G)

    n_valid <- floor(
        n_outer_train * 0.10
    )

    valid_idx <- sample(
        seq_len(n_outer_train),
        size = n_valid,
        replace = FALSE
    )

    valid_ids <- train$ids[valid_idx]

    writeLines(
        valid_ids,
        validation_file
    )
}

valid_idx <- as.integer(valid_idx)

fit_idx <- setdiff(
    seq_len(nrow(train$G)),
    valid_idx
)

cat("validation samples:", length(valid_idx), "\n")
cat("training samples  :", length(fit_idx), "\n")

# ------------------------------------------------------------
# Prepare matrices
# ------------------------------------------------------------

trainMat <- train$G[fit_idx, , drop = FALSE]
trainPheno <- y_outer_train[fit_idx]

validMat <- train$G[valid_idx, , drop = FALSE]
validPheno <- y_outer_train[valid_idx]

testMat <- test$G

cat("DeepGS training matrix :", nrow(trainMat), "x", ncol(trainMat), "\n")
cat("DeepGS validation matrix:", nrow(validMat), "x", ncol(validMat), "\n")
cat("DeepGS test matrix      :", nrow(testMat), "x", ncol(testMat), "\n")

# ------------------------------------------------------------
# DeepGS architecture
# Official DeepGS 1.2 example:
# 1 x 18 convolution
# 8 filters
# 1 x 4 max pooling
# 32 -> 1 fully connected layers
# dropout = 0.2, 0.1, 0.05
# ------------------------------------------------------------

conv_kernel <- c("1*18")
conv_stride <- c("1*1")
conv_num_filter <- c(8)

pool_act_type <- c("relu")
pool_type <- c("max")
pool_kernel <- c("1*4")
pool_stride <- c("1*4")

fullayer_num_hidden <- c(32, 1)
fullayer_act_type <- c("sigmoid")

drop_float <- c(
    0.2,
    0.1,
    0.05
)

cnnFrame <- list(
    conv_kernel = conv_kernel,
    conv_num_filter = conv_num_filter,
    conv_stride = conv_stride,
    pool_act_type = pool_act_type,
    pool_type = pool_type,
    pool_kernel = pool_kernel,
    pool_stride = pool_stride,
    fullayer_num_hidden = fullayer_num_hidden,
    fullayer_act_type = fullayer_act_type,
    drop_float = drop_float
)

markerImage <- paste0(
    "1*",
    ncol(trainMat)
)

# ------------------------------------------------------------
# Trait-specific training rounds
#
# BB  : 6000
# H6  : 200
#
# CHIP and CHIP+QTN use the same rounds within each trait.
# ------------------------------------------------------------

num_round <- if (trait == "BB") 6000 else 200

cat("\n===== DEEPGS PARAMETERS =====\n")
cat("markerImage       :", markerImage, "\n")
cat("num_round         :", num_round, "\n")
cat("array_batch_size  : 30\n")
cat("learning_rate     : 0.01\n")
cat("momentum          : 0.5\n")
cat("weight_decay      : 0.00001\n")
cat("randomseeds       : 0\n")
cat("initializer_idx   : 0.01\n")
cat("eval_metric       : mae\n")
cat("device_type       : cpu\n")

# ------------------------------------------------------------
# Train
# ------------------------------------------------------------

cat("\n===== FIT DEEPGS =====\n")
cat("Training started...\n")

fit <- train_deepGSModel(
    trainMat = trainMat,
    trainPheno = trainPheno,
    validMat = validMat,
    validPheno = validPheno,
    markerImage = markerImage,
    cnnFrame = cnnFrame,
    device_type = "cpu",
    gpuNum = 1,
    eval_metric = "mae",
    num_round = num_round,
    array_batch_size = 30,
    learning_rate = 0.01,
    momentum = 0.5,
    wd = 0.00001,
    randomseeds = 0,
    initializer_idx = 0.01,
    verbose = TRUE
)

cat("\n===== PREDICTION =====\n")

pred <- predict_GSModel(
    GSModel = fit,
    testMat = testMat,
    markerImage = markerImage
)

pred <- as.numeric(pred)

if (length(pred) != length(y_test))
    stop(
        "Prediction length does not match test phenotype length"
    )

if (anyNA(pred))
    stop("NA found in DeepGS predictions")

# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------

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
    mean(
        (pred - y_test)^2
    )
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

# ------------------------------------------------------------
# Save predictions
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Save metrics
# ------------------------------------------------------------

metrics_out <- data.frame(
    trait = trait,
    genotype = geno_type,
    outer = outer,

    n_outer_train = nrow(train$G),
    n_inner_train = nrow(trainMat),
    n_validation = nrow(validMat),
    n_test = nrow(testMat),

    n_markers = ncol(trainMat),

    num_round = num_round,
    array_batch_size = 30,
    learning_rate = 0.01,
    momentum = 0.5,
    wd = 0.00001,
    randomseeds = 0,
    initializer_idx = 0.01,

    Pearson = pearson,
    Spearman = spearman,
    RMSE = rmse,
    MAE = mae,
    Bias = bias
)

write.table(
    metrics_out,
    metrics_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ------------------------------------------------------------
# Save completion flag
# ------------------------------------------------------------

writeLines(
    paste(
        "DeepGS completed",
        trait,
        geno_type,
        paste0("outer_", outer)
    ),
    complete_file
)

cat("\n===== DONE =====\n")
cat("predictions:", pred_file, "\n")
cat("metrics    :", metrics_file, "\n")
cat("validation :", validation_file, "\n")
cat("complete   :", complete_file, "\n")
