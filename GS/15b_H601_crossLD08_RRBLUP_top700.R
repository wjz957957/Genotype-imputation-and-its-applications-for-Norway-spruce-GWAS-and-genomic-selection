#!/usr/bin/env Rscript

suppressPackageStartupMessages(
    library(rrBLUP)
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
    stop("Usage: Rscript 15b_H601_crossLD08_RRBLUP_top700.R WORK OUTER INNER")
}

WORK  <- args[1]
OUTER <- args[2]
INNER <- args[3]

if (!identical(OUTER, "01") || !INNER %in% c("01", "02", "03")) {
    stop("This pilot is restricted to H6 outer_01 inner_01..inner_03")
}

COUNTS <- c(0L, 1L, 2L, 3L, 4L, 5L, 10L, 20L, 30L, 40L,
            50L, 75L, 100L, 150L, 200L, 300L, 400L, 500L,
            600L, 700L)

TRAITS <- c(
    H6 = "Trait-Adj_H6.txt"
)

CV_DIR <- file.path(
    WORK, "02_CV",
    paste0("outer_", OUTER),
    paste0("inner_", INNER)
)

CHIP_DIR <- file.path(
    WORK, "07_chip_RRBLUP",
    paste0("outer_", OUTER),
    paste0("inner_", INNER)
)

CHIP_GENO <- file.path(CHIP_DIR, "genotype")

CHIP_MARKERS <- file.path(
    WORK, "06_chip_clumping",
    paste0("outer_", OUTER),
    paste0("inner_", INNER)
)

OPT_DIR <- file.path(
    WORK, "07_chip_RRBLUP", "optimization",
    paste0("outer_", OUTER)
)

# New marker sets from inner_train-only selected-CHIP vs QTN LD filtering.
EXTRA_DIR <- file.path(
    WORK, "08_imputed_GWAS",
    "08_crosschip_LD_pilot",
    paste0("outer_", OUTER),
    paste0("inner_", INNER),
    "H6", "r2_0p8", "marker_sets"
)

# Both RAW matrices and prediction outputs live under this isolated pilot root.
# Never use the old Step 9a matrix, whose SNPs are from the unfiltered ranking.
BASE <- file.path(
    WORK, "09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_H6_pilot",
    paste0("outer_", OUTER),
    paste0("inner_", INNER)
)
RESULT_BASE <- BASE

IMP_GENO <- file.path(BASE, "genotype")

dir.create(BASE, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_BASE, recursive = TRUE, showWarnings = FALSE)

##################################################
# 1. Read sample / marker IDs
##################################################

read_ids <- function(path, allow_empty = FALSE) {

    if (!file.exists(path)) {
        stop("Missing file: ", path)
    }

    x <- trimws(readLines(path, warn = FALSE))
    x <- x[nzchar(x)]

    if (!allow_empty && length(x) == 0L) {
        stop("Empty ID file: ", path)
    }

    if (anyDuplicated(x)) {
        stop("Duplicate IDs: ", path)
    }

    x
}

TRAIN_IDS <- read_ids(
    file.path(CV_DIR, "inner_train.samples")
)

VALID_IDS <- read_ids(
    file.path(CV_DIR, "inner_valid.samples")
)

OUTER_TEST <- read_ids(
    file.path(
        WORK, "02_CV",
        paste0("outer_", OUTER),
        "outer_test.samples"
    )
)

if (length(intersect(TRAIN_IDS, VALID_IDS)) > 0L ||
    length(intersect(TRAIN_IDS, OUTER_TEST)) > 0L ||
    length(intersect(VALID_IDS, OUTER_TEST)) > 0L) {
    stop("ERROR: Nested CV sample overlap")
}

if (length(TRAIN_IDS) < 500L ||
    length(VALID_IDS) < 250L) {
    stop("Unexpected CV sample counts")
}

cat("====================================\n")
cat("Step 15b: H6 selected-CHIP + cross-LD-filtered QTN RRBLUP\n")
cat("Outer:", OUTER, "\n")
cat("Inner:", INNER, "\n")
cat("Train:", length(TRAIN_IDS), "\n")
cat("Valid:", length(VALID_IDS), "\n")
cat("====================================\n")

##################################################
# 2. Read phenotype: identical to Step 7
##################################################

read_pheno <- function(path) {

    p <- read.table(
        path,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (ncol(p) != 3L) {
        stop("Phenotype must have three columns: ", path)
    }

    if (any(p[[1]] != p[[2]])) {
        stop("Phenotype FID/IID mismatch: ", path)
    }

    ids <- as.character(p[[2]])
    y <- suppressWarnings(as.numeric(p[[3]]))

    if (anyDuplicated(ids)) {
        stop("Duplicate phenotype IDs: ", path)
    }

    names(y) <- ids
    y
}

##################################################
# 3. Read PLINK RAW: identical to Step 7
##################################################

read_raw <- function(path) {

    if (!file.exists(path)) {
        stop("Missing genotype: ", path)
    }

    if (requireNamespace("data.table", quietly = TRUE)) {

        x <- data.table::fread(
            path,
            header = TRUE,
            check.names = FALSE,
            data.table = FALSE,
            na.strings = "NA"
        )

    } else {

        x <- read.table(
            path,
            header = TRUE,
            check.names = FALSE,
            stringsAsFactors = FALSE,
            na.strings = "NA",
            comment.char = ""
        )
    }

    if (ncol(x) < 7L) {
        stop("Invalid PLINK RAW: ", path)
    }

    if (!identical(names(x)[1:2], c("FID", "IID"))) {
        stop("Unexpected RAW columns: ", path)
    }

    ids <- as.character(x$IID)

    if (anyDuplicated(ids)) {
        stop("Duplicate genotype IDs: ", path)
    }

    if (any(as.character(x$FID) != ids)) {
        stop("Genotype FID/IID mismatch: ", path)
    }

    columns <- names(x)[-(1:6)]

    snp_ids <- sub("_[^_]+$", "", columns)
    alleles <- sub("^.*_", "", columns)

    if (anyDuplicated(snp_ids)) {
        stop("Duplicate SNP IDs in RAW: ", path)
    }

    G <- as.matrix(x[, -(1:6), drop = FALSE])

    storage.mode(G) <- "numeric"

    if (any(G < 0 | G > 2, na.rm = TRUE)) {
        stop("Invalid genotype dosage: ", path)
    }

    if (any(is.infinite(G))) {
        stop("Infinite genotype dosage: ", path)
    }

    rownames(G) <- ids

    list(
        ids = ids,
        snp_ids = snp_ids,
        alleles = alleles,
        G = G
    )
}

##################################################
# 4. Check raw sample membership
##################################################

check_samples <- function(obj, expected, label) {

    if (!setequal(obj$ids, expected)) {
        stop("Incorrect sample membership: ", label)
    }

    if (length(obj$ids) != length(expected)) {
        stop("Incorrect sample count: ", label)
    }
}

##################################################
# 5. Align selected SNPs and allele dosage
##################################################

align_markers <- function(
    train,
    valid,
    selected,
    label
) {

    if (!setequal(train$snp_ids, selected)) {
        stop("Training marker mismatch: ", label)
    }

    if (!setequal(valid$snp_ids, selected)) {
        stop("Validation marker mismatch: ", label)
    }

    idx_train <- match(selected, train$snp_ids)
    idx_valid <- match(selected, valid$snp_ids)

    X_train <- train$G[, idx_train, drop = FALSE]
    X_valid <- valid$G[, idx_valid, drop = FALSE]

    a_train <- train$alleles[idx_train]
    a_valid <- valid$alleles[idx_valid]

    parts <- strsplit(selected, ":", fixed = TRUE)

    if (any(lengths(parts) != 4L)) {
        stop("Unexpected SNP ID format: ", label)
    }

    REF <- vapply(
        parts, function(z) z[3], character(1)
    )

    ALT <- vapply(
        parts, function(z) z[4], character(1)
    )

    if (any(!(a_train == REF | a_train == ALT))) {
        stop("Training counted allele mismatch: ", label)
    }

    if (any(!(a_valid == REF | a_valid == ALT))) {
        stop("Validation counted allele mismatch: ", label)
    }

    FLIP <- which(a_train != a_valid)

    if (length(FLIP) > 0L) {
        X_valid[, FLIP] <-
            2 - X_valid[, FLIP, drop = FALSE]
    }

    list(
        train = X_train,
        valid = X_valid,
        n_flipped = length(FLIP)
    )
}

##################################################
# 6. RRBLUP: preserve Step 7 method
##################################################

fit_rrblup <- function(
    X_train,
    X_valid,
    y_train,
    y_valid,
    n_chip
) {

    means <- colMeans(X_train, na.rm = TRUE)

    variance <- colMeans(
        X_train^2, na.rm = TRUE
    ) - means^2

    KEEP_SNP <- (
        is.finite(means) &
        is.finite(variance) &
        variance > 1e-12
    )

    N_USED <- sum(KEEP_SNP)

    if (N_USED == 0L) {
        stop("No informative SNPs")
    }

    N_CHIP_USED <- sum(
        KEEP_SNP[seq_len(n_chip)]
    )

    N_EXTRA_USED <- N_USED - N_CHIP_USED

    X_train <- X_train[, KEEP_SNP, drop = FALSE]
    X_valid <- X_valid[, KEEP_SNP, drop = FALSE]

    means <- means[KEEP_SNP]

    if (anyNA(X_train)) {

        pos <- which(
            is.na(X_train), arr.ind = TRUE
        )

        X_train[pos] <- means[pos[, 2]]
    }

    if (anyNA(X_valid)) {

        pos <- which(
            is.na(X_valid), arr.ind = TRUE
        )

        X_valid[pos] <- means[pos[, 2]]
    }

    FIT <- rrBLUP::mixed.solve(
        y = y_train,
        Z = X_train
    )

    prediction <- as.numeric(
        FIT$beta[1] + X_valid %*% FIT$u
    )

    if (any(!is.finite(prediction))) {
        stop("Non-finite RRBLUP prediction")
    }

    ACC <- cor(prediction, y_valid)

    if (!is.finite(ACC)) {
        stop("Invalid prediction correlation")
    }

    RMSE <- sqrt(
        mean((prediction - y_valid)^2)
    )

    if (!is.finite(RMSE)) {
        stop("Invalid prediction RMSE")
    }

    list(
        accuracy = ACC,
        rmse = RMSE,
        n_used = N_USED,
        n_chip_used = N_CHIP_USED,
        n_extra_used = N_EXTRA_USED
    )
}

##################################################
# 7. Compare with Step 7 baseline
##################################################

check_baseline <- function(
    result,
    old_row,
    n_selected,
    chip_flipped,
    n_train,
    n_valid,
    trait
) {

    if (nrow(old_row) != 1L) {
        stop("Missing or duplicate Step 7 baseline: ", trait)
    }

    int_checks <- c(
        selected = result$n_used <= n_selected,
        n_used = result$n_used == old_row$N_used,
        n_selected = n_selected == old_row$N_selected,
        flipped = chip_flipped == old_row$Alleles_flipped,
        n_train = n_train == old_row$N_train,
        n_valid = n_valid == old_row$N_valid
    )

    if (!all(int_checks)) {
        print(int_checks)
        stop("Step 7 baseline count mismatch: ", trait)
    }

    # Step 7 writes floating-point numbers as text.
    # Allow only small rounding/numerical differences.

    acc_diff <- abs(
        result$accuracy - old_row$Accuracy
    )

    rmse_diff <- abs(
        result$rmse - old_row$RMSE
    )

    acc_tolerance <- 1e-5
    rmse_tolerance <- 1e-5 *
        max(1, abs(old_row$RMSE))

    cat(
        "Baseline comparison:",
        trait,
        "accuracy difference=", acc_diff,
        "RMSE difference=", rmse_diff,
        "\n"
    )

    if (!is.finite(acc_diff) ||
        !is.finite(rmse_diff) ||
        acc_diff > acc_tolerance ||
        rmse_diff > rmse_tolerance) {

        stop(
            "ERROR: extra_0 does not reproduce Step 7: ",
            trait
        )
    }

    cat("PASS: extra_0 reproduces Step 7\n")
}

##################################################
# 8. Validate Step 9a
##################################################

if (!file.exists(
    file.path(BASE, "step9a.complete.flag")
)) {
    stop("Step 9a not completed")
}

if (!file.exists(
    file.path(dirname(EXTRA_DIR), "complete.flag")
)) {
    stop("Step 8.4 not completed")
}

##################################################
# 9. Process H6 only
##################################################

for (TRAIT in names(TRAITS)) {

    cat("\n====================================\n")
    cat("Trait:", TRAIT, "\n")
    cat("====================================\n")

    ########################################
    # 9.1. Phenotype
    ########################################

    PHENO_FILE <- file.path(
        WORK, "00_input", TRAITS[[TRAIT]]
    )

    phenotype <- read_pheno(PHENO_FILE)

    if (!all(
        c(TRAIN_IDS, VALID_IDS) %in% names(phenotype)
    )) {
        stop("Missing phenotype IDs: ", TRAIT)
    }

    ########################################
    # 9.2. Best chip configuration
    ########################################

    BEST_FILE <- file.path(
        OPT_DIR,
        paste0(TRAIT, "_best_marker.txt")
    )

    BEST <- read_ids(BEST_FILE)

    if (length(BEST) != 1L ||
        !startsWith(BEST, "chip_")) {
        stop("Invalid best chip setting: ", BEST_FILE)
    }

    CHIP_SET <- sub("^chip_", "", BEST)

    SELECTED_CHIP <- read_ids(
        file.path(
            CHIP_MARKERS,
            TRAIT,
            "marker_sets",
            paste0(BEST, ".snps")
        )
    )

    cat("Best chip:", BEST, "\n")
    cat("Selected chip SNPs:", length(SELECTED_CHIP), "\n")

    ########################################
    # 9.3. Read existing Step 7 chip RAW
    ########################################

    CHIP_TRAIN <- read_raw(
        file.path(
            CHIP_GENO,
            paste0(
                "chip_", TRAIT, "_",
                CHIP_SET, "_train.raw"
            )
        )
    )

    CHIP_VALID <- read_raw(
        file.path(
            CHIP_GENO,
            paste0(
                "chip_", TRAIT, "_",
                CHIP_SET, "_valid.raw"
            )
        )
    )

    check_samples(
        CHIP_TRAIN, TRAIN_IDS, "chip train"
    )

    check_samples(
        CHIP_VALID, VALID_IDS, "chip valid"
    )

    CHIP <- align_markers(
        CHIP_TRAIN,
        CHIP_VALID,
        SELECTED_CHIP,
        paste0(TRAIT, " chip")
    )

    ########################################
    # 9.4. Read new imputed RAW
    ########################################

    IMP_TRAIN <- read_raw(
        file.path(
            IMP_GENO,
            paste0(
                "imputed_", TRAIT,
                "_top700_train.raw"
            )
        )
    )

    IMP_VALID <- read_raw(
        file.path(
            IMP_GENO,
            paste0(
                "imputed_", TRAIT,
                "_top700_valid.raw"
            )
        )
    )

    check_samples(
        IMP_TRAIN, TRAIN_IDS, "imputed train"
    )

    check_samples(
        IMP_VALID, VALID_IDS, "imputed valid"
    )

    ########################################
    # 9.5. Align imputed sample row order
    #      to existing chip RAW
    ########################################

    train_row_idx <- match(
        CHIP_TRAIN$ids,
        IMP_TRAIN$ids
    )

    valid_row_idx <- match(
        CHIP_VALID$ids,
        IMP_VALID$ids
    )

    if (anyNA(train_row_idx) ||
        anyNA(valid_row_idx)) {
        stop("Cannot align chip/imputed samples")
    }

    IMP_TRAIN$G <- IMP_TRAIN$G[
        train_row_idx, , drop = FALSE
    ]

    IMP_VALID$G <- IMP_VALID$G[
        valid_row_idx, , drop = FALSE
    ]

    IMP_TRAIN$ids <- CHIP_TRAIN$ids
    IMP_VALID$ids <- CHIP_VALID$ids

    ########################################
    # 9.6. Validate all nested SNP sets
    ########################################

    EXTRA_PATH <- EXTRA_DIR

    TOP700 <- read_ids(
        file.path(
            EXTRA_PATH, "extra_700.snps"
        )
    )

    if (length(TOP700) != 700L) {
        stop("Incorrect top700 SNP count")
    }

    if (length(intersect(
        TOP700, SELECTED_CHIP
    )) > 0L) {
        stop("Chip and imputed marker overlap")
    }

    for (N in COUNTS) {

        IDS <- read_ids(
            file.path(
                EXTRA_PATH,
                paste0("extra_", N, ".snps")
            ),
            allow_empty = (N == 0L)
        )

        if (length(IDS) != N) {
            stop("Incorrect extra SNP count: ", N)
        }

        if (!identical(
            IDS,
            TOP700[seq_len(N)]
        )) {
            stop(
                "Extra SNP sets are not nested: ",
                TRAIT, " ", N
            )
        }
    }

    IMP <- align_markers(
        IMP_TRAIN,
        IMP_VALID,
        TOP700,
        paste0(TRAIT, " imputed")
    )

    ########################################
    # 9.7. Align phenotype
    ########################################

    y_train <- unname(
        phenotype[CHIP_TRAIN$ids]
    )

    y_valid <- unname(
        phenotype[CHIP_VALID$ids]
    )

    if (any(!is.finite(y_train)) ||
        any(!is.finite(y_valid))) {
        stop("Invalid phenotype: ", TRAIT)
    }

    ########################################
    # 9.8. Load Step 7 baseline
    ########################################

    OLD_FILE <- file.path(
        CHIP_DIR,
        paste0(TRAIT, "_chip_accuracy.txt")
    )

    OLD <- read.delim(
        OLD_FILE,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    OLD_ROW <- OLD[
        OLD$Marker == BEST,
        ,
        drop = FALSE
    ]

    if (nrow(OLD_ROW) != 1L ||
        OLD_ROW$Trait != TRAIT) {
        stop("Invalid Step 7 baseline: ", TRAIT)
    }

    ########################################
    # 9.9. Fit 20 RRBLUP models
    ########################################

    results <- vector("list", length(COUNTS))

    BASE_ACC <- NA_real_
    BASE_RMSE <- NA_real_

    N_CHIP <- length(SELECTED_CHIP)

    for (i in seq_along(COUNTS)) {

        N <- COUNTS[i]

        cat(
            "\nRunning",
            TRAIT,
            BEST,
            "+ extra_", N,
            "\n",
            sep = " "
        )

        if (N == 0L) {

            X_train <- CHIP$train
            X_valid <- CHIP$valid

        } else {

            X_train <- cbind(
                CHIP$train,
                IMP$train[
                    , seq_len(N), drop = FALSE
                ]
            )

            X_valid <- cbind(
                CHIP$valid,
                IMP$valid[
                    , seq_len(N), drop = FALSE
                ]
            )
        }

        ####################################
        # RRBLUP fit and validation
        ####################################

        FIT <- fit_rrblup(
            X_train = X_train,
            X_valid = X_valid,
            y_train = y_train,
            y_valid = y_valid,
            n_chip = N_CHIP
        )

        ####################################
        # Exact baseline QC at extra_0
        ####################################

        if (N == 0L) {

            check_baseline(
                result = FIT,
                old_row = OLD_ROW,
                n_selected = N_CHIP,
                chip_flipped = CHIP$n_flipped,
                n_train = length(y_train),
                n_valid = length(y_valid),
                trait = TRAIT
            )

            BASE_ACC <- FIT$accuracy
            BASE_RMSE <- FIT$rmse
        }

        ####################################
        # Save statistics
        ####################################

        results[[i]] <- data.frame(
            Outer = OUTER,
            Inner = INNER,
            Trait = TRAIT,
            Best_Chip = BEST,
            N_chip_selected = N_CHIP,
            N_chip_used = FIT$n_chip_used,
            Extra_Count = N,
            N_extra_used = FIT$n_extra_used,
            N_total_selected = N_CHIP + N,
            N_total_used = FIT$n_used,
            N_train = length(y_train),
            N_valid = length(y_valid),
            Chip_alleles_flipped = CHIP$n_flipped,
            Extra_alleles_flipped = if (N == 0L) {
                0L
            } else {
                sum(
                    TOP700[
                        seq_len(N)
                    ] %in% TOP700[
                        which(
                            IMP_TRAIN$alleles[
                                match(
                                    TOP700,
                                    IMP_TRAIN$snp_ids
                                )
                            ] !=
                            IMP_VALID$alleles[
                                match(
                                    TOP700,
                                    IMP_VALID$snp_ids
                                )
                            ]
                        )
                    ]
                )
            },
            Accuracy = FIT$accuracy,
            RMSE = FIT$rmse,
            Delta_Accuracy_vs_0 =
                FIT$accuracy - BASE_ACC,
            Delta_RMSE_vs_0 =
                FIT$rmse - BASE_RMSE
        )

        cat(
            "PASS:",
            "Extra=", N,
            "Used=", FIT$n_used,
            "Accuracy=", FIT$accuracy,
            "Delta=", FIT$accuracy - BASE_ACC,
            "\n"
        )

        rm(X_train, X_valid, FIT)

        gc(verbose = FALSE)
    }

    ########################################
    # 9.10. Save trait results
    ########################################

    RESULT <- do.call(rbind, results)

    if (!identical(
        RESULT$Extra_Count,
        COUNTS
    )) {
        stop("Extra count order mismatch")
    }

    OUTPUT <- file.path(
        RESULT_BASE,
        paste0(
            TRAIT,
            "_chip_plus_imputed_accuracy.tsv"
        )
    )

    write.table(
        RESULT,
        file = paste0(OUTPUT, ".tmp"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    if (!file.rename(
        paste0(OUTPUT, ".tmp"),
        OUTPUT
    )) {
        stop("Failed to publish result: ", OUTPUT)
    }

    cat("SAVED:", OUTPUT, "\n")

    rm(
        CHIP_TRAIN, CHIP_VALID,
        IMP_TRAIN, IMP_VALID,
        CHIP, IMP, RESULT, results
    )

    gc(verbose = FALSE)
}

cat("\nDONE Step9b outer_", OUTER,
    " inner_", INNER, "\n", sep = "")
