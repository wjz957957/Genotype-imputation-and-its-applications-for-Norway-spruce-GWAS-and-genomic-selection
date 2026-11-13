#!/usr/bin/env Rscript

suppressPackageStartupMessages(
    library(rrBLUP)
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript 07_run_chip_rrblup.R WORK OUTER INNER")
}

WORK  <- args[1]
OUTER <- args[2]
INNER <- args[3]

SETS <- c(
    "50", "100", "200", "300", "400", "500",
    "1000", "2000", "3000", "4000", "5000",
    "8000", "10000", "20000", "30000",
    "40000", "LD_ALL", "ALL"
)

TRAITS <- c(
    BB = "Trait-Adj_Budburst.txt",
    H6 = "Trait-Adj_H6.txt"
)

CV_DIR <- file.path(
    WORK, "02_CV",
    paste0("outer_", OUTER),
    paste0("inner_", INNER)
)

BASE <- file.path(
    WORK, "07_chip_RRBLUP",
    paste0("outer_", OUTER),
    paste0("inner_", INNER)
)

GENO_DIR <- file.path(BASE, "genotype")

MARKER_DIR <- file.path(
    WORK, "06_chip_clumping",
    paste0("outer_", OUTER),
    paste0("inner_", INNER)
)

dir.create(BASE, recursive = TRUE, showWarnings = FALSE)

##################################################
# Read sample IDs
##################################################

read_ids <- function(path) {

    if (!file.exists(path)) {
        stop("Missing file: ", path)
    }

    ids <- trimws(readLines(path, warn = FALSE))
    ids <- ids[nzchar(ids)]

    if (length(ids) == 0) {
        stop("Empty ID file: ", path)
    }

    if (anyDuplicated(ids)) {
        stop("Duplicate IDs: ", path)
    }

    ids
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

if (length(intersect(TRAIN_IDS, VALID_IDS)) > 0) {
    stop("ERROR: inner train / valid overlap")
}

if (
    length(intersect(TRAIN_IDS, OUTER_TEST)) > 0 ||
    length(intersect(VALID_IDS, OUTER_TEST)) > 0
) {
    stop("ERROR: outer test leakage")
}

cat("====================================\n")
cat("Outer:", OUTER, "\n")
cat("Inner:", INNER, "\n")
cat("Train samples:", length(TRAIN_IDS), "\n")
cat("Valid samples:", length(VALID_IDS), "\n")
cat("====================================\n")

##################################################
# Read phenotype
##################################################

read_pheno <- function(path) {

    p <- read.table(
        path,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (ncol(p) != 3) {
        stop("Phenotype must have 3 columns: ", path)
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
# Read PLINK raw genotype
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

    if (ncol(x) < 7) {
        stop("Invalid PLINK raw: ", path)
    }

    if (!identical(names(x)[1:2], c("FID", "IID"))) {
        stop("Unexpected PLINK raw columns: ", path)
    }

    ids <- as.character(x$IID)

    if (anyDuplicated(ids)) {
        stop("Duplicate genotype IDs: ", path)
    }

    if (any(as.character(x$FID) != ids)) {
        stop("Genotype FID/IID mismatch: ", path)
    }

    columns <- names(x)[-(1:6)]

    # Example:
    # 12:279019076:T:G_G
    # SNP ID = 12:279019076:T:G
    # Counted allele = G

    snp_ids <- sub("_[^_]+$", "", columns)
    alleles <- sub("^.*_", "", columns)

    if (anyDuplicated(snp_ids)) {
        stop("Duplicate SNP IDs in raw: ", path)
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
# Analyze each trait
##################################################

for (TRAIT in names(TRAITS)) {

    cat("\n====================================\n")
    cat("Trait:", TRAIT, "\n")
    cat("====================================\n")

    PHENO_FILE <- file.path(
        WORK, "00_input", TRAITS[[TRAIT]]
    )

    phenotype <- read_pheno(PHENO_FILE)

    if (!all(c(TRAIN_IDS, VALID_IDS) %in% names(phenotype))) {
        stop("Missing phenotype IDs for ", TRAIT)
    }

    results <- vector("list", length(SETS))

    ##################################################
    # Analyze each SNP set
    ##################################################

    for (i in seq_along(SETS)) {

        SET <- SETS[i]

        cat("\nRunning ", TRAIT, " chip_", SET, "\n", sep = "")

        SNP_FILE <- file.path(
            MARKER_DIR,
            TRAIT,
            "marker_sets",
            paste0("chip_", SET, ".snps")
        )

        SELECTED <- read_ids(SNP_FILE)

        if (SET %in% as.character(c(
            50, 100, 200, 300, 400, 500,
            1000, 2000, 3000, 4000, 5000,
            8000, 10000, 20000, 30000, 40000
        ))) {

            if (length(SELECTED) != as.integer(SET)) {
                stop("Incorrect marker set size: ", SNP_FILE)
            }
        }

        TRAIN_FILE <- file.path(
            GENO_DIR,
            paste0("chip_", TRAIT, "_", SET, "_train.raw")
        )

        VALID_FILE <- file.path(
            GENO_DIR,
            paste0("chip_", TRAIT, "_", SET, "_valid.raw")
        )

        TRAIN <- read_raw(TRAIN_FILE)
        VALID <- read_raw(VALID_FILE)

        ##################################################
        # Check sample membership
        ##################################################

        if (!setequal(TRAIN$ids, TRAIN_IDS)) {
            stop("Incorrect training samples: ", TRAIN_FILE)
        }

        if (!setequal(VALID$ids, VALID_IDS)) {
            stop("Incorrect validation samples: ", VALID_FILE)
        }

        ##################################################
        # Check SNP membership
        ##################################################

        if (!setequal(TRAIN$snp_ids, SELECTED)) {
            stop("Training SNP set mismatch: ", TRAIN_FILE)
        }

        if (!setequal(VALID$snp_ids, SELECTED)) {
            stop("Validation SNP set mismatch: ", VALID_FILE)
        }

        # Always align validation to training SNP order.

        IDX_TRAIN <- match(SELECTED, TRAIN$snp_ids)
        IDX_VALID <- match(SELECTED, VALID$snp_ids)

        X_train <- TRAIN$G[, IDX_TRAIN, drop = FALSE]
        X_valid <- VALID$G[, IDX_VALID, drop = FALSE]

        A1_train <- TRAIN$alleles[IDX_TRAIN]
        A1_valid <- VALID$alleles[IDX_VALID]

        ##################################################
        # Check and align allele dosage directions
        ##################################################

        parts <- strsplit(
            SELECTED, ":", fixed = TRUE
        )

        if (any(lengths(parts) != 4)) {
            stop("Unexpected SNP ID format in ", SNP_FILE)
        }

        REF <- vapply(parts, function(x) x[3], character(1))
        ALT <- vapply(parts, function(x) x[4], character(1))

        if (any(!(A1_train == REF | A1_train == ALT))) {
            stop("Training allele mismatch: ", TRAIN_FILE)
        }

        if (any(!(A1_valid == REF | A1_valid == ALT))) {
            stop("Validation allele mismatch: ", VALID_FILE)
        }

        FLIP <- which(A1_train != A1_valid)

        if (length(FLIP) > 0) {

            X_valid[, FLIP] <-
                2 - X_valid[, FLIP, drop = FALSE]
        }

        ##################################################
        # Align phenotype to genotype sample order
        ##################################################

        y_train <- unname(phenotype[TRAIN$ids])
        y_valid <- unname(phenotype[VALID$ids])

        if (
            any(!is.finite(y_train)) ||
            any(!is.finite(y_valid))
        ) {
            stop("Missing or invalid phenotype: ", TRAIT)
        }

        ##################################################
        # Training-only missing value imputation
        # and zero-variance SNP removal
        ##################################################

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

        if (N_USED == 0) {
            stop("No informative SNPs: ", TRAIT, " ", SET)
        }

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

        ##################################################
        # RRBLUP
        ##################################################

        FIT <- rrBLUP::mixed.solve(
            y = y_train,
            Z = X_train
        )

        # Include model intercept.
        prediction <- as.numeric(
            FIT$beta[1] + X_valid %*% FIT$u
        )

        if (any(!is.finite(prediction))) {
            stop("Non-finite prediction: ", TRAIT, " ", SET)
        }

        ACC <- cor(prediction, y_valid)

        if (!is.finite(ACC)) {
            stop("Invalid prediction correlation: ", TRAIT, " ", SET)
        }

        RMSE <- sqrt(
            mean((prediction - y_valid)^2)
        )

        ##################################################
        # Save result
        ##################################################

        results[[i]] <- data.frame(
            Outer = OUTER,
            Inner = INNER,
            Trait = TRAIT,
            Marker = paste0("chip_", SET),
            N_selected = length(SELECTED),
            N_used = N_USED,
            N_train = length(y_train),
            N_valid = length(y_valid),
            Alleles_flipped = length(FLIP),
            Accuracy = ACC,
            RMSE = RMSE
        )

        cat(
            "PASS:",
            "selected=", length(SELECTED),
            "used=", N_USED,
            "train=", length(y_train),
            "valid=", length(y_valid),
            "flipped=", length(FLIP),
            "cor=", ACC,
            "\n"
        )

        rm(
            TRAIN, VALID,
            X_train, X_valid,
            FIT, prediction
        )

        gc(verbose = FALSE)
    }

    RESULT <- do.call(rbind, results)

    OUT_FILE <- file.path(
        BASE,
        paste0(TRAIT, "_chip_accuracy.txt")
    )

    write.table(
        RESULT,
        file = OUT_FILE,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("SAVED:", OUT_FILE, "\n")
}

cat("\nDONE outer_", OUTER, " inner_", INNER, "\n", sep = "")
