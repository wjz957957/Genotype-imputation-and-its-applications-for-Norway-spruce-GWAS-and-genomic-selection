#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: Rscript 07_summarize_chip_rrblup.R WORK")
}

WORK <- args[1]

BASE <- file.path(WORK, "07_chip_RRBLUP")
CV_DIR <- file.path(WORK, "02_CV")
OUTDIR <- file.path(BASE, "optimization")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUTERS <- sprintf("%02d", 1:10)
INNERS <- sprintf("%02d", 1:3)

TRAITS <- c("BB", "H6")

SETS <- paste0(
    "chip_",
    c(
        "50", "100", "200", "300", "400",
        "500", "1000", "2000", "3000",
        "4000", "5000", "8000",
        "10000", "20000", "30000",
        "40000", "LD_ALL", "ALL"
    )
)

EXPECTED_COLUMNS <- c(
    "Outer", "Inner", "Trait", "Marker",
    "N_selected", "N_used",
    "N_train", "N_valid",
    "Alleles_flipped", "Accuracy", "RMSE"
)

read_ids <- function(path) {

    if (!file.exists(path)) {
        stop("Missing sample file: ", path)
    }

    x <- trimws(readLines(path, warn = FALSE))
    x <- x[nzchar(x)]

    if (length(x) == 0 || anyDuplicated(x)) {
        stop("Empty or duplicate sample IDs: ", path)
    }

    x
}

save_table <- function(x, path) {

    write.table(
        x,
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
}

##################################################
# Load and validate all 60 accuracy files
##################################################

all_results <- list()
k <- 0L

for (outer in OUTERS) {

    outer_train <- read_ids(
        file.path(
            CV_DIR,
            paste0("outer_", outer),
            "outer_train.samples"
        )
    )

    valid_groups <- list()

    for (inner in INNERS) {

        fold_dir <- file.path(
            CV_DIR,
            paste0("outer_", outer),
            paste0("inner_", inner)
        )

        train_ids <- read_ids(
            file.path(fold_dir, "inner_train.samples")
        )

        valid_ids <- read_ids(
            file.path(fold_dir, "inner_valid.samples")
        )

        if (!setequal(train_ids, setdiff(outer_train, valid_ids))) {
            stop(
                "Incorrect inner train partition: outer_",
                outer, " inner_", inner
            )
        }

        if (length(intersect(train_ids, valid_ids)) != 0) {
            stop("Training / validation overlap")
        }

        valid_groups[[inner]] <- valid_ids

        for (trait in TRAITS) {

            accuracy_file <- file.path(
                BASE,
                paste0("outer_", outer),
                paste0("inner_", inner),
                paste0(trait, "_chip_accuracy.txt")
            )

            if (!file.exists(accuracy_file)) {
                stop("Missing accuracy file: ", accuracy_file)
            }

            x <- read.delim(
                accuracy_file,
                header = TRUE,
                sep = "\t",
                check.names = FALSE,
                stringsAsFactors = FALSE,
                colClasses = c(
                    "character", "character",
                    "character", "character",
                    "integer", "integer",
                    "integer", "integer",
                    "integer", "numeric", "numeric"
                )
            )

            if (!identical(names(x), EXPECTED_COLUMNS)) {
                stop("Unexpected columns: ", accuracy_file)
            }

            if (nrow(x) != length(SETS)) {
                stop("Expected 18 results: ", accuracy_file)
            }

            if (anyDuplicated(x$Marker) ||
                !setequal(x$Marker, SETS)) {
                stop("Missing or duplicated SNP sets: ", accuracy_file)
            }

            if (!all(x$Outer == outer) ||
                !all(x$Inner == inner) ||
                !all(x$Trait == trait)) {
                stop("Fold/trait metadata mismatch: ", accuracy_file)
            }

            if (any(x$N_train != length(train_ids)) ||
                any(x$N_valid != length(valid_ids))) {
                stop("Sample count mismatch: ", accuracy_file)
            }

            if (
                anyNA(x) ||
                any(!is.finite(x$Accuracy)) ||
                any(abs(x$Accuracy) > 1) ||
                any(!is.finite(x$RMSE)) ||
                any(x$RMSE < 0)
            ) {
                stop("Invalid prediction statistics: ", accuracy_file)
            }

            if (
                any(x$N_selected <= 0) ||
                any(x$N_used <= 0) ||
                any(x$N_used > x$N_selected) ||
                any(x$Alleles_flipped < 0) ||
                any(x$Alleles_flipped > x$N_selected)
            ) {
                stop("Invalid SNP counts: ", accuracy_file)
            }

            k <- k + 1L
            all_results[[k]] <- x
        }
    }

    # The three validation folds must partition outer training.
    all_valid <- unlist(valid_groups, use.names = FALSE)

    if (anyDuplicated(all_valid) ||
        !setequal(all_valid, outer_train)) {
        stop("Inner validation folds do not partition outer_", outer)
    }
}

ALL <- do.call(rbind, all_results)
rownames(ALL) <- NULL

if (nrow(ALL) != 1080L) {
    stop("Expected 1080 inner validation results")
}

save_table(
    ALL,
    file.path(OUTDIR, "all_inner_accuracies.tsv")
)

##################################################
# Mean inner-validation accuracy per outer/trait
##################################################

all_summaries <- list()
best_settings <- list()

summary_k <- 0L
best_k <- 0L

for (outer in OUTERS) {

    outer_out <- file.path(
        OUTDIR,
        paste0("outer_", outer)
    )

    dir.create(
        outer_out,
        recursive = TRUE,
        showWarnings = FALSE
    )

    for (trait in TRAITS) {

        summary_rows <- vector("list", length(SETS))

        for (j in seq_along(SETS)) {

            marker <- SETS[j]

            x <- ALL[
                ALL$Outer == outer &
                ALL$Trait == trait &
                ALL$Marker == marker,
                ,
                drop = FALSE
            ]

            # Require exactly three distinct inner folds.
            if (
                nrow(x) != 3L ||
                !setequal(x$Inner, INNERS)
            ) {
                stop(
                    "Incomplete inner results: outer_",
                    outer, " ", trait, " ", marker
                )
            }

            x <- x[match(INNERS, x$Inner), , drop = FALSE]

            summary_rows[[j]] <- data.frame(
                Outer = outer,
                Trait = trait,
                Marker = marker,
                N_inner = 3L,
                Mean_Accuracy = mean(x$Accuracy),
                SD_Accuracy = sd(x$Accuracy),
                Mean_RMSE = mean(x$RMSE),
                Mean_N_selected = mean(x$N_selected),
                Min_N_selected = min(x$N_selected),
                Max_N_selected = max(x$N_selected),
                Mean_N_used = mean(x$N_used),
                stringsAsFactors = FALSE
            )
        }

        summary_df <- do.call(rbind, summary_rows)

        # Rank by the mean accuracy across all three inners.
        # Ties follow SETS order (smaller numerical sets first).
        ord <- order(
            -summary_df$Mean_Accuracy,
            match(summary_df$Marker, SETS)
        )

        summary_df <- summary_df[ord, , drop = FALSE]
        rownames(summary_df) <- NULL

        summary_df$Rank <- seq_len(nrow(summary_df))

        summary_df <- summary_df[
            ,
            c("Rank", setdiff(names(summary_df), "Rank"))
        ]

        save_table(
            summary_df,
            file.path(
                outer_out,
                paste0(trait, "_chip_mean_accuracies.tsv")
            )
        )

        best <- summary_df[1, , drop = FALSE]

        # Convenient single-line input for subsequent Step 8.
        writeLines(
            best$Marker,
            file.path(
                outer_out,
                paste0(trait, "_best_marker.txt")
            )
        )

        best_k <- best_k + 1L

        best_settings[[best_k]] <- data.frame(
            Outer = outer,
            Trait = trait,
            Best_Marker = best$Marker,
            Mean_Accuracy = best$Mean_Accuracy,
            SD_Accuracy = best$SD_Accuracy,
            Mean_RMSE = best$Mean_RMSE,
            Mean_N_selected = best$Mean_N_selected,
            Min_N_selected = best$Min_N_selected,
            Max_N_selected = best$Max_N_selected,
            stringsAsFactors = FALSE
        )

        summary_k <- summary_k + 1L
        all_summaries[[summary_k]] <- summary_df

        cat(
            "Outer", outer,
            "Trait", trait,
            "Best", best$Marker,
            "Mean accuracy", round(best$Mean_Accuracy, 6),
            "\n"
        )
    }
}

BEST <- do.call(rbind, best_settings)
SUMMARY <- do.call(rbind, all_summaries)

rownames(BEST) <- NULL
rownames(SUMMARY) <- NULL

stopifnot(
    nrow(BEST) == 20L,
    nrow(SUMMARY) == 360L
)

save_table(
    BEST,
    file.path(OUTDIR, "best_chip_settings.tsv")
)

save_table(
    SUMMARY,
    file.path(OUTDIR, "all_chip_mean_accuracies.tsv")
)

cat("\n====================================\n")
cat("All inner results:", nrow(ALL), "\n")
cat("All averaged settings:", nrow(SUMMARY), "\n")
cat("Best outer-trait settings:", nrow(BEST), "\n")
cat("DONE Step7.3\n")
