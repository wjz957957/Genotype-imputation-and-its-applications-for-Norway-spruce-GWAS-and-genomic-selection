#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
    stop("Usage: Rscript 13c_BB01_crossLD08_summarize.R WORK")
}

WORK <- args[1]

BASE <- file.path(WORK, "09_chip_plus_imputed_RRBLUP_crosschip_r2_0p8_pilot")
OUTDIR <- file.path(BASE, "optimization")

CV_ROOT <- file.path(WORK, "02_CV")

CHIP_ROOT <- file.path(WORK, "07_chip_RRBLUP")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUTERS <- "01"
INNERS <- sprintf("%02d", 1:3)
TRAITS <- "BB"

COUNTS <- as.integer(
    c(0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 75,
      100, 150, 200, 300, 400, 500, 750, 1000)
)

EXPECTED <- c(
    "Outer",
    "Inner",
    "Trait",
    "Best_Chip",
    "N_chip_selected",
    "N_chip_used",
    "Extra_Count",
    "N_extra_used",
    "N_total_selected",
    "N_total_used",
    "N_train",
    "N_valid",
    "Chip_alleles_flipped",
    "Extra_alleles_flipped",
    "Accuracy",
    "RMSE",
    "Delta_Accuracy_vs_0",
    "Delta_RMSE_vs_0"
)

############################################
# Helpers
############################################

read_ids <- function(path) {

    if (!file.exists(path)) {
        stop("Missing sample file: ", path)
    }

    x <- trimws(readLines(path, warn = FALSE))
    x <- x[nzchar(x)]

    if (length(x) == 0L || anyDuplicated(x)) {
        stop("Invalid sample IDs: ", path)
    }

    x
}

save_table <- function(x, path) {

    temporary <- paste0(path, ".tmp.", Sys.getpid())

    write.table(
        x,
        file = temporary,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )

    if (!file.rename(temporary, path)) {
        stop("Cannot publish file: ", path)
    }
}

############################################
# 1. Load best chip settings from Step 7
############################################

BEST_CHIP_FILE <- file.path(
    CHIP_ROOT, "optimization",
    "best_chip_settings.tsv"
)

if (!file.exists(BEST_CHIP_FILE)) {
    stop("Missing best chip settings")
}

BEST_CHIP <- read.delim(
    BEST_CHIP_FILE,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(
        "character", "character", "character",
        rep("numeric", 6)
    )
)

if (!identical(
    names(BEST_CHIP)[1:3],
    c("Outer", "Trait", "Best_Marker")
)) {
    stop("Unexpected Step 7 best chip columns")
}

if (nrow(BEST_CHIP) != 20L ||
    anyDuplicated(
        paste(BEST_CHIP$Outer, BEST_CHIP$Trait)
    )) {
    stop("Invalid Step 7 best chip settings")
}

for (outer in OUTERS) {
    for (trait in TRAITS) {

        z <- BEST_CHIP[
            BEST_CHIP$Outer == outer &
            BEST_CHIP$Trait == trait,
            ,
            drop = FALSE
        ]

        if (nrow(z) != 1L ||
            !startsWith(z$Best_Marker, "chip_")) {
            stop("Invalid chip configuration: ", outer, " ", trait)
        }
    }
}

############################################
# 2. Load and validate BB outer_01 in three inner folds
############################################

all_results <- vector("list", 3L)

k <- 0L

for (outer in OUTERS) {

    outer_dir <- file.path(
        CV_ROOT, paste0("outer_", outer)
    )

    outer_train <- read_ids(
        file.path(outer_dir, "outer_train.samples")
    )

    outer_test <- read_ids(
        file.path(outer_dir, "outer_test.samples")
    )

    valid_groups <- list()

    for (inner in INNERS) {

        fold_dir <- file.path(
            outer_dir, paste0("inner_", inner)
        )

        train_ids <- read_ids(
            file.path(fold_dir, "inner_train.samples")
        )

        valid_ids <- read_ids(
            file.path(fold_dir, "inner_valid.samples")
        )

        if (!setequal(
            train_ids,
            setdiff(outer_train, valid_ids)
        )) {
            stop("Incorrect inner partition: ", outer, " ", inner)
        }

        if (length(intersect(train_ids, valid_ids)) ||
            length(intersect(train_ids, outer_test)) ||
            length(intersect(valid_ids, outer_test))) {
            stop("Nested CV sample leakage: ", outer, " ", inner)
        }

        valid_groups[[inner]] <- valid_ids

        input_dir <- file.path(
            BASE,
            paste0("outer_", outer),
            paste0("inner_", inner)
        )

        if (!file.exists(
            file.path(input_dir, "step9b.complete.flag")
        )) {
            stop("Step 9b incomplete: ", input_dir)
        }

        for (trait in TRAITS) {

            input <- file.path(
                input_dir,
                paste0(
                    trait,
                    "_chip_plus_imputed_accuracy.tsv"
                )
            )

            if (!file.exists(input)) {
                stop("Missing accuracy file: ", input)
            }

            x <- read.delim(
                input,
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE,
                check.names = FALSE,
                colClasses = c(
                    rep("character", 4),
                    rep("integer", 10),
                    rep("numeric", 4)
                )
            )

            if (!identical(names(x), EXPECTED)) {
                stop("Unexpected columns: ", input)
            }

            if (nrow(x) != 20L ||
                !identical(x$Extra_Count, COUNTS)) {
                stop("Incorrect extra SNP settings: ", input)
            }

            if (!all(x$Outer == outer) ||
                !all(x$Inner == inner) ||
                !all(x$Trait == trait)) {
                stop("Fold/trait mismatch: ", input)
            }

            best_row <- BEST_CHIP[
                BEST_CHIP$Outer == outer &
                BEST_CHIP$Trait == trait,
                ,
                drop = FALSE
            ]

            best_marker <- best_row$Best_Marker

            if (!all(x$Best_Chip == best_marker)) {
                stop("Best chip mismatch: ", input)
            }

            if (anyNA(x) ||
                any(!is.finite(x$Accuracy)) ||
                any(abs(x$Accuracy) > 1) ||
                any(!is.finite(x$RMSE)) ||
                any(x$RMSE < 0)) {
                stop("Invalid RRBLUP statistics: ", input)
            }

            if (any(x$N_train != length(train_ids)) ||
                any(x$N_valid != length(valid_ids))) {
                stop("Sample count mismatch: ", input)
            }

            if (length(unique(x$N_chip_selected)) != 1L ||
                length(unique(x$N_chip_used)) != 1L ||
                length(unique(x$Chip_alleles_flipped)) != 1L) {
                stop("Chip baseline varies across extra settings: ", input)
            }

            if (any(
                x$N_total_selected !=
                    x$N_chip_selected + x$Extra_Count
            )) {
                stop("Total selected SNP count mismatch: ", input)
            }

            if (any(
                x$N_total_used !=
                    x$N_chip_used + x$N_extra_used
            )) {
                stop("Total used SNP count mismatch: ", input)
            }

            if (any(x$N_chip_used < 1) ||
                any(x$N_chip_used > x$N_chip_selected) ||
                any(x$N_extra_used < 0) ||
                any(x$N_extra_used > x$Extra_Count)) {
                stop("Invalid SNP count: ", input)
            }

            if (x$N_extra_used[1] != 0L ||
                x$Extra_alleles_flipped[1] != 0L) {
                stop("Invalid extra_0 baseline: ", input)
            }

            baseline_acc <- x$Accuracy[1]
            baseline_rmse <- x$RMSE[1]

            if (any(abs(
                x$Delta_Accuracy_vs_0 -
                (x$Accuracy - baseline_acc)
            ) > 2e-5)) {
                stop("Accuracy delta mismatch: ", input)
            }

            if (any(abs(
                x$Delta_RMSE_vs_0 -
                (x$RMSE - baseline_rmse)
            ) > 2e-5)) {
                stop("RMSE delta mismatch: ", input)
            }

            ################################
            # Recheck against Step 7
            ################################

            old_file <- file.path(
                CHIP_ROOT,
                paste0("outer_", outer),
                paste0("inner_", inner),
                paste0(trait, "_chip_accuracy.txt")
            )

            old <- read.delim(
                old_file,
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE,
                colClasses = c(
                    rep("character", 4),
                    rep("integer", 5),
                    rep("numeric", 2)
                )
            )

            old_row <- old[
                old$Marker == best_marker,
                ,
                drop = FALSE
            ]

            if (nrow(old_row) != 1L) {
                stop("Cannot find Step 7 baseline: ", old_file)
            }

            if (x$N_chip_selected[1] != old_row$N_selected ||
                x$N_chip_used[1] != old_row$N_used ||
                x$Chip_alleles_flipped[1] != old_row$Alleles_flipped ||
                abs(baseline_acc - old_row$Accuracy) > 2e-5 ||
                abs(baseline_rmse - old_row$RMSE) >
                    2e-5 * max(1, abs(old_row$RMSE))) {
                stop("Step 9 baseline differs from Step 7: ", input)
            }

            k <- k + 1L
            all_results[[k]] <- x

            cat(
                "PASS:", outer, inner, trait,
                "chip=", best_marker,
                "baseline=", round(baseline_acc, 6),
                "\n"
            )
        }
    }

    all_valid <- unlist(
        valid_groups, use.names = FALSE
    )

    if (anyDuplicated(all_valid) ||
        !setequal(all_valid, outer_train)) {
        stop("Inner validation folds do not partition outer train: ", outer)
    }
}

ALL <- do.call(rbind, all_results)
rownames(ALL) <- NULL

if (k != 3L || nrow(ALL) != 60L) {
    stop("Expected exactly 60 validated BB outer_01 RRBLUP models")
}

save_table(
    ALL,
    file.path(OUTDIR, "all_inner_accuracies.tsv")
)

############################################
# 3. Mean accuracy across three inner folds
############################################

all_means <- list()
best_rows <- list()

m <- 0L
b <- 0L

for (outer in OUTERS) {

    outer_out <- file.path(
        OUTDIR, paste0("outer_", outer)
    )

    dir.create(
        outer_out, recursive = TRUE,
        showWarnings = FALSE
    )

    for (trait in TRAITS) {

        rows <- vector("list", length(COUNTS))

        for (j in seq_along(COUNTS)) {

            n_extra <- COUNTS[j]

            x <- ALL[
                ALL$Outer == outer &
                ALL$Trait == trait &
                ALL$Extra_Count == n_extra,
                ,
                drop = FALSE
            ]

            if (nrow(x) != 3L ||
                !setequal(x$Inner, INNERS)) {
                stop(
                    "Incomplete inner results: ",
                    outer, " ", trait, " ", n_extra
                )
            }

            x <- x[
                match(INNERS, x$Inner),
                ,
                drop = FALSE
            ]

            rows[[j]] <- data.frame(
                Outer = outer,
                Trait = trait,
                Best_Chip = x$Best_Chip[1],
                Extra_Count = n_extra,
                N_inner = 3L,
                Mean_Accuracy = mean(x$Accuracy),
                SD_Accuracy = sd(x$Accuracy),
                Mean_RMSE = mean(x$RMSE),
                Mean_Delta_Accuracy_vs_0 =
                    mean(x$Delta_Accuracy_vs_0),
                SD_Delta_Accuracy_vs_0 =
                    sd(x$Delta_Accuracy_vs_0),
                Mean_Delta_RMSE_vs_0 =
                    mean(x$Delta_RMSE_vs_0),
                Mean_N_chip_selected =
                    mean(x$N_chip_selected),
                Mean_N_chip_used =
                    mean(x$N_chip_used),
                Mean_N_extra_used =
                    mean(x$N_extra_used),
                stringsAsFactors = FALSE
            )
        }

        summary_df <- do.call(rbind, rows)

        ####################################
        # Rank by mean validation accuracy.
        # Near-exact ties favor fewer SNPs.
        ####################################

        max_acc <- max(summary_df$Mean_Accuracy)

        tied <- which(
            summary_df$Mean_Accuracy >= max_acc - 1e-8
        )

        best_idx <- tied[
            which.min(summary_df$Extra_Count[tied])
        ]

        best <- summary_df[
            best_idx, , drop = FALSE
        ]

        baseline <- summary_df[
            summary_df$Extra_Count == 0L,
            ,
            drop = FALSE
        ]

        if (nrow(baseline) != 1L) {
            stop("Missing extra_0 baseline")
        }

        ranked <- summary_df[
            order(
                -summary_df$Mean_Accuracy,
                summary_df$Extra_Count
            ),
            ,
            drop = FALSE
        ]

        ranked$Rank <- seq_len(nrow(ranked))

        ranked <- ranked[
            ,
            c("Rank", setdiff(names(ranked), "Rank"))
        ]

        rownames(ranked) <- NULL

        save_table(
            ranked,
            file.path(
                outer_out,
                paste0(
                    trait,
                    "_extra_mean_accuracies.tsv"
                )
            )
        )

        writeLines(
            as.character(best$Extra_Count),
            file.path(
                outer_out,
                paste0(trait, "_best_extra.txt")
            )
        )

        b <- b + 1L

        best_rows[[b]] <- data.frame(
            Outer = outer,
            Trait = trait,
            Best_Chip = best$Best_Chip,
            Best_Extra = best$Extra_Count,
            Baseline_Mean_Accuracy =
                baseline$Mean_Accuracy,
            Best_Mean_Accuracy =
                best$Mean_Accuracy,
            Gain_Accuracy =
                best$Mean_Accuracy -
                baseline$Mean_Accuracy,
            Baseline_Mean_RMSE =
                baseline$Mean_RMSE,
            Best_Mean_RMSE =
                best$Mean_RMSE,
            Mean_N_chip_selected =
                best$Mean_N_chip_selected,
            Mean_N_extra_used =
                best$Mean_N_extra_used,
            stringsAsFactors = FALSE
        )

        m <- m + 1L
        all_means[[m]] <- summary_df

        cat(
            "SELECTED:",
            "outer=", outer,
            "trait=", trait,
            "chip=", best$Best_Chip,
            "extra=", best$Extra_Count,
            "baseline=", round(
                baseline$Mean_Accuracy, 6
            ),
            "best=", round(
                best$Mean_Accuracy, 6
            ),
            "gain=", round(
                best$Mean_Accuracy -
                baseline$Mean_Accuracy, 6
            ),
            "\n"
        )
    }
}

MEANS <- do.call(rbind, all_means)
BEST <- do.call(rbind, best_rows)

rownames(MEANS) <- NULL
rownames(BEST) <- NULL

if (nrow(MEANS) != 20L ||
    nrow(BEST) != 1L) {
    stop("Incorrect optimization result counts")
}

save_table(
    MEANS,
    file.path(
        OUTDIR, "all_extra_mean_accuracies.tsv"
    )
)

save_table(
    BEST,
    file.path(
        OUTDIR, "best_extra_settings.tsv"
    )
)

############################################
# 4. Final checks
############################################

if (anyNA(BEST) ||
    any(!BEST$Best_Extra %in% COUNTS) ||
    any(BEST$Gain_Accuracy < -1e-8)) {
    stop("Invalid best extra SNP settings")
}

cat("\n========================================\n")
cat("Step 9c completed\n")
cat("Validated inner models:", nrow(ALL), "\n")
cat("Mean accuracy settings:", nrow(MEANS), "\n")
cat("Best outer-trait settings:", nrow(BEST), "\n")
cat("========================================\n")

print(BEST[
    ,
    c(
        "Outer", "Trait", "Best_Chip",
        "Best_Extra",
        "Baseline_Mean_Accuracy",
        "Best_Mean_Accuracy",
        "Gain_Accuracy"
    )
])

cat("\nDONE Step9c\n")
