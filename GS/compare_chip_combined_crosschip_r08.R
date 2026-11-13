options(stringsAsFactors = FALSE)

traits <- c("BB", "H6")
metrics <- c("Pearson", "Spearman", "RMSE", "MAE", "Bias")

all_results <- list()

for (trait in traits) {

    chip_list <- list()
    combined_list <- list()

    for (o in sprintf("%02d", 1:10)) {

        # CHIP baseline is reused from the original analysis after byte-for-byte
        # verification of all CHIP train/test RAW files.
        chip_file <- file.path(
            "results", "rrBLUP", trait, "chip",
            paste0("outer_", o, "_metrics.tsv")
        )

        # Cross-CHIP LD r2=0.8 combined results are written separately.
        combined_file <- file.path(
            "results_crosschip_r2_0p8", "rrBLUP", trait, "combined",
            paste0("outer_", o, "_metrics.tsv")
        )

        if (!file.exists(chip_file))
            stop("Missing CHIP file: ", chip_file)

        if (!file.exists(combined_file))
            stop("Missing CHIP+QTN file: ", combined_file)

        chip <- read.table(
            chip_file,
            header = TRUE,
            sep = "\t",
            check.names = FALSE
        )

        combined <- read.table(
            combined_file,
            header = TRUE,
            sep = "\t",
            check.names = FALSE
        )

        chip_list[[o]] <- chip
        combined_list[[o]] <- combined
    }

    chip_all <- do.call(rbind, chip_list)
    combined_all <- do.call(rbind, combined_list)

    rownames(chip_all) <- NULL
    rownames(combined_all) <- NULL

    chip_all$outer <- sprintf("%02d", as.integer(chip_all$outer))
    combined_all$outer <- sprintf("%02d", as.integer(combined_all$outer))

    m <- merge(
        chip_all[, c("trait", "outer", "n_markers", metrics)],
        combined_all[, c("trait", "outer", "n_markers", metrics)],
        by = c("trait", "outer"),
        suffixes = c("_CHIP", "_COMBINED"),
        sort = TRUE
    )

    ## calculate paired differences
    for (metric in metrics) {
        m[[paste0("Delta_", metric)]] <-
            m[[paste0(metric, "_COMBINED")]] -
            m[[paste0(metric, "_CHIP")]]
    }

    ## save per-fold results
    outdir <- file.path(
        "results_crosschip_r2_0p8", "rrBLUP", "comparison"
    )

    dir.create(
        outdir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    write.table(
        m,
        file.path(
            outdir,
            paste0(trait, "_CHIP_vs_CHIPplusQTN_per_fold.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ## summary
    summary_list <- list()

    for (metric in metrics) {

        chip_x <- m[[paste0(metric, "_CHIP")]]
        comb_x <- m[[paste0(metric, "_COMBINED")]]
        delta <- m[[paste0("Delta_", metric)]]

        ## paired t-test
        ttest <- t.test(
            comb_x,
            chip_x,
            paired = TRUE
        )

        ## Wilcoxon signed-rank test
        ## suppress warning when all differences are zero
        wilcox_p <- tryCatch(
            {
                suppressWarnings(
                    wilcox.test(
                        comb_x,
                        chip_x,
                        paired = TRUE,
                        exact = FALSE
                    )$p.value
                )
            },
            error = function(e) NA_real_
        )

        summary_list[[metric]] <- data.frame(

            trait = trait,
            metric = metric,

            CHIP_mean = mean(chip_x),
            CHIP_SD = sd(chip_x),

            CHIPplusQTN_mean = mean(comb_x),
            CHIPplusQTN_SD = sd(comb_x),

            Mean_Delta = mean(delta),
            SD_Delta = sd(delta),

            CI95_low = ttest$conf.int[1],
            CI95_high = ttest$conf.int[2],

            paired_t_p = ttest$p.value,
            paired_wilcoxon_p = wilcox_p,

            stringsAsFactors = FALSE
        )
    }

    summary <- do.call(rbind, summary_list)

    write.table(
        summary,
        file.path(
            outdir,
            paste0(trait, "_CHIP_vs_CHIPplusQTN_summary.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    all_results[[trait]] <- summary
}

cat("\n===== CHIP vs CHIP+QTN SUMMARY =====\n\n")

for (trait in traits) {

    cat("\n########## ", trait, " ##########\n", sep = "")

    x <- all_results[[trait]]

    print(
        x[, c(
            "metric",
            "CHIP_mean",
            "CHIP_SD",
            "CHIPplusQTN_mean",
            "CHIPplusQTN_SD",
            "Mean_Delta",
            "paired_t_p",
            "paired_wilcoxon_p"
        )],
        digits = 6,
        row.names = FALSE
    )
}

cat("\n===== DONE =====\n")
cat("Per-fold results:\n")
cat("results_crosschip_r2_0p8/rrBLUP/comparison/*_CHIP_vs_CHIPplusQTN_per_fold.tsv\n")
cat("\nSummary results:\n")
cat("results_crosschip_r2_0p8/rrBLUP/comparison/*_CHIP_vs_CHIPplusQTN_summary.tsv\n")
