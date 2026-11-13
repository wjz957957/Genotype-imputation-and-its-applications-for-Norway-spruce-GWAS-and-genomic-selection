metrics <- c("Pearson", "Spearman", "RMSE", "MAE", "Bias")
traits <- c("BB", "H6")
genotypes <- c("chip", "combined")
folds <- sprintf("%02d", 1:10)

chip_base_dir <- "deepgs_results"
combined_base_dir <- "deepgs_results_crosschip_r2_0p8"
out_dir <- "deepgs_comparison_crosschip_r2_0p8"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

get_col <- function(d, name, required = TRUE) {
    idx <- which(tolower(names(d)) == tolower(name))
    if (length(idx) != 1) {
        if (!required) return(NA_real_)
        stop(
            "Column '", name, "' not found or ambiguous. Available: ",
            paste(names(d), collapse = ", ")
        )
    }
    d[[idx]]
}

# Read and validate all 40 metric files
records <- list()
k <- 0L

for (trait in traits) {
    for (geno in genotypes) {
        for (fold in folds) {

            d <- if (geno == "chip") {
                file.path(chip_base_dir, trait, geno)
            } else {
                file.path(combined_base_dir, trait, geno)
            }
            prefix <- file.path(d, paste0("outer_", fold))

            metric_file <- paste0(prefix, "_metrics.tsv")
            pred_file <- paste0(prefix, "_predictions.tsv")
            complete_file <- paste0(prefix, ".complete")

            if (!file.exists(complete_file))
                stop("Missing completion marker: ", complete_file)

            if (!file.exists(metric_file) ||
                file.info(metric_file)$size == 0)
                stop("Missing/empty metrics: ", metric_file)

            if (!file.exists(pred_file) ||
                file.info(pred_file)$size == 0)
                stop("Missing/empty predictions: ", pred_file)

            x <- read.delim(
                metric_file,
                check.names = FALSE,
                stringsAsFactors = FALSE
            )

            if (nrow(x) != 1)
                stop("Expected exactly one metrics row: ", metric_file)

            row <- data.frame(
                trait = trait,
                genotype = geno,
                outer = fold,
                n_markers = as.numeric(
                    get_col(x, "n_markers", required = FALSE)[1]
                ),
                stringsAsFactors = FALSE
            )

            for (m in metrics) {
                row[[m]] <- as.numeric(get_col(x, m)[1])
            }

            if (any(!is.finite(as.matrix(row[, metrics]))))
                stop("Invalid metric value: ", metric_file)

            k <- k + 1L
            records[[k]] <- row
        }
    }
}

all_metrics <- do.call(rbind, records)

if (nrow(all_metrics) != 40)
    stop("Expected 40 metric records, found ", nrow(all_metrics))

write.table(
    all_metrics,
    file.path(out_dir, "DEEPGS_all_folds.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

# Paired comparison within each outer fold
fold_records <- list()
k <- 0L

for (trait in traits) {
    for (fold in folds) {


        # Use explicit indexing to avoid relying on non-base R syntax
        chip <- all_metrics[
            all_metrics$trait == trait &
            all_metrics$genotype == "chip" &
            all_metrics$outer == fold, , drop = FALSE
        ]
        combined <- all_metrics[
            all_metrics$trait == trait &
            all_metrics$genotype == "combined" &
            all_metrics$outer == fold, , drop = FALSE
        ]

        if (nrow(chip) != 1 || nrow(combined) != 1)
            stop("Missing/duplicated model record: ", trait, " outer_", fold)

        row <- data.frame(
            trait = trait,
            outer = fold,
            n_markers_CHIP = chip$n_markers,
            n_markers_CHIPplusQTN = combined$n_markers
        )

        for (m in metrics) {
            row[[paste0(m, "_CHIP")]] <- chip[[m]]
            row[[paste0(m, "_CHIPplusQTN")]] <- combined[[m]]
            row[[paste0("Delta_", m)]] <- combined[[m]] - chip[[m]]
        }

        k <- k + 1L
        fold_records[[k]] <- row
    }
}

per_fold <- do.call(rbind, fold_records)

write.table(
    per_fold,
    file.path(out_dir, "DEEPGS_CHIP_vs_CHIPplusQTN_per_fold.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

# Summary across the 10 outer folds
safe_p <- function(x, y, method) {
    if (all(abs(x - y) < 1e-12)) return(NA_real_)

    tryCatch({
        if (method == "t") {
            t.test(x, y, paired = TRUE)$p.value
        } else {
            wilcox.test(
                x, y, paired = TRUE, exact = FALSE
            )$p.value
        }
    }, error = function(e) NA_real_)
}

summary_records <- list()
k <- 0L
tol <- 1e-12

for (trait in traits) {
    z <- per_fold[per_fold$trait == trait, , drop = FALSE]

    for (m in metrics) {
        chip <- z[[paste0(m, "_CHIP")]]
        combined <- z[[paste0(m, "_CHIPplusQTN")]]
        delta <- combined - chip

        # Improvement direction depends on the metric
        if (m %in% c("Pearson", "Spearman")) {
            diff <- delta
        } else if (m %in% c("RMSE", "MAE")) {
            diff <- -delta
        } else {
            # Bias is better when its absolute value is closer to zero
            diff <- abs(chip) - abs(combined)
        }

        improved <- sum(diff > tol)
        worse <- sum(diff < -tol)
        tied <- sum(abs(diff) <= tol)

        k <- k + 1L
        summary_records[[k]] <- data.frame(
            trait = trait,
            metric = m,
            CHIP_mean = mean(chip),
            CHIP_SD = sd(chip),
            CHIPplusQTN_mean = mean(combined),
            CHIPplusQTN_SD = sd(combined),
            Mean_Delta = mean(delta),
            SD_Delta = sd(delta),
            folds_improved = improved,
            folds_worse = worse,
            folds_tied = tied,
            paired_t_p = safe_p(combined, chip, "t"),
            paired_wilcoxon_p = safe_p(combined, chip, "wilcox")
        )
    }
}

summary_table <- do.call(rbind, summary_records)

write.table(
    summary_table,
    file.path(out_dir, "DEEPGS_CHIP_vs_CHIPplusQTN_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n===== DEEPGS: CHIP vs CHIP+QTN =====\n")
cat("Metric records validated:", nrow(all_metrics), "\n")
cat("Fold comparisons:", nrow(per_fold), "\n\n")
print(summary_table, row.names = FALSE, digits = 5)

cat("\n===== OUTPUT FILES =====\n")
cat(file.path(out_dir, "DEEPGS_all_folds.tsv"), "\n")
cat(file.path(out_dir, "DEEPGS_CHIP_vs_CHIPplusQTN_per_fold.tsv"), "\n")
cat(file.path(out_dir, "DEEPGS_CHIP_vs_CHIPplusQTN_summary.tsv"), "\n")
