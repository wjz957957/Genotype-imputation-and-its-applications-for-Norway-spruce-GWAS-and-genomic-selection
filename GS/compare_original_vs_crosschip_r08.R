#!/usr/bin/env Rscript
# Comparison of the original CHIP+QTN and cross-CHIP-LD-filtered CHIP+QTN
# on the SAME 10 outer-test folds. Base R only; never reruns any model.
options(stringsAsFactors = FALSE, scipen = 8)

old_root <- "results/rrBLUP/comparison"
new_root <- "results_crosschip_r2_0p8/rrBLUP/comparison"
out_root <- "results_crosschip_r2_0p8/rrBLUP/comparison/direct_old_vs_crosschip"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) out_root <- args[[1L]]
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

read_trait <- function(trait, root) {
  path <- file.path(root, paste0(trait, "_CHIP_vs_CHIPplusQTN_per_fold.tsv"))
  if (!file.exists(path)) stop("Missing input: ", path)
  d <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE, colClasses = c(outer = "character"))
  required <- c("trait", "outer", "n_markers_CHIP", "n_markers_COMBINED")
  for (metric in c("Pearson", "Spearman", "RMSE", "MAE", "Bias")) {
    required <- c(required, paste0(metric, "_CHIP"), paste0(metric, "_COMBINED"))
  }
  miss <- setdiff(required, names(d))
  if (length(miss)) stop(path, " missing columns: ", paste(miss, collapse = ", "))
  if (any(is.na(suppressWarnings(as.integer(d$outer))))) {
    stop("outer must be numeric (01-10) in ", path)
  }
  d$outer <- sprintf("%02d", as.integer(d$outer))
  if (nrow(d) != 10L || anyDuplicated(d$outer) ||
      !setequal(d$outer, sprintf("%02d", 1:10))) {
    stop("Expected one entry for each outer_01-outer_10 in ", path)
  }
  if (any(d$trait != trait)) stop("Trait mismatch in ", path)
  for (cn in setdiff(required, c("trait", "outer"))) {
    d[[cn]] <- suppressWarnings(as.numeric(d[[cn]]))
    if (any(!is.finite(d[[cn]]))) stop("Invalid numeric values in ", path, " column ", cn)
  }
  d[order(d$outer), , drop = FALSE]
}

out_perfold <- list()
out_summary <- list()
out_counts <- list()
out_baseline <- list()
metric_names <- c("Pearson", "Spearman", "RMSE", "MAE", "Bias", "AbsBias")

exact_signflip_p <- function(delta) {
  # Two-sided paired sign-flip reference distribution (2^10 possibilities).
  n <- length(delta)
  if (all(delta == 0)) return(1.0)
  obs <- abs(mean(delta))
  sums <- numeric(2^n)
  for (mask in 0:(2^n - 1L)) {
    signs <- ifelse(bitwAnd(mask, bitwShiftL(1L, 0:(n-1L))) != 0L, 1.0, -1.0)
    sums[mask + 1L] <- abs(mean(signs * delta))
  }
  mean(sums + 1e-12 >= obs)
}

for (trait in c("BB", "H6")) {
  old <- read_trait(trait, old_root)
  new <- read_trait(trait, new_root)
  if (!identical(old$outer, new$outer)) stop("Outer fold pairing differs for ", trait)

  # A hard guard against accidentally comparing different fold splits, traits,
  # or a different CHIP-only baseline.
  baseline_columns <- c("n_markers_CHIP", paste0(c("Pearson", "Spearman", "RMSE", "MAE", "Bias"), "_CHIP"))
  diffs <- sapply(baseline_columns, function(cn) max(abs(old[[cn]] - new[[cn]])))
  base_ok <- all(diffs < 1e-9)
  out_baseline[[trait]] <- data.frame(
    Trait = trait, Check = baseline_columns, MaxAbsDifference = as.numeric(diffs),
    Pass = diffs < 1e-9, row.names = NULL
  )
  if (!base_ok) {
    bad <- paste(baseline_columns[diffs >= 1e-9], collapse = ", ")
    stop("STOP: CHIP baseline differs between old/new for ", trait, ": ", bad)
  }

  counts <- data.frame(
    Trait = trait, Outer = old$outer,
    CHIP_used = old$n_markers_CHIP,
    Old_QTN_used = old$n_markers_COMBINED - old$n_markers_CHIP,
    New_QTN_used = new$n_markers_COMBINED - new$n_markers_CHIP,
    Old_TOTAL = old$n_markers_COMBINED,
    New_TOTAL = new$n_markers_COMBINED
  )
  if (any(counts$Old_QTN_used < 0 | counts$New_QTN_used < 0)) stop("Negative QTN count for ", trait)
  out_counts[[trait]] <- counts

  for (metric in metric_names) {
    if (metric == "AbsBias") {
      v_old <- abs(old$Bias_COMBINED)
      v_new <- abs(new$Bias_COMBINED)
    } else {
      v_old <- old[[paste0(metric, "_COMBINED")]]
      v_new <- new[[paste0(metric, "_COMBINED")]]
    }
    delta <- v_new - v_old
    # Bigger = better for correlations; smaller = better for errors/absolute bias.
    direction <- if (metric %in% c("Pearson", "Spearman")) 1 else if (metric %in% c("RMSE", "MAE", "AbsBias")) -1 else 0
    win <- if (direction == 0) rep(FALSE, length(delta)) else direction * delta > 1e-12
    loss <- if (direction == 0) rep(FALSE, length(delta)) else direction * delta < -1e-12
    ties <- abs(delta) <= 1e-12
    pf <- data.frame(
      Trait = trait, Outer = old$outer, Metric = metric,
      Original = v_old, CrosschipLD = v_new,
      Delta_new_minus_original = delta,
      CHIP_only = if (metric == "AbsBias") abs(old$Bias_CHIP) else old[[paste0(metric, "_CHIP")]],
      Old_QTN_used = counts$Old_QTN_used, New_QTN_used = counts$New_QTN_used,
      New_Better = if (direction == 0) NA else win
    )
    out_perfold[[paste(trait, metric)]] <- pf

    # Correlation of outer folds due to shared training sets means p values/CI
    # should be treated as descriptive/approximate, not new external validation.
    if (all(abs(delta - delta[[1L]]) < 1e-14)) {
      if (abs(delta[[1L]]) < 1e-14) {
        t_p <- 1; ci <- c(0, 0)
      } else {
        t_p <- NA_real_; ci <- c(NA_real_, NA_real_)
      }
    } else {
      tt <- t.test(v_new, v_old, paired = TRUE, conf.level = 0.95)
      t_p <- tt$p.value; ci <- unname(tt$conf.int)
    }
    if (all(abs(delta) < 1e-14)) {
      wilcox_p <- 1
    } else {
      wilcox_p <- suppressWarnings(wilcox.test(v_new, v_old,
                        paired = TRUE, exact = FALSE, correct = FALSE)$p.value)
    }
    out_summary[[paste(trait, metric)]] <- data.frame(
      Trait = trait, Metric = metric, N_folds = length(delta),
      Original_mean = mean(v_old), CrosschipLD_mean = mean(v_new),
      Mean_Delta_new_minus_original = mean(delta),
      SD_Delta = sd(delta), Median_Delta = median(delta),
      CI95_low = ci[1], CI95_high = ci[2],
      Paired_t_p = t_p, Paired_Wilcoxon_p = wilcox_p,
      Exact_signflip_p = exact_signflip_p(delta),
      New_Better_folds = if (direction == 0) NA_integer_ else sum(win),
      New_Worse_folds = if (direction == 0) NA_integer_ else sum(loss),
      Equal_folds = sum(ties)
    )
  }
}

write_tsv <- function(x, name) {
  path <- file.path(out_root, name)
  write.table(x, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "NA")
  message("[WROTE] ", path, " (", nrow(x), " rows)")
}
summary_df <- do.call(rbind, out_summary)
perfold_df <- do.call(rbind, out_perfold)
counts_df <- do.call(rbind, out_counts)
baseline_df <- do.call(rbind, out_baseline)
write_tsv(summary_df, "old_vs_crosschip_summary.tsv")
write_tsv(perfold_df, "old_vs_crosschip_per_fold.tsv")
write_tsv(counts_df, "old_vs_crosschip_marker_counts.tsv")
write_tsv(baseline_df, "old_vs_crosschip_CHIP_baseline_check.tsv")

cat("\n========== PEARSON and RMSE: DIRECT OLD vs NEW =========\n")
cols <- c("Trait", "Metric", "Original_mean", "CrosschipLD_mean",
          "Mean_Delta_new_minus_original", "CI95_low", "CI95_high",
          "Paired_t_p", "Paired_Wilcoxon_p", "Exact_signflip_p",
          "New_Better_folds", "New_Worse_folds")
print(summary_df[summary_df$Metric %in% c("Pearson", "RMSE"), cols], row.names = FALSE, digits = 6)
cat("\n========== QTN USED BY FOLD =========\n")
print(counts_df[, c("Trait", "Outer", "Old_QTN_used", "New_QTN_used")], row.names = FALSE)
cat("\n[PASS] Both traits: 10 matched outer folds; CHIP-only baselines identical.\n")
cat("NOTE: old/new differ in both selected QTN and QTN count. Comparison is between\n")
cat("      complete strategies, NOT an isolated causal effect of LD filtering.\n")
cat("NOTE: tests over 10 overlapping CV folds are approximate, not external validation.\n")
