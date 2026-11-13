#!/usr/bin/env Rscript

# ============================================================
# Identify newly detected association regions after imputation
#
# Definition:
#   An imputed association clump is "newly detected" when the
#   ORIGINAL GWAS contains no genome-wide significant SNP
#   within +/-50 kb of the imputed index SNP.
#
# No conditional association analysis is assumed.
#
# Run from the GWAS root directory after:
#   run_revised_gwas_clumping.sh
#
# Required package:
#   data.table
# ============================================================

suppressPackageStartupMessages(library(data.table))

BASE_DIR <- getwd()
WINDOW_BP <- 50000L

cfg <- data.table(
  Trait = c("BB", "DBH", "FD"),
  Subdir = c("GWAS-Bud", "GWAS-DBH", "GWAS-FD"),
  OriginalGWAS = c(
    "Bud-origin-GWAS.loco_fixed.mlma",
    "DBH-origin-GWAS.loco_fixed.mlma",
    "FD-origin-GWAS.loco_fixed.mlma"
  ),
  OriginalClumped = c(
    "gwas_clumps_revised_origin_Bud.clumped",
    "gwas_clumps_revised_origin_DBH.clumped",
    "gwas_clumps_revised_origin_FD.clumped"
  ),
  ImputedClumped = c(
    "gwas_clumps_revised_imputed_Bud.clumped",
    "gwas_clumps_revised_imputed_DBH.clumped",
    "gwas_clumps_revised_imputed_FD.clumped"
  ),
  OriginalThreshold = c(
    9.1044830475e-07,
    7.5243412439e-07,
    6.7133017361e-07
  ),
  ImputedThreshold = c(
    9.6525655554e-08,
    8.8580079049e-08,
    8.5449489952e-08
  )
)

find_col <- function(nms, candidates) {
  low <- tolower(nms)
  for (cand in candidates) {
    idx <- which(low == tolower(cand))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

normalize_chr <- function(x) {
  as.integer(gsub("[^0-9]", "", as.character(x)))
}

read_original_gwas <- function(path) {
  if (!file.exists(path)) stop("Missing original GWAS file: ", path)

  hdr <- fread(path, nrows = 0)
  nms <- names(hdr)

  chr_col <- find_col(nms, c("Chr", "CHR", "chromosome"))
  bp_col  <- find_col(nms, c("bp", "BP", "pos", "POS"))
  p_col   <- find_col(nms, c("p", "P"))
  snp_col <- find_col(nms, c("SNP_ID", "SNP", "rsid", "ID"))

  if (any(is.na(c(chr_col, bp_col, p_col)))) {
    stop(
      "Could not identify Chr/BP/P columns in: ", path,
      "\nColumns: ", paste(nms, collapse = ", ")
    )
  }

  keep <- unique(na.omit(c(chr_col, bp_col, p_col, snp_col)))
  d <- fread(path, select = keep)

  if (chr_col != "CHR") setnames(d, chr_col, "CHR")
  if (bp_col != "BP") setnames(d, bp_col, "BP")
  if (p_col != "P") setnames(d, p_col, "P")
  if (!is.na(snp_col) && snp_col != "SNP_ID") setnames(d, snp_col, "SNP_ID")
  if (!"SNP_ID" %in% names(d)) d[, SNP_ID := NA_character_]

  d[, CHR := normalize_chr(CHR)]
  d[, BP := as.numeric(BP)]
  d[, P := as.numeric(P)]

  d <- d[
    !is.na(CHR) & CHR >= 1 & CHR <= 12 &
    is.finite(BP) &
    is.finite(P) & P > 0 & P <= 1
  ]

  d
}

read_clumped <- function(path) {
  if (!file.exists(path)) {
    message("No .clumped file found: ", path, " (treated as zero clumps)")
    return(data.table(
      CHR = integer(),
      SNP = character(),
      BP = numeric(),
      P = numeric()
    ))
  }

  d <- fread(path, fill = TRUE)

  required <- c("CHR", "SNP", "BP", "P")
  if (!all(required %in% names(d))) {
    stop(
      "Unexpected .clumped format in: ", path,
      "\nColumns: ", paste(names(d), collapse = ", ")
    )
  }

  d <- d[, .(
    CHR = normalize_chr(CHR),
    SNP = as.character(SNP),
    BP = as.numeric(BP),
    P = as.numeric(P)
  )]

  d <- d[
    !is.na(CHR) & CHR >= 1 & CHR <= 12 &
    is.finite(BP) &
    is.finite(P)
  ]

  d
}

classify_imputed_regions <- function(trait_row) {

  trait <- trait_row$Trait
  subdir <- file.path(BASE_DIR, trait_row$Subdir)

  original_gwas_path <- file.path(subdir, trait_row$OriginalGWAS)
  original_clump_path <- file.path(subdir, trait_row$OriginalClumped)
  imputed_clump_path  <- file.path(subdir, trait_row$ImputedClumped)

  message("\n==================== ", trait, " ====================")

  original_gwas <- read_original_gwas(original_gwas_path)
  original_sig <- original_gwas[P <= trait_row$OriginalThreshold]

  original_clumps <- read_clumped(original_clump_path)
  imputed_clumps  <- read_clumped(imputed_clump_path)

  if (nrow(imputed_clumps) == 0) {
    detail <- data.table(
      Trait = character(),
      Imputed_index_SNP = character(),
      CHR = integer(),
      BP = numeric(),
      Imputed_index_P = numeric(),
      Region_start = numeric(),
      Region_end = numeric(),
      Original_sig_SNPs_within_50kb = integer(),
      Original_min_P_within_50kb = numeric(),
      Closest_original_sig_distance_bp = numeric(),
      Newly_detected = logical()
    )
  } else {

    rows <- vector("list", nrow(imputed_clumps))

    for (i in seq_len(nrow(imputed_clumps))) {

      chr_i <- imputed_clumps$CHR[i]
      bp_i  <- imputed_clumps$BP[i]

      start_i <- max(1, bp_i - WINDOW_BP)
      end_i   <- bp_i + WINDOW_BP

      all_window <- original_gwas[
        CHR == chr_i &
        BP >= start_i &
        BP <= end_i
      ]

      sig_window <- original_sig[
        CHR == chr_i &
        BP >= start_i &
        BP <= end_i
      ]

      min_p_window <- if (nrow(all_window) > 0) {
        min(all_window$P, na.rm = TRUE)
      } else {
        NA_real_
      }

      closest_sig_dist <- if (nrow(sig_window) > 0) {
        min(abs(sig_window$BP - bp_i), na.rm = TRUE)
      } else {
        NA_real_
      }

      rows[[i]] <- data.table(
        Trait = trait,
        Imputed_index_SNP = imputed_clumps$SNP[i],
        CHR = chr_i,
        BP = bp_i,
        Imputed_index_P = imputed_clumps$P[i],
        Region_start = start_i,
        Region_end = end_i,
        Original_sig_SNPs_within_50kb = nrow(sig_window),
        Original_min_P_within_50kb = min_p_window,
        Closest_original_sig_distance_bp = closest_sig_dist,
        Newly_detected = (nrow(sig_window) == 0)
      )
    }

    detail <- rbindlist(rows)
  }

  summary <- data.table(
    Trait = trait,
    Original_genomewide_significant_SNPs = nrow(original_sig),
    Original_clumps = nrow(original_clumps),
    Imputed_clumps = nrow(imputed_clumps),
    Newly_detected_imputed_regions = sum(detail$Newly_detected),
    Imputed_regions_already_detected_originally = sum(!detail$Newly_detected),
    Original_threshold = trait_row$OriginalThreshold,
    Imputed_threshold = trait_row$ImputedThreshold,
    Comparison_window = "+/-50 kb"
  )

  fwrite(
    detail,
    file.path(BASE_DIR, paste0(trait, "_imputed_region_classification.tsv")),
    sep = "\t",
    na = "NA"
  )

  list(summary = summary, detail = detail)
}

all_results <- lapply(seq_len(nrow(cfg)), function(i) {
  classify_imputed_regions(cfg[i])
})

summary_table <- rbindlist(lapply(all_results, `[[`, "summary"))

fwrite(
  summary_table,
  file.path(BASE_DIR, "association_region_summary.tsv"),
  sep = "\t",
  na = "NA"
)

cat("\n============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")
print(summary_table)

cat("\nDefinition used:\n")
cat("An imputed clump is 'newly detected' when the corresponding\n")
cat("original GWAS has no genome-wide significant SNP within +/-50 kb\n")
cat("of the imputed index SNP.\n")

