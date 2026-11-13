#!/usr/bin/env Rscript

# ============================================================
# GWAS QQ plots and genomic inflation factors
#
# Traits:
#   BB
#   DBH
#   FD
#
# Datasets:
#   Original
#   Imputed
#
# Figure layout:
#
#       Original          Imputed
# A     BB                B
# C     DBH               D
# E     FD                F
#
# X-axis:
#   identical across all six panels
#
# Y-axis:
#   identical within each trait pair
#
# Outputs:
#   GWAS_QQ_lambda/lambda_GC_summary.tsv
#   GWAS_QQ_lambda/GWAS_QQ_original_vs_imputed_fixed.png
#   GWAS_QQ_lambda/GWAS_QQ_original_vs_imputed_fixed.pdf
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(data.table)
library(ggplot2)
library(cowplot)
library(ggrastr)

# ============================================================
# 2. Working directory
#
# The SLURM script already changes to the GWAS directory.
# ============================================================

base_dir <- getwd()

cat("Working directory:\n")
cat(base_dir, "\n\n")


outdir <- file.path(
  base_dir,
  "GWAS_QQ_lambda"
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. Define final GWAS result files
#
# Keep this exact order:
#
# 1 BB Original
# 2 BB Imputed
# 3 DBH Original
# 4 DBH Imputed
# 5 FD Original
# 6 FD Imputed
# ============================================================

files <- data.frame(

  Trait = c(
    "BB",
    "BB",
    "DBH",
    "DBH",
    "FD",
    "FD"
  ),

  Dataset = c(
    "Original",
    "Imputed",
    "Original",
    "Imputed",
    "Original",
    "Imputed"
  ),

  File = c(

    "GWAS-Bud/Bud-origin-GWAS.loco_fixed.mlma",

    "GWAS-Bud/Bud-imputed-PC0-GWASloco_chrall_fixed.mlma",

    "GWAS-DBH/DBH-origin-GWAS.loco_fixed.mlma",

    "GWAS-DBH/DBH-imputed-PC0-GWASloco_chrall_fixed.mlma",

    "GWAS-FD/FD-origin-GWAS.loco_fixed.mlma",

    "GWAS-FD/FD-imputed-PC0-GWASloco_chrall_fixed.mlma"
  ),

  stringsAsFactors = FALSE
)


cat("GWAS datasets:\n")
print(files)
cat("\n")


# ============================================================
# 4. Function:
#    find a column name case-insensitively
# ============================================================

find_col <- function(names_vec, target) {

  idx <- which(
    tolower(names_vec) == tolower(target)
  )

  if (length(idx) == 0) {

    return(NA_character_)
  }

  return(
    names_vec[idx[1]]
  )
}


# ============================================================
# 5. Function:
#    generate QQ plotting data
#
# For large imputed datasets:
#
# - retain the 50,000 smallest P values exactly
# - sample up to 200,000 ranks across the full distribution
#
# Lambda GC is NOT calculated from these sampled points.
# Lambda GC is calculated from all valid GWAS SNPs.
# ============================================================

make_qq_data <- function(
    p,
    max_points = 200000,
    keep_tail = 50000) {


  # ----------------------------------------------------------
  # Keep valid P values
  # ----------------------------------------------------------

  p <- p[
    is.finite(p) &
      !is.na(p) &
      p >= 0 &
      p <= 1
  ]


  n <- length(p)


  if (n == 0) {

    stop(
      "No valid P values were found."
    )
  }


  # ----------------------------------------------------------
  # Replace exact P = 0 only for plotting
  #
  # Exact zero can occur because of numerical underflow.
  # ----------------------------------------------------------

  positive_p <- p[
    p > 0
  ]


  if (any(p == 0)) {

    if (length(positive_p) > 0) {

      replacement <- max(
        min(positive_p) / 10,
        1e-300
      )

    } else {

      replacement <- 1e-300
    }


    p[p == 0] <- replacement
  }


  # ----------------------------------------------------------
  # Sort observed P values from smallest to largest
  # ----------------------------------------------------------

  p <- sort(
    p,
    decreasing = FALSE
  )


  # ----------------------------------------------------------
  # Always retain strongest association signals
  # ----------------------------------------------------------

  tail_n <- min(
    keep_tail,
    n
  )


  idx_tail <- seq_len(
    tail_n
  )


  # ----------------------------------------------------------
  # Sample ranks across the complete P-value distribution
  # ----------------------------------------------------------

  sample_n <- min(
    max_points,
    n
  )


  idx_global <- unique(
    round(
      seq(
        from = 1,
        to = n,
        length.out = sample_n
      )
    )
  )


  # ----------------------------------------------------------
  # Combine indices
  # ----------------------------------------------------------

  idx <- sort(
    unique(
      c(
        idx_tail,
        idx_global
      )
    )
  )


  # ----------------------------------------------------------
  # Expected P-values
  # ----------------------------------------------------------

  expected_p <- idx / (n + 1)


  # ----------------------------------------------------------
  # Return QQ coordinates
  # ----------------------------------------------------------

  qq <- data.frame(

    Expected = -log10(
      expected_p
    ),

    Observed = -log10(
      p[idx]
    )
  )


  return(qq)
}


# ============================================================
# 6. Prepare containers
# ============================================================

n_datasets <- nrow(files)


lambda_results <- vector(
  mode = "list",
  length = n_datasets
)


qq_list <- vector(
  mode = "list",
  length = n_datasets
)


lambda_list <- numeric(
  n_datasets
)


# ============================================================
# 7. Main analysis loop
# ============================================================

for (i in seq_len(n_datasets)) {


  trait <- files$Trait[i]

  dataset <- files$Dataset[i]

  infile <- files$File[i]


  cat("\n")
  cat("====================================================\n")
  cat("Dataset ", i, " of ", n_datasets, "\n", sep = "")
  cat("Trait:   ", trait, "\n", sep = "")
  cat("Dataset: ", dataset, "\n", sep = "")
  cat("File:    ", infile, "\n", sep = "")
  cat("====================================================\n")


  # ==========================================================
  # Check file
  # ==========================================================

  if (!file.exists(infile)) {

    stop(
      paste(
        "Input file does not exist:",
        infile
      )
    )
  }


  # ==========================================================
  # Read file header
  # ==========================================================

  header <- fread(
    infile,
    nrows = 0
  )


  header_names <- names(
    header
  )


  cat(
    "Columns: ",
    paste(
      header_names,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )


  # ==========================================================
  # Identify P, beta and SE columns
  # ==========================================================

  p_col <- find_col(
    header_names,
    "p"
  )


  b_col <- find_col(
    header_names,
    "b"
  )


  se_col <- find_col(
    header_names,
    "se"
  )


  if (is.na(p_col)) {

    stop(
      paste(
        "Cannot find P-value column in:",
        infile
      )
    )
  }


  cat(
    "P column: ",
    p_col,
    "\n",
    sep = ""
  )


  if (
    !is.na(b_col) &&
      !is.na(se_col)
  ) {

    cat(
      "Beta column: ",
      b_col,
      "\n",
      sep = ""
    )

    cat(
      "SE column: ",
      se_col,
      "\n",
      sep = ""
    )

  } else {

    cat(
      "Beta/SE columns not found. ",
      "Lambda GC will be calculated from P values.\n",
      sep = ""
    )
  }


  # ==========================================================
  # Read only required columns
  # ==========================================================

  if (
    !is.na(b_col) &&
      !is.na(se_col)
  ) {

    select_cols <- c(
      p_col,
      b_col,
      se_col
    )

  } else {

    select_cols <- p_col
  }


  dat <- fread(
    infile,
    select = select_cols,
    showProgress = TRUE
  )


  # ==========================================================
  # Extract P values
  #
  # IMPORTANT:
  # Correct R syntax is dat[[column_name]]
  # ==========================================================

  p <- dat[[p_col]]


  valid_p <-
    is.finite(p) &
    !is.na(p) &
    p >= 0 &
    p <= 1


  n_gwas <- sum(
    valid_p
  )


  cat(
    "Valid GWAS P values: ",
    n_gwas,
    "\n",
    sep = ""
  )


  # ==========================================================
  # Calculate genomic inflation factor lambda GC
  #
  # Preferred:
  #
  # chi-square = (beta / SE)^2
  #
  # lambda GC =
  # median(observed chi-square) /
  # median(chi-square with df = 1)
  #
  # qchisq(0.5, df = 1) ~ 0.4549364
  # ==========================================================

  if (
    !is.na(b_col) &&
      !is.na(se_col)
  ) {


    beta <- dat[[b_col]]

    se <- dat[[se_col]]


    valid_stat <-
      is.finite(beta) &
      is.finite(se) &
      !is.na(beta) &
      !is.na(se) &
      se > 0


    chi2_obs <- (
      beta[valid_stat] /
        se[valid_stat]
    )^2


    lambda_gc <-
      median(
        chi2_obs,
        na.rm = TRUE
      ) /
      qchisq(
        0.5,
        df = 1
      )


    lambda_method <-
      "(beta/SE)^2"


    n_lambda <-
      length(
        chi2_obs
      )


  } else {


    # --------------------------------------------------------
    # Fallback:
    # derive chi-square statistics from P values
    # --------------------------------------------------------

    p_lambda <- p[
      valid_p &
        p > 0 &
        p <= 1
    ]


    chi2_obs <- qchisq(
      1 - p_lambda,
      df = 1
    )


    lambda_gc <-
      median(
        chi2_obs,
        na.rm = TRUE
      ) /
      qchisq(
        0.5,
        df = 1
      )


    lambda_method <-
      "P-value"


    n_lambda <-
      length(
        chi2_obs
      )
  }


  cat(
    "Lambda GC: ",
    sprintf(
      "%.4f",
      lambda_gc
    ),
    "\n",
    sep = ""
  )


  # ==========================================================
  # Generate QQ plotting data
  # ==========================================================

  qq <- make_qq_data(
    p[valid_p]
  )


  cat(
    "QQ plotting points retained: ",
    nrow(qq),
    "\n",
    sep = ""
  )


  cat(
    "Maximum expected -log10(P): ",
    sprintf(
      "%.3f",
      max(
        qq$Expected,
        na.rm = TRUE
      )
    ),
    "\n",
    sep = ""
  )


  cat(
    "Maximum observed -log10(P): ",
    sprintf(
      "%.3f",
      max(
        qq$Observed,
        na.rm = TRUE
      )
    ),
    "\n",
    sep = ""
  )


  # ==========================================================
  # Store QQ data and lambda
  # ==========================================================

  qq_list[[i]] <- qq

  lambda_list[i] <- lambda_gc


  # ==========================================================
  # Store summary
  # ==========================================================

  lambda_results[[i]] <- data.frame(

    Trait =
      trait,

    Dataset =
      dataset,

    GWAS_tested_SNPs =
      n_gwas,

    Lambda_GC =
      lambda_gc,

    Lambda_method =
      lambda_method,

    N_for_lambda =
      n_lambda,

    stringsAsFactors =
      FALSE
  )


  # ==========================================================
  # Memory cleanup
  # ==========================================================

  rm(
    dat,
    p,
    chi2_obs,
    qq,
    header
  )


  if (exists("beta")) {
    rm(beta)
  }


  if (exists("se")) {
    rm(se)
  }


  if (exists("p_lambda")) {
    rm(p_lambda)
  }


  gc()
}


# ============================================================
# 8. Verify all six datasets
# ============================================================

if (
  length(qq_list) != 6 ||
    any(
      vapply(
        qq_list,
        is.null,
        logical(1)
      )
    )
) {

  stop(
    "ERROR: Not all six QQ datasets were generated."
  )
}


if (
  length(lambda_list) != 6 ||
    any(
      !is.finite(
        lambda_list
      )
    )
) {

  stop(
    "ERROR: Not all six lambda GC values were generated."
  )
}


# ============================================================
# 9. Save lambda GC summary
# ============================================================

lambda_table <- do.call(
  rbind,
  lambda_results
)


lambda_output <- file.path(
  outdir,
  "lambda_GC_summary.tsv"
)


fwrite(
  lambda_table,
  lambda_output,
  sep = "\t"
)


cat("\n")
cat("====================================================\n")
cat("Lambda GC summary\n")
cat("====================================================\n")

print(
  lambda_table
)


cat(
  "\nSaved to: ",
  lambda_output,
  "\n",
  sep = ""
)


# ============================================================
# 10. Determine common X-axis
#
# For ~16 million SNPs:
#
# max expected -log10(P) is approximately 7.2
#
# Therefore this will normally give:
#
# global_xmax = 8
# ============================================================

global_xmax <- ceiling(
  max(
    vapply(
      qq_list,
      function(x) {

        max(
          x$Expected,
          na.rm = TRUE
        )
      },
      numeric(1)
    )
  )
)


# Ensure a reasonable minimum range

global_xmax <- max(
  global_xmax,
  4
)


cat("\n")
cat("====================================================\n")
cat("QQ plot coordinate ranges\n")
cat("====================================================\n")


cat(
  "Common X-axis: 0 - ",
  global_xmax,
  "\n",
  sep = ""
)


# ============================================================
# 11. Determine paired Y-axis ranges
#
# BB Original + BB Imputed:
# same Y-axis
#
# DBH Original + DBH Imputed:
# same Y-axis
#
# FD Original + FD Imputed:
# same Y-axis
# ============================================================

bb_ymax_raw <- ceiling(
  max(
    qq_list[[1]]$Observed,
    qq_list[[2]]$Observed,
    global_xmax,
    na.rm = TRUE
  )
)


dbh_ymax_raw <- ceiling(
  max(
    qq_list[[3]]$Observed,
    qq_list[[4]]$Observed,
    global_xmax,
    na.rm = TRUE
  )
)


fd_ymax_raw <- ceiling(
  max(
    qq_list[[5]]$Observed,
    qq_list[[6]]$Observed,
    global_xmax,
    na.rm = TRUE
  )
)


# Add small headroom

bb_ymax <- bb_ymax_raw * 1.03

dbh_ymax <- dbh_ymax_raw * 1.03

fd_ymax <- fd_ymax_raw * 1.03


cat(
  "BB Y-axis:  0 - ",
  sprintf("%.1f", bb_ymax),
  "\n",
  sep = ""
)


cat(
  "DBH Y-axis: 0 - ",
  sprintf("%.1f", dbh_ymax),
  "\n",
  sep = ""
)


cat(
  "FD Y-axis:  0 - ",
  sprintf("%.1f", fd_ymax),
  "\n",
  sep = ""
)


# ============================================================
# 12. Function:
#     draw one QQ plot
# ============================================================

make_qq_plot <- function(
    qq,
    trait,
    dataset,
    lambda_gc,
    x_max,
    y_max) {


  p <- ggplot(
    qq,
    aes(
      x = Expected,
      y = Observed
    )
  ) +


    # --------------------------------------------------------
    # Null expectation line
    # --------------------------------------------------------

    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.65
    ) +


    # --------------------------------------------------------
    # QQ points
    # --------------------------------------------------------
    ggrastr::rasterise(
      geom_point(
        size = 0.65,
        alpha = 0.55
      ),
      dpi = 600
    ) +


    # --------------------------------------------------------
    # Lambda GC
    # --------------------------------------------------------

    annotate(
      "text",

      x =
        x_max * 0.96,

      y =
        y_max * 0.035,

      label =
        paste0(
          "\u03BBGC = ",
          sprintf(
            "%.3f",
            lambda_gc
          )
        ),

      hjust = 1,

      vjust = 0,

      size = 4,

      family = "Arial"
    ) +


    # --------------------------------------------------------
    # Common X-axis
    # --------------------------------------------------------

    scale_x_continuous(

      limits = c(
        0,
        x_max
      ),

      breaks = seq(
        0,
        x_max,
        by = 2
      ),

      expand = expansion(
        mult = c(
          0.01,
          0.02
        )
      )
    ) +


    # --------------------------------------------------------
    # Trait-paired Y-axis
    # --------------------------------------------------------

    scale_y_continuous(

      limits = c(
        0,
        y_max
      ),

      expand = expansion(
        mult = c(
          0.01,
          0.02
        )
      )
    ) +


    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    labs(

      title =
        paste(
          trait,
          dataset
        ),

      x =
        expression(
          Expected~~-log[10](italic(P))
        ),

      y =
        expression(
          Observed~~-log[10](italic(P))
        )
    ) +


    # --------------------------------------------------------
    # Theme
    # --------------------------------------------------------

    theme_bw() +


    theme(

      text = element_text(
        family = "Arial",
        color = "black"
      ),


      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 5
        )
      ),


      axis.title = element_text(
        size = 11,
        face = "bold"
      ),


      axis.text = element_text(
        size = 9,
        color = "black"
      ),


      panel.border = element_rect(
        color = "black",
        linewidth = 0.8
      ),


      panel.grid.major = element_line(
        color = "#E5E5E5",
        linewidth = 0.35
      ),


      panel.grid.minor =
        element_blank(),


      plot.margin = margin(
        t = 8,
        r = 8,
        b = 8,
        l = 8
      )
    )


  return(p)
}


# ============================================================
# 13. Generate six explicit plots
# ============================================================


# A: BB Original

p_A <- make_qq_plot(

  qq =
    qq_list[[1]],

  trait =
    "BB",

  dataset =
    "Original",

  lambda_gc =
    lambda_list[1],

  x_max =
    global_xmax,

  y_max =
    bb_ymax
)


# B: BB Imputed

p_B <- make_qq_plot(

  qq =
    qq_list[[2]],

  trait =
    "BB",

  dataset =
    "Imputed",

  lambda_gc =
    lambda_list[2],

  x_max =
    global_xmax,

  y_max =
    bb_ymax
)


# C: DBH Original

p_C <- make_qq_plot(

  qq =
    qq_list[[3]],

  trait =
    "DBH",

  dataset =
    "Original",

  lambda_gc =
    lambda_list[3],

  x_max =
    global_xmax,

  y_max =
    dbh_ymax
)


# D: DBH Imputed

p_D <- make_qq_plot(

  qq =
    qq_list[[4]],

  trait =
    "DBH",

  dataset =
    "Imputed",

  lambda_gc =
    lambda_list[4],

  x_max =
    global_xmax,

  y_max =
    dbh_ymax
)


# E: FD Original

p_E <- make_qq_plot(

  qq =
    qq_list[[5]],

  trait =
    "FD",

  dataset =
    "Original",

  lambda_gc =
    lambda_list[5],

  x_max =
    global_xmax,

  y_max =
    fd_ymax
)


# F: FD Imputed

p_F <- make_qq_plot(

  qq =
    qq_list[[6]],

  trait =
    "FD",

  dataset =
    "Imputed",

  lambda_gc =
    lambda_list[6],

  x_max =
    global_xmax,

  y_max =
    fd_ymax
)


# ============================================================
# 14. Combine panels
#
# A BB Original       B BB Imputed
# C DBH Original      D DBH Imputed
# E FD Original       F FD Imputed
# ============================================================

combined_qq <- plot_grid(

  p_A,
  p_B,

  p_C,
  p_D,

  p_E,
  p_F,

  ncol = 2,

  nrow = 3,

  labels = c(
    "A",
    "B",
    "C",
    "D",
    "E",
    "F"
  ),

  label_size = 17,

  label_fontface = "bold",

  label_family = "Arial",

  label_x = 0.01,

  label_y = 0.995,

  hjust = 0,

  vjust = 1,

  align = "hv",

  axis = "tblr",

  rel_widths = c(
    1,
    1
  ),

  rel_heights = c(
    1,
    1,
    1
  )
)


# ============================================================
# 15. Define output filenames
# ============================================================

png_file <- file.path(
  outdir,
  "GWAS_QQ_original_vs_imputed_fixed.png"
)


pdf_file <- file.path(
  outdir,
  "GWAS_QQ_original_vs_imputed_fixed.pdf"
)


# ============================================================
# 16. Save PNG
# ============================================================

ggsave(

  filename =
    png_file,

  plot =
    combined_qq,

  width =
    10.5,

  height =
    13.5,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


# ============================================================
# 17. Save PDF
# ============================================================

ggsave(
  filename = pdf_file,
  plot = combined_qq,
  width = 10.5,
  height = 13.5,
  units = "in",
  bg = "white",
  device = cairo_pdf
)

# ============================================================
# 18. Final report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("GWAS QQ analysis completed successfully\n")
cat("====================================================\n")


cat(
  "Common X-axis: 0 - ",
  global_xmax,
  "\n",
  sep = ""
)


cat(
  "BB Y-axis:     0 - ",
  sprintf(
    "%.1f",
    bb_ymax
  ),
  "\n",
  sep = ""
)


cat(
  "DBH Y-axis:    0 - ",
  sprintf(
    "%.1f",
    dbh_ymax
  ),
  "\n",
  sep = ""
)


cat(
  "FD Y-axis:     0 - ",
  sprintf(
    "%.1f",
    fd_ymax
  ),
  "\n",
  sep = ""
)


cat("\nOutput files:\n")


cat(
  "1. ",
  lambda_output,
  "\n",
  sep = ""
)


cat(
  "2. ",
  png_file,
  "\n",
  sep = ""
)


cat(
  "3. ",
  pdf_file,
  "\n",
  sep = ""
)


cat("\nLambda GC values:\n")


for (
  i in seq_len(
    nrow(
      lambda_table
    )
  )
) {

  cat(

    lambda_table$Trait[i],
    " ",
    lambda_table$Dataset[i],
    ": lambda GC = ",
    sprintf(
      "%.3f",
      lambda_table$Lambda_GC[i]
    ),
    "\n",
    sep = ""
  )
}


cat("\n")
cat("====================================================\n")
cat("Finished.\n")
cat("====================================================\n")
