#!/usr/bin/env Rscript

# ============================================================
# Revised Manhattan plots for Reviewer #9
#
# Uses dataset-specific Bonferroni thresholds based on
# LD-pruned marker counts (--indep-pairwise 1000kb 1 0.8).
#
# Output:
#   Main Figure 3: BB original + imputed
#   Supplementary Figure S1: DBH original + imputed
#   Supplementary Figure S2: FD original + imputed
#
# Only ONE horizontal line is drawn in each panel:
# the actual dataset-specific genome-wide threshold.
#
# Background points may be downsampled for plotting only;
# ALL genome-wide significant SNPs are always retained.
# Set MAX_BACKGROUND = Inf to plot every background SNP.
# ============================================================

library(data.table)
library(ggplot2)
library(cowplot)

BASE_DIR <- getwd()
OUTDIR <- file.path(BASE_DIR, "GWAS_Manhattan_revised_thresholds")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Final GWAS files and dataset-specific thresholds
# ------------------------------------------------------------------

cfg <- data.frame(
  Trait = c("BB", "BB", "DBH", "DBH", "FD", "FD"),
  Dataset = c("Original", "Imputed", "Original", "Imputed", "Original", "Imputed"),
  File = c(
    "GWAS-Bud/Bud-origin-GWAS.loco_fixed.mlma",
    "GWAS-Bud/Bud-imputed-PC0-GWASloco_chrall_fixed.mlma",
    "GWAS-DBH/DBH-origin-GWAS.loco_fixed.mlma",
    "GWAS-DBH/DBH-imputed-PC0-GWASloco_chrall_fixed.mlma",
    "GWAS-FD/FD-origin-GWAS.loco_fixed.mlma",
    "GWAS-FD/FD-imputed-PC0-GWASloco_chrall_fixed.mlma"
  ),
  Threshold = c(
    9.1044830475e-07,
    9.6525655554e-08,
    7.5243412439e-07,
    8.8580079049e-08,
    6.7133017361e-07,
    8.5449489952e-08
  ),
  stringsAsFactors = FALSE
)

# Plot all SNPs by default to preserve the original Manhattan-plot display.
# If runtime/file size becomes excessive, you may set this to e.g. 1000000;
# all genome-wide significant SNPs will still be retained.
MAX_BACKGROUND <- Inf

# Original manuscript-like chromosome palette
CHR_COLORS <- c(
  "#6BAED6", "#2171B5",
  "#FDAE6B", "#E6550D",
  "#66C2A5", "#238B45",
  "#9E9AC8", "#756BB1",
  "#FD8D3C", "#E6550D",
  "#41AB5D", "#006D2C"
)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

find_col <- function(nms, candidates) {
  low <- tolower(nms)
  for (cand in candidates) {
    idx <- which(low == tolower(cand))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

read_gwas <- function(path, threshold) {
  if (!file.exists(path)) stop("Missing GWAS file: ", path)

  hdr <- fread(path, nrows = 0)
  nms <- names(hdr)

  chr_col <- find_col(nms, c("Chr", "CHR", "chromosome"))
  bp_col  <- find_col(nms, c("bp", "BP", "pos", "POS"))
  p_col   <- find_col(nms, c("p", "P"))

  if (any(is.na(c(chr_col, bp_col, p_col)))) {
    stop(
      "Could not identify chromosome/BP/P columns in ", path,
      "\nColumns: ", paste(nms, collapse = ", ")
    )
  }

  d <- fread(
    path,
    select = c(chr_col, bp_col, p_col),
    showProgress = TRUE
  )
  setnames(d, c(chr_col, bp_col, p_col), c("CHR", "BP", "P"))

  d <- d[
    is.finite(P) & !is.na(P) &
    P > 0 & P <= 1 &
    !is.na(CHR) & !is.na(BP)
  ]

  d[, CHR := as.integer(gsub("[^0-9]", "", as.character(CHR)))]
  d <- d[CHR %between% c(1L, 12L)]
  d[, BP := as.numeric(BP)]
  d <- d[is.finite(BP)]

  d[, LOGP := -log10(P)]
  d[, Significant := P <= threshold]

  d
}

prepare_plot_data <- function(d, max_background = 1000000) {

  # chromosome lengths based on observed max BP in this GWAS
  chr_info <- d[, .(chr_len = max(BP, na.rm = TRUE)), by = CHR]
  setorder(chr_info, CHR)

  chr_info[, offset := c(0, head(cumsum(chr_len), -1))]
  chr_info[, center := offset + chr_len / 2]

  d <- merge(d, chr_info[, .(CHR, offset)], by = "CHR", all.x = TRUE)
  d[, BPcum := BP + offset]

  sig <- d[Significant == TRUE]
  bg  <- d[Significant == FALSE]

  if (is.finite(max_background) && nrow(bg) > max_background) {
    set.seed(20261009)
    bg <- bg[sample.int(nrow(bg), max_background)]
  }

  list(
    plot_data = rbindlist(list(bg, sig), use.names = TRUE),
    chr_info = chr_info
  )
}

make_panel <- function(path, threshold, title, point_size = 0.45) {

  d <- read_gwas(path, threshold)
  prep <- prepare_plot_data(d, MAX_BACKGROUND)
  plot_data <- prep$plot_data
  chr_info <- prep$chr_info

  threshold_y <- -log10(threshold)

  # color by chromosome; significant SNPs highlighted in red
  plot_data[, PointColor := ifelse(
    Significant,
    "#FF0000",
    CHR_COLORS[CHR]
  )]

  p <- ggplot(
    plot_data,
    aes(x = BPcum, y = LOGP)
  ) +
    geom_point(
      aes(color = PointColor),
      size = point_size,
      alpha = 0.85,
      show.legend = FALSE
    ) +
    scale_color_identity() +
    geom_hline(
      yintercept = threshold_y,
      linetype = "dashed",
      linewidth = 0.55,
      color = "black"
    ) +
    scale_x_continuous(
      breaks = chr_info$center,
      labels = chr_info$CHR,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    labs(
      title = title,
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic(base_family = "Times New Roman") +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 11
      ),
      axis.title.y = element_text(
        size = 11
      ),
      axis.text.x = element_text(
        size = 9
      ),
      axis.text.y = element_text(
        size = 9
      ),
      axis.line = element_line(
        linewidth = 0.6,
        color = "black"
      ),
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        color = "#E6E6E6"
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 8, 5, 8)
    )

  rm(d, prep, plot_data)
  gc()

  p
}

save_pair <- function(trait, file_prefix, width = 14, height = 9) {

  rows <- cfg[cfg$Trait == trait, , drop = FALSE]
  rows <- rows[match(c("Original", "Imputed"), rows$Dataset), ]

  p1 <- make_panel(
    rows$File[1],
    rows$Threshold[1],
    paste0(trait, "-original"),
    point_size = ifelse(trait == "BB", 0.55, 0.45)
  )

  p2 <- make_panel(
    rows$File[2],
    rows$Threshold[2],
    paste0(trait, "-imputed"),
    point_size = 0.38
  )

  combined <- plot_grid(
    p1, p2,
    ncol = 1,
    labels = c("A", "B"),
    label_size = 18,
    label_fontface = "bold",
    label_family = "Times New Roman",
    label_x = 0.005,
    label_y = 0.995,
    hjust = 0,
    vjust = 1,
    align = "v"
  )

  png_file <- file.path(
    OUTDIR,
    paste0(file_prefix, "_revised_threshold.png")
  )

  pdf_file <- file.path(
    OUTDIR,
    paste0(file_prefix, "_revised_threshold.pdf")
  )

  ggsave(
    png_file,
    combined,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    pdf_file,
    combined,
    width = width,
    height = height,
    bg = "white",
    device = cairo_pdf
  )

  cat("\nSaved:\n", png_file, "\n", pdf_file, "\n", sep = "")
}

# ------------------------------------------------------------------
# Generate revised figures
#
# Figure 3 = BB
# Supplementary Fig. S1 = DBH
# Supplementary Fig. S2 = FD
# ------------------------------------------------------------------

save_pair("BB",  "Figure3_BB",             width = 14, height = 9)
save_pair("DBH", "Supplementary_Fig_S1_DBH", width = 14, height = 9)
save_pair("FD",  "Supplementary_Fig_S2_FD",  width = 14, height = 9)

# ------------------------------------------------------------------
# Write the thresholds used in the plots
# ------------------------------------------------------------------

threshold_summary <- copy(as.data.table(cfg))
threshold_summary[, MinusLog10Threshold := -log10(Threshold)]

fwrite(
  threshold_summary[, .(
    Trait,
    Dataset,
    Threshold,
    MinusLog10Threshold
  )],
  file.path(OUTDIR, "Manhattan_thresholds_used.tsv"),
  sep = "\t"
)

cat("\nAll revised Manhattan plots completed successfully.\n")

