# visualize_results.R
# Boxplots of the twelve-Case simulation study in the exact style of Manuscript Figures 1 and 2, written as one multi-page PDF (one figure per page, no titles,
# ready to drop into the paper):
#   run_simulation_visualization.pdf: 7 pages,
#       1 AUC, 2 Accuracy, 3 F1-Score (sample classification), 4 Variable-selection TPR, 5 FPR, 6 selection F1, 7 Computation time

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

# Resolve the project root and anchor every path to it.  Nothing here calls
# setwd(), so sourcing this file leaves the caller's working directory alone.
.root <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())

## shared style
MODEL_ORDER <- c("RoLLR", "S-RoLLR", "MOKE", "GMM-NB",
                 "L1-RL", "SVM", "RF", "LightGBM")

.set1 <- c(RoLLR = "#E41A1C", `S-RoLLR` = "#377EB8", MOKE = "#4DAF4A",
           `GMM-NB` = "#984EA3", `L1-RL` = "#FF7F00", SVM = "#A65628",
           RF = "#F781BF", LightGBM = "#999999")
.blend_over_white <- function(hex, a = 0.8) {
  rgb <- col2rgb(hex)
  out <- round(a * rgb + (1 - a) * 255)
  apply(out, 2, function(c) sprintf("#%02X%02X%02X", c[1], c[2], c[3]))
}
MODEL_COLORS <- setNames(.blend_over_white(.set1, 0.8), names(.set1))

BASE_SIZE <- 14

manuscript_theme <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      legend.position = "none",
      strip.text      = element_text(face = "bold", colour = "black"),
      axis.text.x     = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y     = element_text(colour = "black"),
      axis.title.y    = element_text(face = "bold"))
}

# One boxplot page in the manuscript style (no title).
make_box <- function(df, yvar, facet_var, ylab, facet_ncol, log_y = FALSE) {
  df <- df[!is.na(df[[yvar]]), ]
  df$model <- factor(df$model, levels = MODEL_ORDER)
  df <- droplevels(df[!is.na(df$model), ])

  p <- ggplot(df, aes(x = model, y = .data[[yvar]], fill = model)) +
    geom_boxplot(colour = "black", linewidth = 0.4, median.linewidth = 0.6,
                 outlier.size = 0.55, outlier.colour = "grey25",
                 outlier.shape = 19) +
    scale_fill_manual(values = MODEL_COLORS, drop = TRUE) +
    facet_wrap(vars(.data[[facet_var]]), ncol = facet_ncol,
               scales = "free_y") +
    labs(x = NULL, y = ylab) +
    manuscript_theme()

  if (log_y)
    p <- p + scale_y_log10()
  p
}

## Simulation study: run_simulation_visualization.pdf  (7 pages, 10 x 8)
sim <- readRDS(file.path(.root, "results", "simulation_results.rds"))

# Facet order = Case NUMBER (a plain character column would sort Case10/11/12
# straight after Case1). ALL_SCENARIOS from the data generator is already in
# natural order, so it doubles as the factor-level order.
source(file.path(.root, "R", "03_data_gen.R"))
sim$scenario <- factor(sim$scenario, levels = ALL_SCENARIOS)
stopifnot(!any(is.na(sim$scenario)))     # every scenario name must be known

sel <- subset(sim, model %in% c("S-RoLLR", "L1-RL"))
sim_pages <- list(
  make_box(sim, "auc",     "scenario", "AUC",                 3),
  make_box(sim, "acc",     "scenario", "Accuracy",            3),
  make_box(sim, "f1",      "scenario", "F1-Score",            3),
  make_box(sel, "sel_tpr", "scenario", "True positive rate",  3),
  make_box(sel, "sel_fpr", "scenario", "False positive rate", 3),
  make_box(sel, "sel_f1",  "scenario", "Selection F1-Score",  3),
  make_box(sim, "seconds", "scenario", "Seconds (log scale)", 3, log_y = TRUE))

pdf(file.path(.root, "docs", "run_simulation_visualization.pdf"),
    width = 10, height = 8, onefile = TRUE)
for (p in sim_pages) print(p)
invisible(dev.off())
cat("Wrote docs/run_simulation_visualization.pdf (", length(sim_pages), "pages, 10x8 )\n")
