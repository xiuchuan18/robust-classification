# analyze_extension.R ---------------------------------------------------------
# Report for Section 6 of the main article, whose tables live in Appendix G of
# the supplement: turns the three raw result objects written by run_extension.R
# into the manuscript's tables.
#
#   results/Sec6_continuous_raw.rds -> Table G.4 (AUC, selection F1, Maximum
#                                      Performance Gap, median time) and
#                                      Table G.3 (transformation sampling
#                                      variability of the three competitive
#                                      variants)
#   results/Sec6_discrete_raw.rds   -> Table G.5
#   results/Sec6_mixed_raw.rds      -> Table G.6, plus the directional test
#                                      behind the "gain rises monotonically with
#                                      rho" claim of Appendix G.3
#
# Every printed table uses the manuscript's setting labels (C1-C6, D1-D5,
# E1-E5) and the manuscript's variant names; the code ids stored in the raw
# files are mapped by paper_mech() / paper_arm() in R/14_ext_experiments.R.
#
# Reviewer 1's third comment asks for variable-selection accuracy and computing
# time alongside classification, so every selection metric that is recorded is
# also reported, and selection gets its own Performance-Gap treatment rather
# than only an AUC one.
#
# Usage (from the project root):
#       Rscript analyze_extension.R
#       source("analyze_extension.R")

PROJ <- normalizePath(Sys.getenv("EXT_WD", getwd()), winslash = "/")
source(file.path(PROJ, "R", "00c_setup_ext.R")); ext_boot(PROJ)
suppressMessages({library(dplyr); library(tidyr); library(tibble)})
filter <- dplyr::filter; select <- dplyr::select; summarise <- dplyr::summarise

# label rows/columns exactly as the manuscript does
.relabel <- function(m) {
  rownames(m) <- paper_mech(rownames(m))
  colnames(m) <- paper_arm(colnames(m))
  m
}

report <- function(exp, rds, table_no) {
  if (!file.exists(rds)) { cat("\n[missing]", rds, "\n"); return(invisible(NULL)) }
  use_experiment(exp)
  res <- readRDS(rds)
  cfg <- attr(res, "config")
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat(sprintf("MANUSCRIPT TABLE %s  [%s]  (n=%d p=%d |I|=%d reps=%d, rows=%d, failures=%d)\n",
              table_no, basename(rds), cfg$n, cfg$p, cfg$p_relevant, cfg$n_reps,
              nrow(res), sum(!res$ok)))
  res <- res %>% mutate(mechanism = factor(mechanism, levels = ALL_MECHANISMS),
                        arm = factor(arm, levels = ARM_ORDER))
  s <- res %>% filter(ok) %>% group_by(mechanism, arm) %>%
    summarise(auc = mean(auc), tpr = mean(tpr), fpr = mean(fpr),
              sel_f1 = mean(sel_f1), n_sel = mean(n_sel),
              sec = median(seconds), .groups = "drop")
  M <- function(col) {
    m <- s %>% select(mechanism, arm, all_of(col)) %>%
      pivot_wider(names_from = arm, values_from = all_of(col)) %>%
      column_to_rownames("mechanism"); .relabel(as.matrix(m))
  }
  A <- M("auc"); TPR <- M("tpr"); FPR <- M("fpr")
  SF <- M("sel_f1"); NSEL <- M("n_sel"); SEC <- M("sec")

  # Performance Gap (Appendix G.1): the shortfall of a variant relative to the
  # best variant available in that setting. Computed for BOTH headline criteria.
  R  <- sweep(-A,  1, apply(A,  1, max, na.rm = TRUE), "+")
  RS <- sweep(-SF, 1, apply(SF, 1, max, na.rm = TRUE), "+")
  stopifnot(all(R >= -1e-9, na.rm = TRUE), all(RS >= -1e-9, na.rm = TRUE))

  cat("\n-- CLASSIFICATION: AUC --\n");                 print(round(A, 4))
  cat("\n-- CLASSIFICATION: Performance Gap --\n");     print(round(R, 4))
  cat("\n-- SELECTION: F1 vs the true signal set --\n");print(round(SF, 3))
  cat("\n-- SELECTION: Performance Gap in F1 --\n");    print(round(RS, 3))
  cat("\n-- SELECTION: true positive rate --\n");       print(round(TPR, 3))
  cat("\n-- SELECTION: false positive rate --\n");      print(round(FPR, 3))
  cat("\n-- SELECTION: mean model size (true |I| = ", cfg$p_relevant, ") --\n",
      sep = "");                                        print(round(NSEL, 1))
  cat("\n-- COST: median seconds per fit-and-predict cycle --\n")
  print(round(SEC, 2))

  cat("\n-- summary across the settings (the last rows of the manuscript table) --\n")
  print(round(data.frame(
    worst_auc            = apply(A,  2, min),
    auc_max_perf_gap     = apply(R,  2, max),
    auc_mean_perf_gap    = apply(R,  2, mean),
    worst_selF1          = apply(SF, 2, min),
    selF1_max_perf_gap   = apply(RS, 2, max),
    selF1_mean_perf_gap  = apply(RS, 2, mean),
    worst_tpr            = apply(TPR, 2, min),
    largest_fpr          = apply(FPR, 2, max),
    med_sec              = apply(SEC, 2, median)), 4))

  invisible(list(A = A, R = R, SF = SF, RS = RS, TPR = TPR, FPR = FPR,
                 NSEL = NSEL, SEC = SEC, res = res))
}

e1 <- report("exp1", ext_res(Sys.getenv("RDS1", "Sec6_continuous_raw.rds")), "G.4")
e2 <- report("exp2", ext_res(Sys.getenv("RDS2", "Sec6_discrete_raw.rds")),   "G.5")
e3 <- report("exp3", ext_res(Sys.getenv("RDS3", "Sec6_mixed_raw.rds")),      "G.6")

## --- Table G.3: transformation sampling variability --------------------------
# sd across the Monte Carlo replicates of T(x) at fixed evaluation points,
# averaged over the grid.  The grid is fixed a priori (built once under
# GRID_SEED in R/13_ext_harness.R), so this is a genuine sampling sd of the
# induced transformation, not a bootstrap proxy.  The `identity` arm is exactly
# zero by construction, since T(x) = x carries no estimated parameter.
if (!is.null(e1)) {
  use_experiment("exp1")
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("MANUSCRIPT TABLE G.3  transformation sampling variability, S-bar(m)\n")
  gc_ <- grep("^g[0-9]+$", names(e1$res), value = TRUE)
  ts <- e1$res %>% filter(ok) %>% select(mechanism, arm, all_of(gc_)) %>%
    pivot_longer(all_of(gc_), names_to = "gpt", values_to = "tval") %>%
    filter(is.finite(tval)) %>%
    group_by(mechanism, arm, gpt) %>% summarise(sd = sd(tval), .groups = "drop") %>%
    group_by(mechanism, arm) %>% summarise(t_sd = mean(sd), .groups = "drop")
  TS <- ts %>% pivot_wider(names_from = arm, values_from = t_sd) %>%
    column_to_rownames("mechanism")
  TS <- TS[, intersect(EXP1_ARMS, colnames(TS)), drop = FALSE]
  TS <- .relabel(as.matrix(TS))
  print(round(TS, 3))
  cat("\nmax_m S-bar(m)  (the first column of Table G.3):\n")
  print(round(sort(apply(TS, 2, max)), 3))
  cat("\nave_m S-bar(m)  (the second column of Table G.3):\n")
  print(round(sort(apply(TS, 2, mean)), 3))
}

## --- Table G.6 follow-up: the directional test -------------------------------
if (!is.null(e3)) {
  use_experiment("exp3")
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("APPENDIX G.3 DIRECTIONAL TEST (mixed data)\n")
  cat("H: the advantage of the ADAPTIVE rule over each single-family variant\n")
  cat("   increases with rho, the fraction of the discrete block that is NOT a\n")
  cat("   one-parameter exponential family (i.e. unordered nominal), because for\n")
  cat("   Bernoulli/Poisson columns the LLR is affine in x and a continuous\n")
  cat("   plug-in therefore loses nothing.\n")
  cat("Paired per-replicate differences: every arm sees identical data and\n")
  cat("identical CV folds within a replicate.\n\n")
  w <- e3$res %>% filter(ok) %>%
    select(mechanism, rep, arm, auc, sel_f1, fpr, n_sel) %>%
    pivot_wider(names_from = arm, values_from = c(auc, sel_f1, fpr, n_sel))
  for (rival in c("tpsc", "binned", "identity")) {
    cat(sprintf("  adaptive - %-8s (AUC)     ", ARM_LABEL[[rival]]))
    o <- w %>% group_by(mechanism) %>%
      summarise(d = mean(auc_auto - .data[[paste0("auc_", rival)]]),
                se = sd(auc_auto - .data[[paste0("auc_", rival)]]) / sqrt(dplyr::n()),
                .groups = "drop")
    cat(paste(sprintf("%-14s %+.4f (%.4f)", paper_mech(o$mechanism), o$d, o$se),
              collapse = "   "), "\n")
  }
  for (rival in c("tpsc", "binned")) {
    cat(sprintf("  adaptive - %-8s (sel F1)  ", ARM_LABEL[[rival]]))
    o <- w %>% group_by(mechanism) %>%
      summarise(d = mean(sel_f1_auto - .data[[paste0("sel_f1_", rival)]]),
                se = sd(sel_f1_auto - .data[[paste0("sel_f1_", rival)]]) / sqrt(dplyr::n()),
                .groups = "drop")
    cat(paste(sprintf("%-14s %+.4f (%.4f)", paper_mech(o$mechanism), o$d, o$se),
              collapse = "   "), "\n")
  }
  cat("\n  Reading: a positive and INCREASING sequence across E1-E5 supports the\n")
  cat("  hypothesis; a flat sequence would mean the mixed-data experiment is not\n")
  cat("  earning its place as a separate table.\n")
}
