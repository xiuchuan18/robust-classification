# analyze_extension_smalln.R --------------------------------------------------
# Report for the Appendix G.2 sample-size sweep produced by
# run_extension_smalln.R (results/Sec6_smalln_raw.rds).
#
# Reports classification AND selection against n, and the paired per-replicate
# differences that the shared-data design permits.  The paired column
# `vs_skewt` on setting C3 is the quantity quoted in Appendix G.2: 0.061 (0.002)
# at n = 80, decreasing monotonically to 5.6e-5 (7.7e-5) at n = 800.
#
# Kept separate from the runner so the analysis is never coupled to a six-hour
# process.  Tolerates result files that lack tpr / n_sel (produced before those
# columns were added to the harness).
#
# Usage (from the project root):
#       Rscript analyze_extension_smalln.R
#       source("analyze_extension_smalln.R")

PROJ <- normalizePath(Sys.getenv("EXT_WD", getwd()), winslash = "/")
source(file.path(PROJ, "R", "00c_setup_ext.R")); ext_boot(PROJ)
use_experiment("exp1")
suppressMessages({library(dplyr); library(tidyr); library(tibble)})
filter <- dplyr::filter; select <- dplyr::select; summarise <- dplyr::summarise

rds <- ext_res(Sys.getenv("RDS", "Sec6_smalln_raw.rds"))
if (!file.exists(rds))
  stop("missing ", rds, "\n  Produce it with:  Rscript run_extension_smalln.R",
       call. = FALSE)

res <- readRDS(rds)
cfg <- attr(res, "config")
cat(sprintf("Appendix G.2 sample-size sweep [%s]: p=%d |I|=%d test_n=%d reps=%d rows=%d failures=%d\n",
            basename(rds), cfg$p, cfg$p_relevant, cfg$test_n, cfg$n_reps,
            nrow(res), sum(!res$ok)))
have <- intersect(c("auc", "sel_f1", "tpr", "fpr"), names(res))
cat("metrics available:", paste(have, collapse = ", "), "\n")

res   <- res %>% filter(ok) %>% mutate(arm = factor(arm, levels = EXP1_ARMS))
mechs <- unique(res$mechanism)

s <- res %>% group_by(mechanism, n, arm) %>%
  summarise(across(all_of(have), mean), .groups = "drop")

for (m in mechs) for (v in have) {
  cat(sprintf("\n=== setting %s : %s (rows = training sample size n) ===\n",
              paper_mech(m), v))
  tab <- s %>% filter(mechanism == m) %>% select(n, arm, all_of(v)) %>%
    pivot_wider(names_from = arm, values_from = all_of(v)) %>%
    column_to_rownames("n")
  colnames(tab) <- paper_arm(colnames(tab))
  print(round(tab, 4))
}

## --- paired differences ------------------------------------------------------
# Every arm sees identical training data and identical CV folds within a
# replicate, so the per-replicate difference removes the data-to-data variation
# and makes even a 0.001 AUC gap estimable.
cat("\n", strrep("=", 74), "\n", sep = "")
cat("PAIRED DIFFERENCES, TPSC minus rival (positive = TPSC better); se in the\n")
cat("adjacent column.\n")
for (v in intersect(c("auc", "sel_f1"), have)) {
  w <- res %>% select(mechanism, n, rep, arm, all_of(v)) %>%
    pivot_wider(names_from = arm, values_from = all_of(v))
  cat(sprintf("\n-- %s --\n", v))
  out <- w %>% group_by(mechanism, n) %>%
    summarise(vs_skewt = mean(tpsc - skewt),
              se_skewt = sd(tpsc - skewt) / sqrt(dplyr::n()),
              vs_mixture = mean(tpsc - mixture),
              se_mixture = sd(tpsc - mixture) / sqrt(dplyr::n()),
              vs_studentt = mean(tpsc - studentt),
              se_studentt = sd(tpsc - studentt) / sqrt(dplyr::n()),
              .groups = "drop") %>%
    mutate(mechanism = paper_mech(mechanism))
  print(as.data.frame(out), digits = 3, row.names = FALSE)
}
cat("\nReading: on setting C3 the skew-t model is CORRECTLY SPECIFIED, so any gap\n")
cat("there is pure estimation variance, not misspecification. A gap that shrinks\n")
cat("monotonically to zero as n grows is the signature of that variance.\n")
