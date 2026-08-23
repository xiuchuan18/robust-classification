# audit_extension_timing.R ----------------------------------------------------
# Reproducibility check behind the "Median time (s)" row of manuscript Table G.4.
#
# The timings were collected on a shared desktop that can suspend mid-run, so
# the question is whether any wall-clock interruption leaked into the reported
# cost. Two different quantities are at stake and they behave differently:
#
#  (a) the per-arm `seconds` column is a proc.time() elapsed difference measured
#      AROUND one fit-and-predict cycle. proc.time() elapsed follows the system
#      clock, so a suspend that happens mid-fit IS charged to whichever arm was
#      in flight. With NCORE workers, at most NCORE cells can be inflated per
#      suspend.
#
#  (b) the run-level "elapsed" printed by run_extension.R is a wall-clock total
#      for the whole loop, so it would necessarily INCLUDE any suspend.
#
# The script compares the run-level total against the ideal packing of the sum
# of per-arm times over NCORE workers, and lists the largest per-arm cells. On
# the archived runs the two agree to within 2.6, 1.1 and 7.5 minutes and no
# single cell exceeds 110 s, so neither quantity is contaminated; the manuscript
# nevertheless quotes per-arm MEDIANS over 300 replicates, which a handful of
# inflated cells could not move in any case.
#
# Usage (from the project root):
#       Rscript audit_extension_timing.R

PROJ <- normalizePath(Sys.getenv("EXT_WD", getwd()), winslash = "/")
suppressMessages(library(dplyr))
source(file.path(PROJ, "R", "13_ext_harness.R"))        # ARM_LABEL
source(file.path(PROJ, "R", "14_ext_experiments.R"))   # paper_mech() / paper_arm()
.res <- function(f) file.path(PROJ, "results", f)

NCORE <- 14                             # cores used for the archived runs

# The small-n sweep is excluded here: its wall-clock log covers the original
# run, whose job list was a superset of the C3 sweep reported in Section 6, so
# `reported` and the retained rows are not comparable for that experiment.
files <- c(continuous = "Sec6_continuous_raw.rds",
           discrete   = "Sec6_discrete_raw.rds",
           mixed      = "Sec6_mixed_raw.rds")
reported <- c(continuous = 342.9, discrete = 31.5, mixed = 35.8)   # minutes

cat(sprintf("%-11s %10s %10s %10s %8s %8s %9s %9s\n", "experiment", "sum_core_h",
            "ideal_min", "report_min", "gap_min", "max_s", "n>600s", "n>3600s"))
for (nm in names(files)) {
  fp <- .res(files[nm])
  if (!file.exists(fp)) { cat(sprintf("%-11s [missing %s]\n", nm, files[nm])); next }
  r <- readRDS(fp)
  s <- r$seconds[r$ok & is.finite(r$seconds)]
  sum_h <- sum(s) / 3600
  ideal <- sum(s) / NCORE / 60          # perfect packing, no overhead, no sleep
  cat(sprintf("%-11s %10.2f %10.1f %10.1f %8.1f %8.0f %9d %9d\n",
              nm, sum_h, ideal, reported[nm], reported[nm] - ideal, max(s),
              sum(s > 600), sum(s > 3600)))
}

cat("\n-- the largest per-arm timings, by experiment --\n")
for (nm in names(files)) {
  fp <- .res(files[nm]); if (!file.exists(fp)) next
  r <- readRDS(fp)
  d <- r[r$ok & is.finite(r$seconds), ]
  d <- d[order(-d$seconds), ][1:5, ]
  cat(sprintf("\n%s:\n", nm))
  print(data.frame(setting = paper_mech(d$mechanism), arm = paper_arm(d$arm),
                   rep = d$rep,
                   seconds = round(d$seconds, 1)), row.names = FALSE)
}

cat("\n-- median vs mean vs max per arm (continuous run): the gap between the\n")
cat("   median and the max is what makes the Table G.4 median robust --\n")
fp <- .res(files[["continuous"]])
if (file.exists(fp)) {
  r <- readRDS(fp)
  print(r %>% filter(ok) %>% mutate(arm = paper_arm(arm)) %>% group_by(arm) %>%
          summarise(median_s = round(median(seconds), 2),
                    mean_s   = round(mean(seconds), 2),
                    p99_s    = round(quantile(seconds, 0.99), 1),
                    max_s    = round(max(seconds), 1), .groups = "drop") %>%
          as.data.frame(), row.names = FALSE)
}
