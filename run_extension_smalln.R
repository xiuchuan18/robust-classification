# run_extension_smalln.R ------------------------------------------------------
# Appendix G.2, the sample-size sweep: the STATISTICAL cost of stage-1
# flexibility.  Saves results/Sec6_smalln_raw.rds; report it with
# analyze_extension_smalln.R.
#
# Table G.4 shows that TPSC, skew-t and the BIC mixture are indistinguishable in
# accuracy at n = 800, so the choice among them rests on cost.  Table G.4's timing
# row measures the computational cost; this sweep measures the statistical one:
# how well each working density can be estimated when stage 1 sees few
# observations per class.  It produces the paragraph in Appendix G.2 that quotes
# the paired AUC difference 0.061 (0.002) at n = 80 falling to 5.6e-5 (7.7e-5)
# at n = 800.
#
# Two deliberate design differences from the main run (run_extension.R):
#   * the sweep runs on setting C3, which is skew-t's OWN correctly specified
#     mechanism.  That is what makes it decisive: if skew-t loses to TPSC on the
#     mechanism it is correctly specified for, the loss can only be estimation
#     variance, not misspecification.  (Set MECHS to add "C6: Contam", where no
#     family is correct, as a cross-check; only C3 is reported in Section 6.)
#   * the test set is a FRESH sample of size TEST_N drawn from the same
#     mechanism, rather than a 30% holdout.  At n = 80 a holdout would contain
#     24 points and the AUC noise would swamp the effect being measured.
#
# Seed contract as in R/13_ext_harness.R: one job = one (setting, n, replicate);
# the seed drives generation and fold assignment, and all arms then see
# identical training data and identical folds.
#
# Usage (from the project root):
#       Rscript run_extension_smalln.R
#       Sys.setenv(N_REPS = "300"); source("run_extension_smalln.R")

PROJ <- normalizePath(Sys.getenv("EXT_WD", getwd()), winslash = "/")
source(file.path(PROJ, "R", "00c_setup_ext.R")); ext_boot(PROJ)
use_experiment("exp1")

N_REPS <- as.integer(Sys.getenv("N_REPS", "300"))
NGRID  <- as.integer(strsplit(Sys.getenv("NGRID", "80,150,300,800"), ",")[[1]])
MECHS  <- strsplit(Sys.getenv("MECHS", "C3: SkewT"), "\\|")[[1]]
PP     <- as.integer(Sys.getenv("PP",   "100"))
PREL   <- as.integer(Sys.getenv("PREL", "20"))
TEST_N <- as.integer(Sys.getenv("TEST_N", "3000"))
NCORE  <- as.integer(Sys.getenv("NCORE", as.character(max(1, detectCores() - 2))))
OUT    <- Sys.getenv("OUT", "Sec6_smalln_raw.rds")
OUT    <- if (basename(OUT) == OUT) ext_res(OUT) else OUT
SEED_1B <- 771000000L

jobs <- expand.grid(mech = MECHS, n = NGRID, rep = seq_len(N_REPS),
                    stringsAsFactors = FALSE)
cat(sprintf("small-n sweep: %d settings x %d sample sizes x %d reps = %d jobs | arms: %s\n",
            length(MECHS), length(NGRID), N_REPS, nrow(jobs),
            paste(EXP1_ARMS, collapse = ", ")))

run_smalln_job <- function(mech, n, rep_id, mech_index) {
  seed <- SEED_1B + 1000000L * mech_index + 1000L * n + rep_id
  set.seed(seed)
  d  <- generate_cont_data(mech, n, PP, PREL)          # training sample
  dt <- generate_cont_data(mech, TEST_N, PP, PREL)     # fresh large test sample
  ytr <- as.numeric(as.character(d$y)); yte <- as.numeric(as.character(dt$y))
  # 3 folds when the training sample is small, so each fold keeps a usable size
  nf <- if (n < 150L) 3L else 5L
  foldid <- sample(rep(seq_len(nf), length.out = length(ytr)))

  do.call(rbind, lapply(EXP1_ARMS, function(arm) {
    tryCatch({
      t0 <- proc.time()[["elapsed"]]
      m  <- fit_srollr_any(d$X, ytr, policy = arm, foldid = foldid)
      pr <- predict_srollr_plugin(m, dt$X)
      sec <- proc.time()[["elapsed"]] - t0
      cm <- classification_metrics(pr, yte)
      sm <- selection_metrics(srollr_plugin_selected(m), d$signal, PP)
      data.frame(mechanism = mech, n = n, rep = rep_id, arm = arm,
                 auc = cm$auc, tpr = sm$tpr, fpr = sm$fpr, sel_f1 = sm$f1,
                 n_sel = sm$n_selected, seconds = sec,
                 ok = TRUE, seed = seed, stringsAsFactors = FALSE)
    }, error = function(e)
      data.frame(mechanism = mech, n = n, rep = rep_id, arm = arm,
                 auc = NA_real_, tpr = NA_real_, fpr = NA_real_,
                 sel_f1 = NA_real_, n_sel = NA_integer_,
                 seconds = NA_real_, ok = FALSE, seed = seed,
                 stringsAsFactors = FALSE))
  }))
}

cl <- makeCluster(NCORE); registerDoParallel(cl)
clusterExport(cl, c("PROJ", "PP", "PREL", "TEST_N", "jobs", "MECHS",
                    "SEED_1B", "run_smalln_job"), envir = environment())
cpp_exports <- c("d_tpsc_cpp", "tpsc_nll_cpp", "tpsc_transform_cpp",
                 "tpsc_rowloglik_cpp")

t0 <- proc.time()[["elapsed"]]
results <- foreach(k = seq_len(nrow(jobs)), .combine = "rbind",
                   .packages = c("Rcpp", "glmnet", "pROC", "sn", "MASS"),
                   .noexport = cpp_exports) %dopar% {
  if (!exists("EXT_READY", envir = globalenv())) {
    source(file.path(PROJ, "R", "00c_setup_ext.R")); ext_boot(PROJ)
    use_experiment("exp1")
    assign("EXT_READY", TRUE, envir = globalenv())
  }
  j <- jobs[k, ]
  run_smalln_job(j$mech, j$n, j$rep, match(j$mech, MECHS))
}
stopCluster(cl); registerDoSEQ()

cat(sprintf("elapsed %.1f min | rows %d / %d | failures %d\n",
            (proc.time()[["elapsed"]] - t0) / 60, nrow(results),
            nrow(jobs) * length(EXP1_ARMS), sum(!results$ok)))
attr(results, "config") <- list(n_grid = NGRID, mechs = MECHS, p = PP,
                                p_relevant = PREL, test_n = TEST_N,
                                n_reps = N_REPS)
saveRDS(results, OUT); cat("saved:", OUT, "\n")
cat("Report it with:  Rscript analyze_extension_smalln.R\n")
