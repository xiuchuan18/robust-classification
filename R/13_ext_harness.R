# 13_ext_harness.R
# Harness for the working-density study of Section 6. Defines the deterministic seed scheme, the fixed T-evaluation grid, and the function that executes one
# replicate of one mechanism across all arms.
#
# * One job = one (mechanism, replicate). Its seed is
#       SEED_BASE + 1e6 * mechanism_index + replicate
#   so seeds are unique, far apart, and independent of how foreach schedules the jobs across workers. Re-running any single job in isolation reproduces it.
# * The seed is consumed by the data generation and the train/test split only. Every arm then sees identical training data, an identical split, and an
#   identical vector of cross-validation fold labels (`foldid`, derived from the same stream BEFORE the arm loop starts). The comparison is therefore fully
#   paired, which is what makes small regrets estimable at moderate replication.
# * No stage-1 fitter consumes randomness (all EM/optimizer starts are deterministic functions of the data), and cv.glmnet is given `foldid`, so
#   the arm loop is deterministic. Consequently the arm order cannot affect any result.
# * The T-evaluation grids are precomputed ONCE in the master under grid_seed and exported to the workers, so no worker ever calls set.seed() itself.

SEED_BASE <- 20260728L
GRID_SEED <- 424242L

ARM_LABEL <- c(identity = "L1-RL", gaussian = "Gaussian", studentt = "Student-t",
               snorm = "skew-normal", skewt = "skew-t", tpsc = "TPSC",
               mixture = "mixture (BIC)", poisson = "Poisson",
               nbinom = "negative binomial", empirical = "empirical pmf",
               binned = "binned", auto = "adaptive")

## T-evaluation grid for the transform-stability diagnostic
# A fixed grid per mechanism, computed once from a large reference sample so that the sd of T(x) taken across replicates is a genuine sampling standard
# deviation of the transform at fixed evaluation points.
TGRID_NPT <- 25L

build_tgrids <- function(p, p_relevant, npt = TGRID_NPT, n_ref = 4000L) {
  set.seed(GRID_SEED)
  out <- lapply(ALL_MECHANISMS, function(m) {
    v <- generate_any(m, n_ref, p, p_relevant)$X[, 1]
    u <- unique(v)
    if (all(abs(v - round(v)) < 1e-8) && length(u) <= npt) sort(u)
    else as.numeric(seq(quantile(v, 0.01), quantile(v, 0.99), length.out = npt))
  })
  names(out) <- ALL_MECHANISMS
  out
}

## one job = one (mechanism, replicate), all arms
run_plugin_job <- function(mech, rep_id, n, p, p_relevant, tgrid,
                           mech_index, train_frac = 0.7, nfolds = 5L) {
  seed <- SEED_BASE + 1000000L * as.integer(mech_index) + as.integer(rep_id)
  set.seed(seed)

  d   <- generate_any(mech, n, p, p_relevant)
  idx <- sample(n, floor(train_frac * n))
  ytr <- as.numeric(as.character(d$y[idx]))
  yte <- as.numeric(as.character(d$y[-idx]))
  Xtr <- d$X[idx, , drop = FALSE]; Xte <- d$X[-idx, , drop = FALSE]
  # shared CV folds, so every arm is cross-validated identically
  foldid <- sample(rep(seq_len(nfolds), length.out = length(ytr)))

  blank <- function(arm, msg = NA_character_) data.frame(
    mechanism = mech, block = mech_block(mech), rep = rep_id, arm = arm,
    auc = NA_real_, acc = NA_real_, f1 = NA_real_, seconds = NA_real_,
    n_sel = NA_integer_, tpr = NA_real_, fpr = NA_real_, sel_f1 = NA_real_,
    t_absmax = NA_real_, ok = FALSE, error = msg, stringsAsFactors = FALSE)

  rows <- lapply(arm_coverage(mech), function(arm) {
    tryCatch({
      t0 <- proc.time()[["elapsed"]]
      m  <- fit_srollr_any(Xtr, ytr, policy = arm, foldid = foldid)
      pr <- predict_srollr_plugin(m, Xte)
      secs <- proc.time()[["elapsed"]] - t0
      cm <- classification_metrics(pr, yte)
      sm <- selection_metrics(srollr_plugin_selected(m), d$signal, p)
      Tte <- apply_plugin_transform(Xte, m$transform)
      r <- data.frame(
        mechanism = mech, block = mech_block(mech), rep = rep_id, arm = arm,
        auc = cm$auc, acc = cm$acc, f1 = cm$f1, seconds = secs,
        n_sel = sm$n_selected, tpr = sm$tpr, fpr = sm$fpr, sel_f1 = sm$f1,
        t_absmax = max(abs(Tte)), ok = TRUE, error = NA_character_,
        stringsAsFactors = FALSE)
      # T(x) on the fixed grid for column 1, used for the transform-stability diagnostic (sd across replicates). Refit on the TRAIN split of column 1.
      tg <- m$transform$fns[[1]](tgrid)
      attr(r, "tgrid_values") <- as.numeric(tg)
      r
    }, error = function(e) blank(arm, conditionMessage(e)))
  })

  pad <- function(v) { out <- rep(NA_real_, TGRID_NPT)
    if (!is.null(v)) out[seq_len(min(length(v), TGRID_NPT))] <-
      v[seq_len(min(length(v), TGRID_NPT))]
    out }
  tg_mat <- do.call(rbind, lapply(rows, function(r) pad(attr(r, "tgrid_values"))))
  colnames(tg_mat) <- paste0("g", seq_len(TGRID_NPT))
  res <- do.call(rbind, lapply(rows, function(r) { attr(r, "tgrid_values") <- NULL; r }))
  res$seed <- seed
  cbind(res, tg_mat)
}
