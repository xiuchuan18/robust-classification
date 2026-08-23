# 08_exp_E1.R
# E1: how Psi is estimated -- initialization, optimizer, and the parameter box.
#
#   E1a  three estimators M0-M2 on common data: does M0 find the same optimum as a near-global search (M1), and can a different optimizer (M2) do better?
#   E1b  the shipped estimator under six different boxes: how sensitive are the fitted Likelihood and the fitted Distribution to each face of the box?
#
# Both blocks use common random numbers: for a given (config, n, rep) every method and every box sees the identical sample, so differences are attributable to the method/box alone.

E1_N    <- c(20, 50, 100, 400)
E1_REPS <- 500L

## E1a: estimator comparison
run_E1a <- function(configs = E_CONFIGS, n_grid = E1_N, reps = E1_REPS,
                    fitters = FITTERS, seed = 20260717L, quad_m = 801L,
                    verbose = TRUE) {
  out <- list(); k <- 0L
  for (cf in names(configs)) {
    P0 <- unname(configs[[cf]]); P1 <- class1_par(P0)
    for (n in n_grid) {
      for (r in seq_len(reps)) {
        set.seed(.key_seed(seed, cf, n, r))
        x0 <- rtpsc(n, P0[1], P0[2], P0[3], P0[4])
        x1 <- rtpsc(n, P1[1], P1[2], P1[3], P1[4])
        for (mn in names(fitters)) {
          f0 <- fitters[[mn]](x0); f1 <- fitters[[mn]](x1)
          k <- k + 1L
          out[[k]] <- data.frame(
            config = cf, n = n, rep = r, method = mn,
            nll0 = f0$nll, nll1 = f1$nll,
            theta_hat = f0$par[1], w_hat = f0$par[2],
            sigma_hat = f0$par[3], delta_hat = f0$par[4],
            tmse = tmse_tpsc(P0, P1, f0$par, f1$par, m = quad_m),
            kl0  = kl_tpsc(P0, f0$par, m = quad_m),
            elapsed = f0$elapsed + f1$elapsed,
            n_starts = f0$n_starts, n_err = f0$n_err + f1$n_err,
            fallback = f0$fallback | f1$fallback,
            stringsAsFactors = FALSE)
        }
      }
      if (verbose) message(sprintf("E1a  %-28s n = %4d  done", cf, n))
    }
  }
  raw <- do.call(rbind, out[seq_len(k)])
  key  <- paste(raw$config, raw$n, raw$rep, sep = "|")
  ref0 <- tapply(raw$nll0, key, min); ref1 <- tapply(raw$nll1, key, min)
  raw$nll_gap <- (raw$nll0 - ref0[key]) + (raw$nll1 - ref1[key])
  raw
}

summarise_E1a <- function(raw, tol = 1e-4) {
  agg <- function(d) with(d, data.frame(
    config = config[1], n = n[1], method = method[1], reps = nrow(d),
    pct_at_best   = 100 * mean(nll_gap <= tol),
    gap_median    = median(nll_gap),
    gap_q95       = unname(quantile(nll_gap, 0.95)),
    tmse_median   = median(tmse),
    kl0_median    = median(kl0),
    delta_median  = median(delta_hat),
    pct_fallback  = 100 * mean(fallback),
    pct_start_err = 100 * mean(n_err > 0),
    ms_per_fit    = 1000 * mean(elapsed),
    stringsAsFactors = FALSE))
  sm <- do.call(rbind, lapply(
    split(raw, list(raw$config, raw$n, raw$method), drop = TRUE), agg))
  sm[order(sm$config, sm$n, sm$method), ]
}

## E1b: box sensitivity
# The shipped box is  w in [0.05, 0.95], sigma > 1e-6, delta in [0.5, 100].
E1_BOXES <- list(
  "shipped:  w[.05,.95] d[.5,100] s>1e-6" = list(lo = c(-Inf, 0.05, 1e-6, 0.5),
                                                 up = c(Inf, 0.95, Inf, 100)),
  "d_max=30                             " = list(lo = c(-Inf, 0.05, 1e-6, 0.5),
                                                 up = c(Inf, 0.95, Inf, 30)),
  "d_max=500                            " = list(lo = c(-Inf, 0.05, 1e-6, 0.5),
                                                 up = c(Inf, 0.95, Inf, 500)),
  "d_min=0.2                            " = list(lo = c(-Inf, 0.05, 1e-6, 0.2),
                                                 up = c(Inf, 0.95, Inf, 100)),
  "w[.01,.99]                           " = list(lo = c(-Inf, 0.01, 1e-6, 0.5),
                                                 up = c(Inf, 0.99, Inf, 100)),
  "sigma > 1e-3                         " = list(lo = c(-Inf, 0.05, 1e-3, 0.5),
                                                 up = c(Inf, 0.95, Inf, 100))
)

run_E1b <- function(configs = E_CONFIGS, n_grid = E1_N, reps = E1_REPS,
                    boxes = E1_BOXES, seed = 20260717L, quad_m = 801L,
                    verbose = TRUE) {
  out <- list(); k <- 0L
  for (cf in names(configs)) {
    P0 <- unname(configs[[cf]])
    for (n in n_grid) {
      for (r in seq_len(reps)) {
        # identical seeds to run_E1a => identical class-0 samples => the blocks are directly comparable
        set.seed(.key_seed(seed, cf, n, r))
        x0 <- rtpsc(n, P0[1], P0[2], P0[3], P0[4])
        for (bn in names(boxes)) {
          bx <- boxes[[bn]]
          f0 <- fit_m0(x0, bx$lo, bx$up)
          bf <- boundary_flags(f0$par, bx$lo, bx$up)
          k <- k + 1L
          out[[k]] <- data.frame(
            config = cf, n = n, rep = r, box = bn,
            nll0 = f0$nll, delta_hat = f0$par[4], w_hat = f0$par[2],
            sigma_hat = f0$par[3],
            kl0 = kl_tpsc(P0, f0$par, m = quad_m),
            n_err = f0$n_err, fallback = f0$fallback, elapsed = f0$elapsed,
            t(bf), stringsAsFactors = FALSE)
        }
      }
      if (verbose) message(sprintf("E1b  %-28s n = %4d  done", cf, n))
    }
  }
  do.call(rbind, out[seq_len(k)])
}

summarise_E1b <- function(raw, ref_box = names(E1_BOXES)[1]) {
  key  <- paste(raw$config, raw$n, raw$rep, sep = "|")
  rb   <- raw[raw$box == ref_box, ]
  kr   <- paste(rb$config, rb$n, rb$rep, sep = "|")
  nref <- setNames(rb$nll0, kr); kref <- setNames(rb$kl0, kr)
  raw$d_nll <- raw$nll0 - nref[key]
  raw$d_kl0 <- raw$kl0  - kref[key]
  agg <- function(d) with(d, data.frame(
    config = config[1], n = n[1], box = box[1], reps = nrow(d),
    # likelihood level
    nll_median    = median(nll0),
    d_nll_median  = median(d_nll),
    d_nll_q95     = unname(quantile(abs(d_nll), 0.95)),
    # distribution level (against the truth)
    kl0_median    = median(kl0),
    d_kl0_median  = median(d_kl0),
    d_kl0_q95     = unname(quantile(abs(d_kl0), 0.95)),
    # mechanism: which face of the box is active, and how delta responds
    delta_median  = median(delta_hat),
    pct_delta_up  = 100 * mean(delta_up),
    pct_w_bound   = 100 * mean(w_lo | w_up),
    pct_sigma_lo  = 100 * mean(sigma_lo),
    # health
    pct_start_err = 100 * mean(n_err > 0),
    pct_fallback  = 100 * mean(fallback),
    ms_per_fit    = 1000 * mean(elapsed),
    stringsAsFactors = FALSE))
  sm <- do.call(rbind, lapply(
    split(raw, list(raw$config, raw$n, raw$box), drop = TRUE), agg))
  sm[order(sm$config, sm$n, sm$box), ]
}

run_E1 <- function(...) {
  a <- run_E1a(...); saveRDS(a, .res("E1a_raw.rds"))
  b <- run_E1b(...); saveRDS(b, .res("E1b_raw.rds"))
  sa <- summarise_E1a(a); sb <- summarise_E1b(b)
  write.csv(sa, .res("E1a_summary.csv"), row.names = FALSE)
  write.csv(sb, .res("E1b_summary.csv"), row.names = FALSE)
  list(E1a = sa, E1b = sb)
}
