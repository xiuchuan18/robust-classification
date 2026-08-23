# 09_exp_E2.R
# E2: finite-sample behavior of the TPSC maximum-likelihood estimator.

E2_N    <- c(10, 20, 30, 50, 100, 200, 500, 1000)
E2_REPS <- 1000L

## the benchmark the paper quotes
ORACLE_BENCHMARK <- qchisq(0.5, 4) / qchisq(0.5, 3)

## E2a: parameter recovery
run_E2a <- function(configs = E_CONFIGS, n_grid = E2_N, reps = E2_REPS,
                    seed = 20260715L, quad_m = 801L, verbose = TRUE) {
  set.seed(seed)
  out <- vector("list", length(configs) * length(n_grid) * reps); k <- 0L
  for (cf in names(configs)) {
    P0 <- unname(configs[[cf]])
    for (n in n_grid) {
      for (r in seq_len(reps)) {
        x  <- rtpsc(n, P0[1], P0[2], P0[3], P0[4])
        ft <- fit_m0(x)
        bf <- boundary_flags(ft$par)
        k  <- k + 1L
        out[[k]] <- data.frame(
          config = cf, n = n, rep = r,
          theta_hat = ft$par[1], w_hat = ft$par[2],
          sigma_hat = ft$par[3], delta_hat = ft$par[4],
          theta0 = P0[1], w0 = P0[2], sigma0 = P0[3], delta0 = P0[4],
          nll = ft$nll, fallback = ft$fallback, n_err = ft$n_err,
          elapsed = ft$elapsed,
          kl = kl_tpsc(P0, ft$par, m = quad_m),
          t(bf), stringsAsFactors = FALSE)
      }
      if (verbose) message(sprintf("E2a  %-28s n = %5d  done", cf, n))
    }
  }
  do.call(rbind, out[seq_len(k)])
}

summarise_E2a <- function(raw) {
  agg <- function(d) {
    n1 <- d$n[1]
    with(d, data.frame(
      config = config[1], n = n1, reps = nrow(d),
      theta_bias = mean(theta_hat - theta0),
      theta_rmse = sqrt(mean((theta_hat - theta0)^2)),
      w_bias     = mean(w_hat - w0),
      w_rmse     = sqrt(mean((w_hat - w0)^2)),
      sigma_bias = mean(sigma_hat - sigma0),
      sigma_rmse = sqrt(mean((sigma_hat - sigma0)^2)),
      delta_median = median(delta_hat),
      delta_iqr    = IQR(delta_hat),
      kl_median    = median(kl),
      pct_delta_up = 100 * mean(delta_up),
      pct_delta_lo = 100 * mean(delta_lo),
      pct_w_bound  = 100 * mean(w_lo | w_up),
      pct_fallback = 100 * mean(fallback),
      pct_start_err = 100 * mean(n_err > 0),
      ms_per_fit   = 1000 * mean(elapsed),
      stringsAsFactors = FALSE))
  }
  sm <- do.call(rbind, lapply(split(raw, list(raw$config, raw$n), drop = TRUE), agg))
  sm[order(sm$config, sm$n), ]
}


## E2a paper tables: full-grid convergence rates
#
# Both blocks are reported together as supplement Table F.4:
#   parameter block:  log(RMSE_j) = a_j + b_j log(n), j in {theta,w,sigma}; theory b=-1/2.
#   divergence block: log(median KL) = a + b log(n);                        theory b=-1.
#
# delta is deliberately excluded from the rate table because it is weakly identified in the light-tailed regime and is therefore not uniformly covered
# by the regular root-n approximation.  Its median and IQR remain archived in E2a_summary.csv.
.fit_loglog_rate <- function(n, y, metric) {
  ok <- is.finite(n) & is.finite(y) & n > 0 & y > 0
  if (sum(ok) < 3L) {
    stop("Need at least three positive finite grid points for ", metric,
         "; found ", sum(ok), ".", call. = FALSE)
  }
  fit <- lm(log(y[ok]) ~ log(n[ok]))
  data.frame(
    metric    = metric,
    slope     = unname(coef(fit)[2]),
    intercept = unname(coef(fit)[1]),
    r2        = summary(fit)$r.squared,
    n_grid    = sum(ok),
    n_min     = min(n[ok]),
    n_max     = max(n[ok]),
    stringsAsFactors = FALSE)
}

E2a_parameter_rate <- function(sm) {
  required <- c("config", "n", "theta_rmse", "w_rmse", "sigma_rmse")
  missing  <- setdiff(required, names(sm))
  if (length(missing))
    stop("E2a summary is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  out <- lapply(split(sm, sm$config, drop = TRUE), function(d) {
    d <- d[order(d$n), , drop = FALSE]
    th <- .fit_loglog_rate(d$n, d$theta_rmse, "theta_rmse")
    ww <- .fit_loglog_rate(d$n, d$w_rmse,     "w_rmse")
    ss <- .fit_loglog_rate(d$n, d$sigma_rmse, "sigma_rmse")
    data.frame(
      config = d$config[1],
      theta_slope = th$slope, theta_r2 = th$r2,
      w_slope     = ww$slope, w_r2     = ww$r2,
      sigma_slope = ss$slope, sigma_r2 = ss$r2,
      theory_slope = -0.5,
      n_grid = th$n_grid,
      n_min  = th$n_min,
      n_max  = th$n_max,
      stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans[order(ans$config), ]
}

E2a_kl_rate <- function(sm) {
  required <- c("config", "n", "kl_median")
  missing  <- setdiff(required, names(sm))
  if (length(missing))
    stop("E2a summary is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  out <- lapply(split(sm, sm$config, drop = TRUE), function(d) {
    d <- d[order(d$n), , drop = FALSE]
    rr <- .fit_loglog_rate(d$n, d$kl_median, "kl_median")
    data.frame(
      config = d$config[1],
      kl_slope = rr$slope,
      kl_r2 = rr$r2,
      theory_slope = -1,
      n_grid = rr$n_grid,
      n_min  = rr$n_min,
      n_max  = rr$n_max,
      stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans[order(ans$config), ]
}

## E2b: transform fidelity + delta-oracle
# Two independent samples of size n from Psi0 (class 0) and Psi1 (class 1, mode shifted by MODE_SHIFT * sigma).  Both classes are fitted twice: by M0 (delta
# estimated) and with delta fixed at the truth.  If the two agree, the weak identification of delta demonstrably does not propagate to the classifier's inputs.
run_E2b <- function(configs = E_CONFIGS, n_grid = E2_N, reps = E2_REPS,
                    seed = 20260716L, quad_m = 801L, verbose = TRUE) {
  set.seed(seed)
  out <- vector("list", length(configs) * length(n_grid) * reps); k <- 0L
  for (cf in names(configs)) {
    P0 <- unname(configs[[cf]]); P1 <- class1_par(P0)
    for (n in n_grid) {
      for (r in seq_len(reps)) {
        x0 <- rtpsc(n, P0[1], P0[2], P0[3], P0[4])
        x1 <- rtpsc(n, P1[1], P1[2], P1[3], P1[4])
        f0 <- fit_m0(x0);                  f1 <- fit_m0(x1)
        o0 <- fit_oracle_delta(x0, P0[4]); o1 <- fit_oracle_delta(x1, P1[4])
        k  <- k + 1L
        out[[k]] <- data.frame(
          config = cf, n = n, rep = r,
          tmse_mle    = tmse_tpsc(P0, P1, f0$par, f1$par, m = quad_m),
          tmse_oracle = tmse_tpsc(P0, P1, o0$par, o1$par, m = quad_m),
          delta_hat0  = f0$par[4], delta0 = P0[4],
          kl0_mle     = kl_tpsc(P0, f0$par, m = quad_m),
          kl0_oracle  = kl_tpsc(P0, o0$par, m = quad_m),
          fallback    = f0$fallback | f1$fallback,
          stringsAsFactors = FALSE)
      }
      if (verbose) message(sprintf("E2b  %-28s n = %5d  done", cf, n))
    }
  }
  do.call(rbind, out[seq_len(k)])
}

summarise_E2b <- function(raw) {
  agg <- function(d) with(d, data.frame(
    config = config[1], n = n[1], reps = nrow(d),
    tmse_mle_median    = median(tmse_mle),
    tmse_oracle_median = median(tmse_oracle),
    tmse_ratio         = median(tmse_mle) / median(tmse_oracle),
    kl0_mle_median     = median(kl0_mle),
    kl0_oracle_median  = median(kl0_oracle),
    kl_oracle_ratio    = median(kl0_mle) / median(kl0_oracle),
    oracle_benchmark   = ORACLE_BENCHMARK,
    delta_hat_iqr      = IQR(delta_hat0),
    stringsAsFactors = FALSE))
  sm <- do.call(rbind, lapply(split(raw, list(raw$config, raw$n), drop = TRUE), agg))
  sm[order(sm$config, sm$n), ]
}

## E2b paper table: one-point delta-oracle comparison
# The oracle comparison is an asymptotic efficiency diagnostic, not a rate regression.  It is therefore reported at one pre-specified large grid point.
# With the default grid this is n=1000.  E2b_summary.csv still archives every n-value, including the important small-sample behavior.
E2_ORACLE_N <- max(E2_N)

E2b_oracle_table <- function(sm, n0 = E2_ORACLE_N) {
  required <- c("config", "n", "tmse_mle_median", "tmse_oracle_median",
                "tmse_ratio", "kl0_mle_median", "kl0_oracle_median",
                "kl_oracle_ratio", "oracle_benchmark")
  missing <- setdiff(required, names(sm))
  if (length(missing))
    stop("E2b summary is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(n0) != 1L || !is.finite(n0))
    stop("n0 must be one finite sample size.", call. = FALSE)
  if (!(n0 %in% sm$n))
    stop("Requested oracle grid point n0=", n0,
         " is not present. Available values: ",
         paste(sort(unique(sm$n)), collapse = ", "), call. = FALSE)

  ans <- sm[sm$n == n0,
            c("config", "n", "reps",
              "tmse_mle_median", "tmse_oracle_median", "tmse_ratio",
              "kl0_mle_median", "kl0_oracle_median", "kl_oracle_ratio",
              "oracle_benchmark", "delta_hat_iqr"),
            drop = FALSE]
  ans$kl_ratio_minus_benchmark <- ans$kl_oracle_ratio - ans$oracle_benchmark
  ans <- ans[order(ans$config), , drop = FALSE]
  rownames(ans) <- NULL
  ans
}

run_E2 <- function(..., oracle_n = E2_ORACLE_N) {
  # The simulation and per-grid summaries are unchanged from the original E2.
  a <- run_E2a(...); saveRDS(a, .res("E2a_raw.rds"))
  b <- run_E2b(...); saveRDS(b, .res("E2b_raw.rds"))
  sa <- summarise_E2a(a)
  sb <- summarise_E2b(b)

  # Paper-facing tables under the final three-table narrative.
  tab1 <- E2a_parameter_rate(sa)       # full grid; theory slope = -1/2
  tab2 <- E2a_kl_rate(sa)              # full grid; theory slope = -1
  tab3 <- E2b_oracle_table(sb, oracle_n) # one large-n point; benchmark = 1.4187

  write.csv(sa, .res("E2a_summary.csv"), row.names = FALSE)
  write.csv(sb, .res("E2b_summary.csv"), row.names = FALSE)

  # Export the three manuscript tables directly.
  write.csv(tab1, .res("E2a_parameter_rate.csv"), row.names = FALSE)
  write.csv(tab2, .res("E2a_kl_rate.csv"), row.names = FALSE)
  write.csv(tab3, .res("E2b_oracle.csv"), row.names = FALSE)

  list(
    E2a = sa,
    E2b = sb,
    table1_parameter_rate = tab1,
    table2_kl_rate = tab2,
    table3_oracle = tab3
  )
}
