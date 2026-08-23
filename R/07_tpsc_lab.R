# 07_tpsc_lab.R --------------------------------------------------------------
# Shared infrastructure for the TPSC-estimation experiments E1-E2.
#
# Every object below concerns the Stage-1 problem only, i.e. estimating Psi = (theta, w, sigma, delta) from a
# single sample, and the quality of the induced log-likelihood-ratio transform T(x) = log f_TPSC(x | Psi_1) - log f_TPSC(x | Psi_0).
#
# Parameter order is (theta, w, sigma, delta) == (mode, weight, scale, df), identical to 01_core_tpsc.R and src/tpsc.cpp.
# Density convention (Fernandez & Steel, 1998), matching src/tpsc.cpp exactly:
#     f(x) = 2w      / sigma1 * g((x-theta)/sigma1; delta)   if x <  theta
#          = 2(1-w)  / sigma2 * g((x-theta)/sigma2; delta)   if x >= theta
#     sigma1 = sigma*sqrt(w/(1-w)),  sigma2 = sigma*sqrt((1-w)/w)
# with g(.; delta) the standard Student-t_delta density.

# Write outputs next to the project root (set by R/00b_setup_lab.R), so that results land in the same place regardless of the caller's working directory.
.res <- function(f) {
  root <- if (exists("ROOT", inherits = TRUE)) get("ROOT", inherits = TRUE) else getwd()
  dir.create(file.path(root, "results"), showWarnings = FALSE)
  file.path(root, "results", f)
}

## 1. Reference distribution functions

# Always computed on the log scale (dt(log = TRUE)) so that the far tails (|x| ~ 1e137, delta ~ 100) stay finite instead of underflowing to 0.
# NOTE: unlike src/tpsc.cpp this reference is NOT floored at log(1e-100), the floor is a deliberate implementation guard.
ldtpsc_c <- function(v, dtheta, w, sigma, delta) {
  stopifnot(w > 0, w < 1, sigma > 0, delta > 0)
  s1 <- sigma * sqrt(w / (1 - w))
  s2 <- sigma * sqrt((1 - w) / w)
  z  <- v - dtheta
  left <- z < 0
  out <- numeric(length(z))
  out[left]  <- log(2 * w)       - log(s1) + dt(z[left]  / s1, delta, log = TRUE)
  out[!left] <- log(2 * (1 - w)) - log(s2) + dt(z[!left] / s2, delta, log = TRUE)
  out
}

ldtpsc <- function(x, theta, w, sigma, delta) ldtpsc_c(x, theta, w, sigma, delta)
dtpsc  <- function(x, theta, w, sigma, delta) exp(ldtpsc(x, theta, w, sigma, delta))

# cdf:  F(x) = 2 w Pt(z1)                     , x <  theta
#            = w + (1-w) (2 Pt(z2) - 1)       , x >= theta
ptpsc <- function(q, theta, w, sigma, delta) {
  s1 <- sigma * sqrt(w / (1 - w))
  s2 <- sigma * sqrt((1 - w) / w)
  left <- q < theta
  out <- numeric(length(q))
  out[left]  <- 2 * w * pt((q[left] - theta) / s1, delta)
  out[!left] <- w + (1 - w) * (2 * pt((q[!left] - theta) / s2, delta) - 1)
  out
}

# rng: with prob w take theta - sigma1 |t_delta|, else theta + sigma2 |t_delta|.
rtpsc <- function(n, theta, w, sigma, delta) {
  stopifnot(w > 0, w < 1, sigma > 0, delta > 0)
  s1 <- sigma * sqrt(w / (1 - w))
  s2 <- sigma * sqrt((1 - w) / w)
  left <- runif(n) < w
  u    <- abs(rt(n, df = delta))
  ifelse(left, theta - s1 * u, theta + s2 * u)
}

# Half-scales and a single "effective width": a strongly skewed TPSC is wider on one side by a factor sqrt((1-w)/w), which `sigma` alone does not convey.
quad_scale <- function(P) {
  w <- min(max(P[2], 1e-12), 1 - 1e-12)
  max(P[3], 1e-300) * max(sqrt(w / (1 - w)), sqrt((1 - w) / w))
}

## 2. Deterministic quadrature on R (heavy-tail safe)
# Every integrand in E1/E2 is a TPSC density (or a function of two of them) times a smooth factor. The rule below is built for exactly that shape:
#   * the two outer half-lines are handled by an exp-sinh (double-exponential) rule,  v = scale * exp(pi/2 * sinh t),  which pushes the endpoint out
#     doubly exponentially and is therefore insensitive to how heavy the tail is (a Cauchy-matched  x = tan(u)  substitution, by contrast, leaves an
#     endpoint singularity of order (pi/2 - u)^(delta-1) and loses ~5 digits at delta = 0.5);
#   * the modes are used as PANEL BREAKPOINTS rather than being integrated across, so no rule ever straddles a kink; the finite panels between
#     consecutive modes get composite Simpson, where the integrand is smooth.
#
# T = 6 is the largest safe half-range: exp(pi/2 * sinh 7) overflows to Inf. At t = 6 the nodes already sit at |v| ~ 1e137 * scale, far past the tail of every TPSC member with delta >= 0.5.

.expsinh <- function(scale, m = 401L, T = 6) {
  t <- seq(-T, T, length.out = m)
  h <- t[2] - t[1]
  u <- scale * exp(pi / 2 * sinh(t))
  w <- h * u * (pi / 2) * cosh(t)
  bad <- !is.finite(u) | !is.finite(w)
  u[bad] <- scale; w[bad] <- 0
  list(u = u, w = w)
}

.simpson <- function(a, b, m = 401L) {
  m <- as.integer(m); if (m %% 2L == 0L) m <- m + 1L
  x <- seq(a, b, length.out = m); h <- x[2] - x[1]
  list(x = x, w = c(1, rep(c(4, 2), length.out = m - 2L), 1) * h / 3)
}

# Nodes/weights for integral over R  of an integrand with kinks at `breaks`. All features are in whatever frame the caller uses (see ldtpsc_c).
quad_nodes_breaks <- function(breaks, scale, m_tail = 401L, m_panel = 401L) {
  b  <- sort(unique(breaks[is.finite(breaks)]))
  if (!length(b)) b <- 0
  es <- .expsinh(scale, m_tail)
  x  <- c(b[1] - es$u, b[length(b)] + es$u)
  w  <- c(es$w, es$w)
  if (length(b) > 1L) {
    for (i in seq_len(length(b) - 1L)) {
      # keep the panel mesh fine relative to `scale` even if the modes are far apart
      mp <- min(4001L, max(m_panel, 2L * ceiling(20 * (b[i + 1] - b[i]) / scale) + 1L))
      sp <- .simpson(b[i], b[i + 1], mp)
      x  <- c(x, sp$x); w <- c(w, sp$w)
    }
  }
  list(x = x, w = w)
}

quad_nodes      <- function(center = 0, scale = 1, m = 401L)
  quad_nodes_breaks(center, scale, m, m)
quad_nodes_tpsc <- function(P, m = 401L) quad_nodes_breaks(P[1], quad_scale(P), m, m)

# Integrate a log-integrand safely:  sum(w * exp(lh)) with exp(-Inf) -> 0.
.quad_exp <- function(nd, lh) sum(nd$w * ifelse(is.finite(lh), exp(lh), 0))

## 3. Divergences and transform-space error
# All three work in the frame centred at a reference mode, so that neither a large |theta| nor a tiny sigma can cause cancellation.

# KL( f(.|P0) || f(.|P1) ) = E_{P0}[ log f0 - log f1 ].
kl_tpsc <- function(P0, P1, m = 401L) {
  ref <- P0[1]; d1 <- P1[1] - P0[1]
  nd  <- quad_nodes_breaks(c(0, d1), max(quad_scale(P0), quad_scale(P1)), m, m)
  l0  <- ldtpsc_c(nd$x, 0,  P0[2], P0[3], P0[4])
  l1  <- ldtpsc_c(nd$x, d1, P1[2], P1[3], P1[4])
  d   <- l0 - l1
  f0  <- ifelse(is.finite(l0), exp(l0), 0)
  sum(nd$w * f0 * ifelse(is.finite(d), d, 0))
}

# T(x; A, B) = log f(x|B) - log f(x|A)   (class 1 params B, class 0 params A).
Tfun <- function(x, PA, PB) {
  ldtpsc(x, PB[1], PB[2], PB[3], PB[4]) - ldtpsc(x, PA[1], PA[2], PA[3], PA[4])
}

# E_{X ~ (f0 + f1)/2} [ (That(X) - Tstar(X))^2 ], by the same rule.
# P0/P1 are the TRUE class-0/1 parameters, H0/H1 their estimates; the integrand
# has kinks at all four modes.
tmse_tpsc <- function(P0, P1, H0, H1, m = 401L) {
  ref <- 0.5 * (P0[1] + P1[1])
  dd  <- c(P0[1], P1[1], H0[1], H1[1]) - ref
  scl <- max(quad_scale(P0), quad_scale(P1), quad_scale(H0), quad_scale(H1),
             abs(P1[1] - P0[1]), 1e-300)
  nd  <- quad_nodes_breaks(dd, scl, m, m)
  l0  <- ldtpsc_c(nd$x, dd[1], P0[2], P0[3], P0[4])
  l1  <- ldtpsc_c(nd$x, dd[2], P1[2], P1[3], P1[4])
  h0  <- ldtpsc_c(nd$x, dd[3], H0[2], H0[3], H0[4])
  h1  <- ldtpsc_c(nd$x, dd[4], H1[2], H1[3], H1[4])
  dens <- 0.5 * (ifelse(is.finite(l0), exp(l0), 0) + ifelse(is.finite(l1), exp(l1), 0))
  err  <- (l1 - l0) - (h1 - h0)              # Tstar - That
  err[!is.finite(err)] <- 0
  sum(nd$w * dens * err^2)
}

## 4. Bounds, initialisations, boundary bookkeeping
# The bounds actually used by fit_tpsc_mle() in 01_core_tpsc.R.  Reported here explicitly because R1(4) asks for "the range of the parameters".
LAB_LOWER <- c(-Inf, 0.05, 1e-6, 0.5)     # (theta, w, sigma, delta)
LAB_UPPER <- c( Inf, 0.95,  Inf, 100)
# src/tpsc.cpp additionally returns +Inf outside
#   w in (1e-6, 1-1e-6), sigma > 1e-8, delta in (0.1, 1000).
CPP_HARD  <- list(w = c(1e-6, 1 - 1e-6), sigma = 1e-8, delta = c(0.1, 1000))

.robust_scale <- function(data) {
  md <- mad(data, constant = 1.4826)
  if (!is.finite(md) || md <= 1e-6) md <- max(sd(data), 1e-3)
  if (!is.finite(md) || md <= 1e-6) md <- 1e-3
  md
}

# M0's four deterministic starts (verbatim from fit_tpsc_mle()).
inits_m0 <- function(data) {
  med <- median(data); md <- .robust_scale(data)
  list(c(med, 0.3, md,        3),
       c(med, 0.7, md,       10),
       c(med, 0.5, md * 0.5,  2),
       c(med, 0.5, md * 1.5, 15))
}

# Latin hypercube starts on (theta, w, log sigma, log delta).
.lhs <- function(n, k) vapply(seq_len(k), function(j) (sample(n) - runif(n)) / n,
                              numeric(n))

# M1's Latin-hypercube multi-start (near-global reference).
inits_m1 <- function(data, n_start = 50L, lower = LAB_LOWER, upper = LAB_UPPER) {
  md <- .robust_scale(data)
  q  <- quantile(data, c(0.1, 0.9), names = FALSE, type = 7)
  U  <- .lhs(n_start, 4L)
  wlo <- max(lower[2], 0.10); whi <- min(upper[2], 0.90)
  dlo <- max(lower[4], 1e-3); dhi <- min(upper[4], 60)
  lapply(seq_len(n_start), function(i)
    c(q[1] + U[i, 1] * (q[2] - q[1]),
      wlo  + U[i, 2] * (whi - wlo),
      md * exp(-1.5 + 3 * U[i, 3]),
      exp(log(dlo) + U[i, 4] * (log(dhi) - log(dlo)))))
}

boundary_flags <- function(par, lower = LAB_LOWER, upper = LAB_UPPER, tol = 1e-6) {
  c(w_lo     = as.integer(par[2] <= lower[2] + tol),
    w_up     = as.integer(par[2] >= upper[2] - tol),
    sigma_lo = as.integer(par[3] <= lower[3] * (1 + 1e-3)),
    delta_lo = as.integer(par[4] <= lower[4] + tol),
    delta_up = as.integer(par[4] >= upper[4] * (1 - 1e-6)))
}

## 6. Objective wrappers
SAFE_BIG <- 1e12   # finite stand-in for +Inf so that L-BFGS-B never errors

.nll_safe <- function(params, data) {
  v <- tpsc_nll_cpp(params, data)
  if (!is.finite(v)) SAFE_BIG else v
}

## 7. Timing
.now <- function() as.numeric(Sys.time())

# All fitters return a named list:
#   par        numeric(4)
#   nll        UNPENALISED negative log-likelihood at par (comparable across methods)
#   n_starts   number of starts attempted
#   n_err      number of starts that threw
#   fallback   TRUE if no start produced a finite value and a canned Psi was returned
#   elapsed    wall-clock seconds
#   lower/upper  THE BOX ACTUALLY USED.  Callers must read the box from here and
#              never assume LAB_LOWER/LAB_UPPER -- E1b fits under six different
#              boxes, and judging a fit against the wrong one silently
#              mis-reports boundary hits.

.fit_core <- function(data, inits, lower, upper, fn, method = "L-BFGS-B",
                      maxit = 500) {
  stopifnot(length(lower) == 4L, length(upper) == 4L, all(lower <= upper))
  best <- NULL; best_val <- Inf; n_err <- 0L
  for (init in inits) {
    init <- pmin(pmax(init, lower), upper)
    res <- tryCatch(
      optim(par = init, fn = fn, data = data, method = method,
            lower = lower, upper = upper,
            control = list(maxit = maxit, pgtol = 1e-8)),
      error = function(e) NULL)
    if (is.null(res)) { n_err <- n_err + 1L; next }
    if (is.finite(res$value) && res$value < best_val) {
      best_val <- res$value; best <- res
    }
  }
  list(best = best, n_err = n_err)
}

.canned <- function(data) {
  c(median(data), 0.5, .robust_scale(data), 5)
}

.wrap <- function(data, inits, lower, upper, fn, method = "L-BFGS-B") {
  t0   <- .now()
  data <- data[is.finite(data)]        # same guard as fit_tpsc_mle()
  r    <- .fit_core(data, inits, lower, upper, fn, method)
  fb   <- is.null(r$best)
  par  <- if (fb) .canned(data) else as.numeric(r$best$par)
  list(par = par, nll = .nll_safe(par, data), n_starts = length(inits),
       n_err = r$n_err, fallback = fb, elapsed = .now() - t0,
       lower = lower, upper = upper)
}

# M0: the shipped estimator (4 deterministic starts, L-BFGS-B)
fit_m0 <- function(data, lower = LAB_LOWER, upper = LAB_UPPER) {
  t0   <- .now()
  data <- data[is.finite(data)]
  if (length(data) < 5)
    return(list(par = as.numeric(c(median(data), 0.5, max(sd(data), 1e-3), 5)),
                nll = NA_real_, n_starts = 0L, n_err = 0L, fallback = TRUE,
                elapsed = .now() - t0, lower = lower, upper = upper))
  .wrap(data, inits_m0(data), lower, upper, tpsc_nll_cpp)
}

# M1: 50 Latin-hypercube starts
fit_m1 <- function(data, lower = LAB_LOWER, upper = LAB_UPPER, n_start = 50L)
  .wrap(data, inits_m1(data, n_start, lower, upper), lower, upper, tpsc_nll_cpp)

# M2: Nelder-Mead on an unconstrained reparametrisation
# theta = a1; w = wlo + (wup-wlo) plogis(a2); sigma = exp(a3);
# delta = dlo + (dup-dlo) plogis(a4).
.nm_to_par <- function(a, lower, upper) {
  dup <- if (is.finite(upper[4])) upper[4] else 100
  c(a[1],
    lower[2] + (upper[2] - lower[2]) * plogis(a[2]),
    exp(a[3]),
    lower[4] + (dup - lower[4]) * plogis(a[4]))
}
.par_to_nm <- function(p, lower, upper) {
  dup <- if (is.finite(upper[4])) upper[4] else 100
  c(p[1],
    qlogis(min(max((p[2] - lower[2]) / (upper[2] - lower[2]), 1e-4), 1 - 1e-4)),
    log(max(p[3], 1e-8)),
    qlogis(min(max((p[4] - lower[4]) / (dup - lower[4]), 1e-4), 1 - 1e-4)))
}
fit_m2 <- function(data, lower = LAB_LOWER, upper = LAB_UPPER) {
  t0 <- .now()
  obj <- function(a, data) .nll_safe(.nm_to_par(a, lower, upper), data)
  best <- NULL; best_val <- Inf; n_err <- 0L
  for (init in inits_m0(data)) {
    a0 <- .par_to_nm(pmin(pmax(init, lower), upper), lower, upper)
    res <- tryCatch(optim(a0, obj, data = data, method = "Nelder-Mead",
                          control = list(maxit = 2000, reltol = 1e-10)),
                    error = function(e) NULL)
    if (is.null(res)) { n_err <- n_err + 1L; next }
    if (is.finite(res$value) && res$value < best_val) { best_val <- res$value; best <- res }
  }
  fb  <- is.null(best)
  par <- if (fb) .canned(data) else .nm_to_par(best$par, lower, upper)
  list(par = par, nll = .nll_safe(par, data), n_starts = 4L, n_err = n_err,
       fallback = fb, elapsed = .now() - t0, lower = lower, upper = upper)
}

# delta-oracle: delta fixed at its true value, only (theta, w, sigma) fit
fit_oracle_delta <- function(data, delta0, lower = LAB_LOWER, upper = LAB_UPPER) {
  t0  <- .now()
  obj <- function(p3, data) .nll_safe(c(p3, delta0), data)
  best <- NULL; best_val <- Inf; n_err <- 0L
  for (init in inits_m0(data)) {
    res <- tryCatch(optim(par = init[1:3], fn = obj, data = data, method = "L-BFGS-B",
                          lower = lower[1:3], upper = upper[1:3],
                          control = list(maxit = 500, pgtol = 1e-8)),
                    error = function(e) NULL)
    if (is.null(res)) { n_err <- n_err + 1L; next }
    if (is.finite(res$value) && res$value < best_val) { best_val <- res$value; best <- res }
  }
  fb  <- is.null(best)
  par <- if (fb) .canned(data) else c(as.numeric(best$par), delta0)
  list(par = par, nll = .nll_safe(par, data), n_starts = 4L, n_err = n_err,
       fallback = fb, elapsed = .now() - t0, lower = lower, upper = upper)
}


FITTERS <- list(
  "M0: 4 deterministic starts (shipped)" = fit_m0,
  "M1: 50 LHS random starts"             = fit_m1,
  "M2: Nelder-Mead (reparametrised)"     = fit_m2
)

.key_seed <- function(...) {
  z <- utf8ToInt(paste(c(...), collapse = "|"))
  h <- 0
  for (i in seq_along(z)) h <- (h * 131 + z[i]) %% 2147483647
  as.integer(h)
}

E_CONFIGS <- list(
  "S1 near-Gaussian symmetric" = c(theta = 0, w = 0.50, sigma = 1, delta = 30),
  "S2 symmetric heavy tail"    = c(theta = 0, w = 0.50, sigma = 1, delta = 2),
  "S3 skewed moderate tail"    = c(theta = 0, w = 0.30, sigma = 1, delta = 5),
  "S4 skewed Cauchy tail"      = c(theta = 0, w = 0.25, sigma = 1, delta = 1)
)

MODE_SHIFT <- 0.5
class1_par <- function(P0) { P1 <- P0; P1[1] <- P0[1] + MODE_SHIFT * P0[3]; P1 }
