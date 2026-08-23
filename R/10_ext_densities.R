# 10_ext_densities.R
# The S-RoLLR framework is plug-in in its first stage: any per-feature,
# per-class working density (or pmf) f_{cj} yields the transform
#
#     T_ij = log f_{1j}(x_ij) - log f_{0j}(x_ij),
#
# which the second stage feeds to L1-penalised logistic regression.
#
# Depends on R/01_core_tpsc.R (fit_tpsc_mle, d_tpsc).
#
# Manuscript: Section 6, working densities of supplement Table G.2 (Gaussian, Student-t, skew-normal, skew-t, TPSC, BIC mixture) plus the `identity` reference (L1-RL).

LOG_FLOOR <- log(1e-100)          # same density floor as the Rcpp core

.safe_log <- function(v) pmax(ifelse(is.finite(v), v, LOG_FLOOR), LOG_FLOOR)

## 1. Gaussian
.pfit_gaussian <- function(x) list(m = mean(x), s = max(sd(x), 1e-6))
.plog_gaussian <- function(x, a) dnorm(x, a$m, a$s, log = TRUE)

## 2. symmetric Student-t (location-scale)
.pfit_t <- function(x) {
  nll <- function(par) {
    mu <- par[1]; ls <- par[2]; ldf <- par[3]
    s <- exp(ls); df <- exp(ldf)
    if (!is.finite(s) || !is.finite(df) || df > 1000) return(1e10)
    -sum(dt((x - mu) / s, df, log = TRUE) - ls)
  }
  md <- mad(x, constant = 1.4826); if (!is.finite(md) || md <= 1e-6) md <- max(sd(x), 1e-3)
  best <- NULL; bv <- Inf
  for (d0 in c(3, 10)) {
    r <- tryCatch(optim(c(median(x), log(md), log(d0)), nll,
                        method = "Nelder-Mead", control = list(maxit = 500)),
                  error = function(e) NULL)
    if (!is.null(r) && is.finite(r$value) && r$value < bv) { bv <- r$value; best <- r }
  }
  if (is.null(best)) list(m = median(x), s = md, df = 5)
  else list(m = best$par[1], s = exp(best$par[2]), df = min(exp(best$par[3]), 1000))
}
.plog_t <- function(x, a) dt((x - a$m) / a$s, a$df, log = TRUE) - log(a$s)

## 3. skew-normal (Azzalini)
.pfit_snorm <- function(x) {
  r <- tryCatch(sn::selm.fit(x = matrix(1, length(x), 1), y = x, family = "SN")$param$dp,
                error = function(e) NULL)
  if (is.null(r) || !all(is.finite(r)))
    r <- c(mean(x), max(sd(x), 1e-6), 0)
  list(xi = r[1], om = max(r[2], 1e-6), al = max(min(r[3], 50), -50))
}
.plog_snorm <- function(x, a) sn::dsn(x, a$xi, a$om, a$al, log = TRUE)

## 4. skew-t (Azzalini-Capitanio)
.pfit_st <- function(x) {
  r <- tryCatch(sn::st.mple(y = x)$dp, error = function(e) NULL)
  if (is.null(r) || !all(is.finite(r)) || r[4] <= 0)
    r <- c(mean(x), max(sd(x), 1e-6), 0, 10)
  list(xi = r[1], om = max(r[2], 1e-6), al = max(min(r[3], 50), -50),
       nu = max(min(r[4], 1e3), 0.3))
}
.plog_st <- function(x, a) sn::dst(x, a$xi, a$om, a$al, a$nu, log = TRUE)

## 5. finite location-scale mixture, K and nu chosen by BIC
MIX_K_GRID  <- 1:3
MIX_DF_GRID <- c(5, Inf)

.mix_ldens <- function(x, m, v, df) {
  if (is.infinite(df)) dnorm(x, m, sqrt(v), log = TRUE)
  else dt((x - m) / sqrt(v), df, log = TRUE) - 0.5 * log(v)
}

# n*K matrix of log(pr_k) + log f_k(x)
.mix_lmat <- function(x, a) {
  K <- length(a$pr)
  matrix(vapply(seq_len(K), function(k)
    log(a$pr[k]) + .mix_ldens(x, a$mu[k], a$vr[k], a$df), numeric(length(x))),
    nrow = length(x), ncol = K)
}
.mix_lse <- function(L) {
  mx <- apply(L, 1, max)
  mx + log(rowSums(exp(L - mx)))
}

# EM for a K-component mixture at fixed df, returns the fit, its log-likelihood and its free-parameter count (for BIC).
.pfit_mixK <- function(x, K, df, iter = 80, tol = 1e-7) {
  n <- length(x); sx <- max(sd(x), 1e-3); vfloor <- (sx / 50)^2
  mu <- as.numeric(quantile(x, seq_len(K) / (K + 1), names = FALSE))
  if (K > 1 && min(diff(sort(mu))) < 1e-8)  
    mu <- mu + sx * 0.2 * (seq_len(K) - (K + 1) / 2)
  a <- list(pr = rep(1/K, K), mu = mu,
            vr = rep(max(var(x) / K, vfloor), K), df = df)
  old <- -Inf; ll <- -Inf
  for (it in seq_len(iter)) {
    L <- .mix_lmat(x, a); lse <- .mix_lse(L); ll <- sum(lse)
    if (!is.finite(ll)) break
    W  <- exp(L - lse)
    # the t E-step also needs the latent scale weights u_ik
    U <- if (is.infinite(df)) matrix(1, n, K) else
      matrix(vapply(seq_len(K), function(k)
        (df + 1) / (df + (x - a$mu[k])^2 / a$vr[k]), numeric(n)), n, K)
    nk <- colSums(W)
    a$pr <- pmax(nk / n, 1e-3); a$pr <- a$pr / sum(a$pr)
    for (k in seq_len(K)) {
      wk <- W[, k] * U[, k]; sw <- max(sum(wk), 1e-8)
      a$mu[k] <- sum(wk * x) / sw
      a$vr[k] <- max(sum(wk * (x - a$mu[k])^2) / max(nk[k], 1e-8), vfloor)
    }
    if (abs(ll - old) < tol * (abs(old) + 1)) break
    old <- ll
  }
  a$ll <- ll; a$npar <- 3 * K - 1        # K means, K variances, K-1 weights
  a
}

.pfit_mixture <- function(x) {
  n <- length(x); best <- NULL; best_bic <- Inf
  for (K in MIX_K_GRID) for (df in MIX_DF_GRID) {
    a <- tryCatch(.pfit_mixK(x, K, df), error = function(e) NULL)
    if (is.null(a) || !is.finite(a$ll)) next
    bic <- -2 * a$ll + a$npar * log(n)
    if (bic < best_bic) { best_bic <- bic; best <- a }
  }
  if (is.null(best))
    best <- list(pr = 1, mu = mean(x), vr = max(var(x), 1e-6), df = Inf,
                 ll = NA_real_, npar = 2)
  best$bic <- best_bic
  best
}

.plog_mixture <- function(x, a) .mix_lse(.mix_lmat(x, a))

## TPSC
.pfit_tpsc <- function(x)
  if (var(x) < 1e-10) c(mean(x), 0.5, 1e-3, 5) else fit_tpsc_mle(x)
.plog_tpsc <- function(x, a)
  .safe_log(log(pmax(d_tpsc(x, a[1], a[2], a[3], a[4]), 1e-100)))


PLUGIN_FAMILIES <- list(
  gaussian = list(fit = .pfit_gaussian, logf = .plog_gaussian),
  studentt = list(fit = .pfit_t,        logf = .plog_t),
  snorm    = list(fit = .pfit_snorm,    logf = .plog_snorm),
  skewt    = list(fit = .pfit_st,       logf = .plog_st),
  tpsc     = list(fit = .pfit_tpsc,     logf = .plog_tpsc),
  mixture  = list(fit = .pfit_mixture,  logf = .plog_mixture)
)

CONT_FAMILY_ORDER <- c("gaussian", "studentt", "snorm", "skewt", "tpsc", "mixture")

fit_llr_family <- function(x0, x1, fam) {
  if (identical(fam, "identity")) return(function(x) x)
  spec <- PLUGIN_FAMILIES[[fam]]
  if (is.null(spec)) stop("unknown plug-in family: ", fam)
  a0 <- spec$fit(x0); a1 <- spec$fit(x1)
  function(x) .safe_log(spec$logf(x, a1)) - .safe_log(spec$logf(x, a0))
}

## S-RoLLR with an arbitrary stage-1 working density.
fit_plugin_transform <- function(X, y, family) {
  X <- as.matrix(X); y <- as.numeric(as.character(y)); p <- ncol(X)
  if (length(family) == 1L) family <- rep(family, p)
  stopifnot(length(family) == p)
  cl <- sort(unique(y))
  fns <- lapply(seq_len(p), function(j)
    fit_llr_family(X[y == cl[1], j], X[y == cl[2], j], family[j]))
  list(fns = fns, family = family, classes = cl, cols = colnames(X))
}

apply_plugin_transform <- function(X, tr) {
  X <- as.matrix(X)
  Tr <- vapply(seq_len(ncol(X)), function(j) tr$fns[[j]](X[, j]), numeric(nrow(X)))
  if (is.null(dim(Tr))) Tr <- matrix(Tr, nrow = nrow(X))
  Tr[!is.finite(Tr)] <- 0
  colnames(Tr) <- if (!is.null(tr$cols)) tr$cols else paste0("V", seq_len(ncol(X)))
  Tr
}

fit_srollr_plugin <- function(X, y, family = "tpsc", nfolds = 5,
                              glmnet_args = list()) {
  y  <- as.numeric(as.character(y))
  tr <- fit_plugin_transform(X, y, family)
  Tr <- apply_plugin_transform(X, tr)
  cv_args <- modifyList(
    list(x = Tr, y = y, family = "binomial", alpha = 1, nfolds = nfolds,
         type.measure = "deviance", standardize = FALSE), glmnet_args)
  cvfit <- do.call(cv.glmnet, cv_args)
  list(transform = tr, cvfit = cvfit, classes = sort(unique(y)))
}

predict_srollr_plugin <- function(model, X_new, s = "lambda.min") {
  Tr <- apply_plugin_transform(X_new, model$transform)
  as.numeric(predict(model$cvfit, newx = Tr, s = s, type = "response"))
}

srollr_plugin_selected <- function(model, s = "lambda.min") {
  co <- as.matrix(coef(model$cvfit, s = s))
  setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
}
