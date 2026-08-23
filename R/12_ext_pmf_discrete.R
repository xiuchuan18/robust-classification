# 12_ext_pmf_discrete.R
# Discrete plug-in pmfs. Estimation procedures are those listed in supplement Table G.2.

.fit_llr_bernoulli <- function(x0, x1, alpha = 1) {
  p0 <- (sum(x0 == 1) + alpha) / (length(x0) + 2*alpha)
  p1 <- (sum(x1 == 1) + alpha) / (length(x1) + 2*alpha)
  l1 <- log(p1); l0 <- log(p0); c1 <- log1p(-p1); c0 <- log1p(-p0)
  function(x) ifelse(x == 1, l1 - l0, c1 - c0)
}

# Laplace-smoothed pmf over the pooled support.
.fit_llr_empirical <- function(x0, x1, alpha = 0.5) {
  lev <- sort(unique(c(x0, x1)))
  q0 <- (tabulate(match(x0, lev), length(lev)) + alpha) / (length(x0) + alpha*length(lev))
  q1 <- (tabulate(match(x1, lev), length(lev)) + alpha) / (length(x1) + alpha*length(lev))
  llr <- log(q1) - log(q0)
  function(x) { i <- match(x, lev); as.numeric(ifelse(is.na(i), 0, llr[i])) }
}

.fit_llr_poisson <- function(x0, x1) {
  l0 <- max(mean(x0), 1e-6); l1 <- max(mean(x1), 1e-6)
  function(x) x * (log(l1) - log(l0)) - (l1 - l0)   # = dpois(.,l1,log) - dpois(.,l0,log)
}

# Negative binomial, method-of-moments.
.fit_nb_par <- function(x) {
  m <- mean(x); v <- var(x)
  if (!is.finite(v) || v <= m) list(size = Inf, mu = max(m, 1e-6))
  else list(size = m^2/(v - m), mu = max(m, 1e-6))
}
.fit_llr_nbinom <- function(x0, x1) {
  a0 <- .fit_nb_par(x0); a1 <- .fit_nb_par(x1)
  lp <- function(x, a) if (is.infinite(a$size)) dpois(x, a$mu, log = TRUE)
                       else dnbinom(x, size = a$size, mu = a$mu, log = TRUE)
  function(x) pmax(lp(x, a1), LOG_FLOOR) - pmax(lp(x, a0), LOG_FLOOR)
}

.DISC_FITTERS <- list(bernoulli = .fit_llr_bernoulli, empirical = .fit_llr_empirical,
                      poisson = .fit_llr_poisson, nbinom = .fit_llr_nbinom)

detect_plugin_types <- function(X, max_levels = 10, count_family = "poisson") {
  apply(as.matrix(X), 2, function(v) {
    v <- v[is.finite(v)]; u <- unique(v)
    if (!all(abs(v - round(v)) < 1e-8)) return("tpsc")
    if (length(u) == 2) return("bernoulli")
    if (length(u) <= max_levels) return("empirical")
    if (min(v) >= 0) return(count_family)
    "tpsc"
  })
}

## one LLR closure for any family name, continuous or discrete
fit_llr_any <- function(x0, x1, fam) {
  if (!is.null(.DISC_FITTERS[[fam]])) .DISC_FITTERS[[fam]](x0, x1)
  else fit_llr_family(x0, x1, fam)
}

resolve_families <- function(X, policy) {
  p <- ncol(X)
  if (identical(policy, "auto")) return(detect_plugin_types(X))
  if (length(policy) == 1L) rep(policy, p) else policy
}

fit_mixed_plugin <- function(X, y, policy = "tpsc") {
  X <- as.matrix(X); y <- as.numeric(as.character(y))
  fams <- resolve_families(X, policy)
  cl <- sort(unique(y))
  fns <- lapply(seq_len(ncol(X)), function(j)
    fit_llr_any(X[y == cl[1], j], X[y == cl[2], j], fams[j]))
  list(fns = fns, family = fams, classes = cl, cols = colnames(X))
}

fit_srollr_any <- function(X, y, policy = "tpsc", foldid = NULL, nfolds = 5) {
  y  <- as.numeric(as.character(y))
  tr <- fit_mixed_plugin(X, y, policy)
  Tr <- apply_plugin_transform(X, tr)
  args <- list(x = Tr, y = y, family = "binomial", alpha = 1,
               type.measure = "deviance", standardize = FALSE)
  if (is.null(foldid)) args$nfolds <- nfolds else args$foldid <- foldid
  cvfit <- do.call(cv.glmnet, args)
  list(transform = tr, cvfit = cvfit, classes = sort(unique(y)))
}
