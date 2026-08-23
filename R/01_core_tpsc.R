# 01_core_tpsc.R
# The TPSC log-likelihood-ratio classifier.
# Density / likelihood / transform are evaluated by the Rcpp core in
# src/tpsc.cpp:  d_tpsc_cpp, tpsc_nll_cpp, tpsc_transform_cpp, tpsc_rowloglik_cpp.

## TPSC density (R wrapper, identical to the Rcpp density)
d_tpsc <- function(x, theta, w, sigma, delta) {
  d_tpsc_cpp(x, w, theta, sigma, delta)
}

## per-feature maximum-likelihood fit
TPSC_LOWER <- c(-Inf, 0.05, 1e-6, 0.5)
TPSC_UPPER <- c( Inf, 0.95, Inf, 100)

fit_tpsc_mle <- function(data) {
  data <- data[is.finite(data)]
  if (length(data) < 5)
    return(as.numeric(c(median(data), 0.5, max(sd(data), 1e-3), 5)))

  med <- median(data)
  md  <- mad(data, constant = 1.4826)
  if (!is.finite(md) || md <= 1e-6) md <- max(sd(data), 1e-3)

  inits <- list(c(med, 0.3, md,       3),
                c(med, 0.7, md,      10),
                c(med, 0.5, md * 0.5, 2),
                c(med, 0.5, md * 1.5,15))

  best <- NULL; best_val <- Inf
  for (init in inits) {
    res <- tryCatch(
      optim(par = init, fn = tpsc_nll_cpp, data = data,
            method = "L-BFGS-B", lower = TPSC_LOWER, upper = TPSC_UPPER,
            control = list(maxit = 500, pgtol = 1e-8)),
      error = function(e) NULL)
    if (!is.null(res) && is.finite(res$value) && res$value < best_val) {
      best_val <- res$value; best <- res
    }
  }
  if (is.null(best)) as.numeric(c(med, 0.5, md, 5)) else as.numeric(best$par)
}

## fit class-conditional TPSC parameters for every feature
# Returns: classes (sorted), and a p * 4 parameter matrix per class.
fit_tpsc_params <- function(X, y) {
  X <- as.matrix(X)
  classes <- sort(unique(y))
  p <- ncol(X)
  params <- lapply(classes, function(cls) {
    Xc <- X[y == cls, , drop = FALSE]
    P  <- matrix(NA_real_, nrow = p, ncol = 4)
    for (j in seq_len(p)) {
      xj <- Xc[, j]
      if (var(xj, na.rm = TRUE) < 1e-10)
        P[j, ] <- c(mean(xj, na.rm = TRUE), 0.5, 1e-3, 5)
      else
        P[j, ] <- fit_tpsc_mle(xj)
    }
    P
  })
  names(params) <- as.character(classes)
  list(classes = classes, params = params)
}

## RoLLR: hard classifier
fit_rollr <- function(X, y) {
  y <- as.numeric(as.character(y))
  fit_tpsc_params(X, y)
}

# class-probability version (softmax of the two class log-likelihoods), and the 0.5 threshold reproduces the hard rule sum_j T_j > 0.
predict_rollr_prob <- function(model, X_new) {
  X_new <- as.matrix(X_new)
  cl <- as.character(model$classes)
  ll0 <- tpsc_rowloglik_cpp(X_new, model$params[[cl[1]]])
  ll1 <- tpsc_rowloglik_cpp(X_new, model$params[[cl[2]]])
  m   <- pmax(ll0, ll1)
  p1  <- exp(ll1 - m) / (exp(ll0 - m) + exp(ll1 - m))
  bad <- !is.finite(ll0) | !is.finite(ll1)
  p1[bad] <- 0.5
  p1
}

predict_rollr <- function(model, X_new) {
  p1 <- predict_rollr_prob(model, X_new)
  ifelse(p1 >= 0.5, model$classes[2], model$classes[1])
}

## feature transform  T_ij = log f_{1j}(x) - log f_{0j}(x)
tpsc_transform <- function(X, model) {
  X  <- as.matrix(X)
  cl <- as.character(model$classes)
  Tr <- tpsc_transform_cpp(X, model$params[[cl[1]]], model$params[[cl[2]]])
  colnames(Tr) <- colnames(X)
  Tr
}

## S-RoLLR: L1 logistic on the transformed features
fit_srollr <- function(X, y, nfolds = 5, parallel = FALSE, glmnet_args = list()) {
  y  <- as.numeric(as.character(y))
  tp <- fit_tpsc_params(X, y)
  Tr <- tpsc_transform(X, tp)

  cv_args <- modifyList(
    list(x = Tr, y = y, family = "binomial", alpha = 1,
         nfolds = nfolds, type.measure = "deviance", standardize = FALSE),
    glmnet_args)
  cv_args$parallel <- parallel
  cvfit <- do.call(cv.glmnet, cv_args)

  list(tpsc = tp, cvfit = cvfit, classes = sort(unique(y)))
}

predict_srollr <- function(model, X_new, s = "lambda.min") {
  Tr <- tpsc_transform(X_new, model$tpsc)
  as.numeric(predict(model$cvfit, newx = Tr, s = s, type = "response"))
}

# nonzero (selected) features of a fitted S-RoLLR model
srollr_selected <- function(model, s = "lambda.min") {
  co <- as.matrix(coef(model$cvfit, s = s))
  setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
}
