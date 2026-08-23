# 02_competitors.R
# Competitor classifiers implemented in-house:
#   MOKE   : mode-based unimodal classifier (Xiong et al., 2025), vectorized.
#   GMM-NB : Gaussian-mixture naive Bayes (mclust).
# The other competitors (L1-RL, SVM, RF, LightGBM) are called directly from their C/Fortran packages in 06_fit_all.R.

## MOKE
.kde_mode <- function(x, bw = NULL) {
  if (length(x) < 2) return(mean(x))
  d <- if (is.null(bw)) density(x) else density(x, bw = bw)
  d$x[which.max(d$y)]
}
.rot_bw <- function(x) {
  n <- length(x); if (n < 2) return(1e-3)
  s <- sd(x, na.rm = TRUE); iqr <- IQR(x, na.rm = TRUE)
  if (is.na(s) || s == 0) s <- 1e-3
  if (is.na(iqr) || iqr == 0) iqr <- 1.34 * s
  bw <- 1.06 * min(s, iqr / 1.34, na.rm = TRUE) * n^(-0.2)
  if (is.na(bw) || bw == 0) 0.1 * n^(-0.2) else bw
}
.scale_fit   <- function(X) {
  ctr <- colMeans(X, na.rm = TRUE)
  scl <- apply(X, 2, sd, na.rm = TRUE); scl[scl == 0 | is.na(scl)] <- 1
  list(center = ctr, scale = scl)
}
.scale_apply <- function(X, sc) sweep(sweep(X, 2, sc$center, "-"), 2, sc$scale, "/")

fit_moke <- function(X, y, scale_data = TRUE) {
  X <- as.matrix(X); y <- as.factor(y); classes <- levels(y)
  sc <- if (scale_data) .scale_fit(X) else NULL
  Xs <- if (scale_data) .scale_apply(X, sc) else X

  # per-feature bandwidths from class-1 observations
  sigma1 <- vapply(seq_len(ncol(Xs)),
                   function(j) .rot_bw(Xs[y == classes[1], j]), numeric(1))
  modes <- lapply(classes, function(cls)
    vapply(seq_len(ncol(Xs)),
           function(j) .kde_mode(Xs[y == cls, j], .rot_bw(Xs[y == cls, j])),
           numeric(1)))
  names(modes) <- classes

  # choose theta = sigma2/sigma1 maximizing training accuracy (vectorized grid)
  A <- sweep(Xs, 2, modes[[classes[1]]], "-"); A <- sweep(A, 2, sigma1, "/")
  psi1 <- rowSums(dnorm(A))                                # fixed across theta
  B <- sweep(Xs, 2, modes[[classes[2]]], "-")             # /(theta*sigma1) later
  Bs <- sweep(B, 2, sigma1, "/")
  theta_grid <- seq(1/3, 3, by = 0.05)
  ycl <- as.character(y)
  acc <- vapply(theta_grid, function(th) {
    psi2 <- rowSums(dnorm(Bs / th))
    mean(ifelse(psi1 - psi2 > 0, classes[1], classes[2]) == ycl)
  }, numeric(1))
  theta <- theta_grid[which.max(acc)]

  list(classes = classes, scale = sc, scale_data = scale_data,
       modes = modes, sigma1 = sigma1, sigma2 = theta * sigma1)
}

.moke_scores <- function(model, X_test) {
  Xs <- if (model$scale_data) .scale_apply(as.matrix(X_test), model$scale) else as.matrix(X_test)
  A <- sweep(Xs, 2, model$modes[[model$classes[1]]], "-"); A <- sweep(A, 2, model$sigma1, "/")
  Bd <- sweep(Xs, 2, model$modes[[model$classes[2]]], "-"); Bd <- sweep(Bd, 2, model$sigma2, "/")
  rowSums(dnorm(A)) - rowSums(dnorm(Bd))                   # psi1 - psi2
}
predict_moke      <- function(model, X_test)
  ifelse(.moke_scores(model, X_test) > 0, model$classes[1], model$classes[2])
predict_moke_prob <- function(model, X_test) {            # P(class 2)
  1 / (1 + exp(.moke_scores(model, X_test)))
}

## GMM-NB
# `G` is the candidate number of mixture components per feature, selected by BIC (mclust). The default 1:9 is mclust's own default. It is exposed as an explicit
# argument so the configuration is transparent and reproducible: capping it (e.g. G = 1:3) is a modelling choice that speeds up the BIC search substantially and should be justified on statistical grounds, not used silently for speed.
fit_gmm_nb <- function(X, y, G = 1:9) {
  X <- as.matrix(X); classes <- sort(unique(y))
  gmm <- lapply(classes, function(cls) {
    Xc <- X[y == cls, , drop = FALSE]
    lapply(seq_len(ncol(X)),
           function(j) mclust::Mclust(Xc[, j], G = G, verbose = FALSE))
  })
  names(gmm) <- as.character(classes)
  list(classes = classes, gmm = gmm)
}

# vectorized prediction: one dens() call per feature over all test rows.
predict_gmm_nb_prob <- function(model, X_test) {          # P(class 2)
  X_test <- as.matrix(X_test); cl <- as.character(model$classes)
  loglik <- sapply(cl, function(c) {
    LL <- vapply(seq_len(ncol(X_test)), function(j) {
      g <- model$gmm[[c]][[j]]
      log(pmax(mclust::dens(modelName = g$modelName, data = X_test[, j],
                            parameters = g$parameters), 1e-100))
    }, numeric(nrow(X_test)))
    rowSums(LL)
  })
  m   <- apply(loglik, 1, max)
  bad <- !is.finite(m)
  e   <- exp(loglik - m); p2 <- (e / rowSums(e))[, 2]
  p2[bad] <- 0.5
  p2
}
