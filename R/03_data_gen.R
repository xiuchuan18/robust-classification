# 03_data_gen.R
# Twelve data-generation mechanisms for the simulation study of Section 4;
# the settings themselves are tabulated in supplement Table E.1.

.ar1_corr <- function(size, rho) rho^abs(outer(seq_len(size), seq_len(size), "-"))

## helpers for the mean-matched shape-signal Cases 10-12
# Half-t mean E|T_delta| (delta > 1); used to centre skewed TPSC draws to mean 0.
.halft_mean <- function(delta) sqrt(delta/pi) * gamma((delta-1)/2) / gamma(delta/2)
# E[X] - theta for a raw two-piece-t (theta, w, sigma, delta) sample.
.tpsc_meanshift <- function(w, sigma, delta) {
  s1 <- sigma*sqrt(w/(1-w)); s2 <- sigma*sqrt((1-w)/w)
  .halft_mean(delta) * ((1-w)*s2 - w*s1)
}

# Returns list(X, y, n_signal) where n_signal = number of informative features.
generate_data <- function(scenario, n, p, p_relevant) {
  X <- matrix(rnorm(n * p), n, p)
  n_signal <- p_relevant

  if (scenario == "Case1: Normal") {
    y <- rbinom(n, 1, 0.5); m <- 0.5; sdv <- 2
    for (j in 1:p_relevant) {
      X[y == 0, j] <- rnorm(sum(y == 0), -m, sdv)
      X[y == 1, j] <- rnorm(sum(y == 1),  m, sdv)
    }
  } else if (scenario == "Case2: Undefined Moments") {
    y <- rbinom(n, 1, 0.5)
    for (j in 1:p_relevant) {
      X[y == 0, j] <- rcauchy(sum(y == 0), 0,   1)
      X[y == 1, j] <- rcauchy(sum(y == 1), 0.5, 1)
    }
  } else if (scenario == "Case3: Skewed LogNormal") {
    y <- rbinom(n, 1, 0.5)
    for (j in 1:p_relevant) {
      X[y == 0, j] <- rlnorm(sum(y == 0), 0,   1)
      X[y == 1, j] <- rlnorm(sum(y == 1), 0.5, 1)
    }
  } else if (scenario == "Case4: Normal vs. HeavyTail") {
    y <- rbinom(n, 1, 0.5)
    for (j in 1:p_relevant) {
      X[y == 0, j] <- rnorm(sum(y == 0), 0, 1)
      X[y == 1, j] <- rt(sum(y == 1), df = 2) + 0.5
    }
  } else if (scenario == "Case5: Outliers") {
    y <- rbinom(n, 1, 0.5); orate <- 0.1
    for (j in 1:p_relevant) {
      mu0 <- runif(1, -0.5, 0); mu1 <- runif(1, 0, 0.5); sc <- runif(1, 0.5, 1.5)
      n0 <- sum(y == 0); n0c <- round(n0 * (1 - orate))
      X[y == 0, j] <- sample(c(rnorm(n0c, mu0, sc), rnorm(n0 - n0c, mu0, sc * 10)))
      n1 <- sum(y == 1); n1c <- round(n1 * (1 - orate))
      X[y == 1, j] <- sample(c(rnorm(n1c, mu1, sc), rnorm(n1 - n1c, mu1, sc * 10)))
    }
  } else if (scenario == "Case6: Unequal Covariance") {
    y <- sort(rbinom(n, 1, 0.5)); n0 <- sum(y == 0); n1 <- sum(y == 1)
    pn <- p - p_relevant
    cs <- .ar1_corr(p_relevant, 0.2)
    S0 <- as.matrix(Matrix::bdiag(cs,     diag(pn)))
    S1 <- as.matrix(Matrix::bdiag(cs * 2, diag(pn)))
    mu0 <- c(rep(0,   p_relevant), rep(0, pn))
    mu1 <- c(rep(0.3, p_relevant), rep(0, pn))
    X <- rbind(MASS::mvrnorm(n0, mu0, S0), MASS::mvrnorm(n1, mu1, S1))
    idx <- sample(n); X <- X[idx, ]; y <- y[idx]
  } else if (scenario == "Case7: Signal Correlation") {
    y <- sort(rbinom(n, 1, 0.5)); n0 <- sum(y == 0); n1 <- sum(y == 1)
    pn <- p - p_relevant
    cs <- .ar1_corr(p_relevant, 0.8)
    S0 <- as.matrix(Matrix::bdiag(cs, diag(pn))); S1 <- S0
    mu0 <- c(rep( 0.5, p_relevant), rep(0, pn))
    mu1 <- c(rep(-0.5, p_relevant), rep(0, pn))
    X <- rbind(MASS::mvrnorm(n0, mu0, S0), MASS::mvrnorm(n1, mu1, S1))
    idx <- sample(n); X <- X[idx, ]; y <- y[idx]
  } else if (scenario == "Case8: Nonlinear") {
    f <- 0
    for (j in 1:8)   f <- f + 0.4 * (X[, j]^2 - 1)
    for (j in 9:14)  f <- f + 0.5 * X[, j] * X[, j + 6]
    for (j in 15:20) f <- f + 0.5 * sin(3 * pi * (X[, j - 1] + X[, j]))
    y <- rbinom(n, 1, pnorm(f))
    n_signal <- 20
  } else if (scenario == "Case9: HighDim HeavyTail") {
    y <- rbinom(n, 1, 0.5)
    for (j in 1:p_relevant) {
      X[y == 0, j] <- rnorm(sum(y == 0), 0, 1)
      X[y == 1, j] <- rt(sum(y == 1), df = 2) + 0.5
    }
  } else if (scenario == "Case10: Mirrored Skew") {
    y <- rbinom(n, 1, 0.5)
    w <- 0.3; sg <- 1.5; dl <- 4
    ms <- .tpsc_meanshift(w, sg, dl)              # E[X] of the mode-0 TPSC draw
    s1 <- sg * sqrt(w / (1 - w)); s2 <- sg * sqrt((1 - w) / w)
    rtpsc0 <- function(m) {                       # mode-0 TPSC, centred to mean 0
      left <- runif(m) < w; z <- abs(rt(m, dl))
      ifelse(left, -z * s1, z * s2) - ms
    }
    for (j in 1:p_relevant) {
      X[y == 0, j] <- -rtpsc0(sum(y == 0))        # left-skew mirror, mean 0
      X[y == 1, j] <-  rtpsc0(sum(y == 1))        # right-skew, mean 0
    }
  } else if (scenario == "Case11: Light Bimodal" ||
             scenario == "Case12: Heavy Bimodal") {
    heavy <- (scenario == "Case12: Heavy Bimodal")
    y <- rbinom(n, 1, 0.5); n0 <- sum(y == 0); n1 <- sum(y == 1)
    rho <- 0.6; df <- if (heavy) 6 else 200
    slo <- if (heavy) 0.48 else 0.50; shi <- if (heavy) 0.58 else 0.65
    m1lo <- if (heavy) -1.6 else -1.5; m1hi <- if (heavy) -1.2 else -1.1
    cs  <- .ar1_corr(p_relevant, rho)
    pw  <- runif(p_relevant, 0.58, 0.66)         # weight of the (heavier, left) mode
    m1v <- runif(p_relevant, m1lo, m1hi)         # heavy-mode centre (< 0)
    sv  <- runif(p_relevant, slo, shi)           # component sd
    m2v <- -pw*m1v/(1 - pw)                        # light-mode centre so that mean = 0
    gen_cls <- function(nc, mir) {                # one class's informative block
      Gm <- MASS::mvrnorm(nc, rep(0, p_relevant), cs)   # correlated mode selection
      Ez <- MASS::mvrnorm(nc, rep(0, p_relevant), cs)   # correlated component noise
      M  <- matrix(0, nc, p_relevant)
      for (j in 1:p_relevant) {
        left <- pnorm(Gm[, j]) < pw[j]
        z    <- if (df >= 100) Ez[, j] else qt(pnorm(Ez[, j]), df)   # normal or t_df
        xj   <- ifelse(left, m1v[j], m2v[j]) + sv[j]*z
        M[, j] <- if (mir) -xj else xj
      }
      M
    }
    X[y == 0, 1:p_relevant] <- gen_cls(n0, FALSE)   # right-skew bimodal
    X[y == 1, 1:p_relevant] <- gen_cls(n1, TRUE)    # left-skew mirror
  } else stop("unknown scenario: ", scenario)

  list(X = X, y = as.factor(y), n_signal = n_signal)
}

ALL_SCENARIOS <- c(
  "Case1: Normal", "Case2: Undefined Moments", "Case3: Skewed LogNormal",
  "Case4: Normal vs. HeavyTail", "Case5: Outliers", "Case6: Unequal Covariance",
  "Case7: Signal Correlation", "Case8: Nonlinear", "Case9: HighDim HeavyTail",
  "Case10: Mirrored Skew", "Case11: Light Bimodal",
  "Case12: Heavy Bimodal")

# The mean-matched shape-signal Cases: class 1 is the mirror of class 0, so the
# class means are equal and only higher-order structure discriminates. Cases
# 11/12 are the multimodal members (light- and heavy-tailed components).
MEANMATCHED_SCENARIOS <- c("Case10: Mirrored Skew",
  "Case11: Light Bimodal", "Case12: Heavy Bimodal")
MULTIMODAL_SCENARIOS  <- c("Case11: Light Bimodal",
  "Case12: Heavy Bimodal")
