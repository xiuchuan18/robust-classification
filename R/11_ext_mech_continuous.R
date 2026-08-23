# 11_ext_mech_continuous.R
# Continuous generative mechanisms C1-C6 of Section 6. Setting labels here are exactly the labels used in the manuscript.
#
#   C1 Normal      -> gaussian is the oracle
#   C2 SkewNormal  -> snorm    is the oracle
#   C3 SkewT       -> skewt    is the oracle
#   C4 MirrorTPSC  -> tpsc     is the oracle  (mean-matched: shape-only signal)
#   C5 Bimodal     -> mixture  is the oracle  (mean-matched, unimodal families misspecified)
#   C6 Contam      -> NO family is correct    (10% scale-inflated contamination)
#
# The p - |I| noise features are drawn from the same law as the signal features but with NO class difference, so that every family faces noise it can model
# equally well, otherwise the false-positive rates would not be comparable.

.pl_halft_mean <- function(delta) sqrt(delta/pi) * gamma((delta-1)/2) / gamma(delta/2)

# mode-0 two-piece scale Student-t draw, centred to mean 0
.pl_rtpsc0 <- function(m, w = 0.3, sg = 1.5, dl = 4) {
  s1 <- sg * sqrt(w / (1 - w)); s2 <- sg * sqrt((1 - w) / w)
  ms <- .pl_halft_mean(dl) * ((1 - w) * s2 - w * s1)
  left <- runif(m) < w; z <- abs(rt(m, dl))
  ifelse(left, -z * s1, z * s2) - ms
}

CONT_SCENARIOS <- c("C1: Normal", "C2: SkewNormal", "C3: SkewT", "C4: MirrorTPSC", "C5: Bimodal", "C6: Contam")

CONT_ORACLE <- c("C1: Normal" = "gaussian", "C2: SkewNormal" = "snorm",
                 "C3: SkewT" = "skewt",     "C4: MirrorTPSC" = "tpsc",
                 "C5: Bimodal" = "mixture", "C6: Contam" = NA)

CONT_SHIFT <- c("C1: Normal" = 0.55, "C2: SkewNormal" = 0.40,
                "C3: SkewT" = 0.75,  "C6: Contam" = 0.55)

generate_cont_data <- function(scenario, n, p, p_relevant, shift = NULL) {
  if (is.null(shift)) shift <- unname(CONT_SHIFT[scenario])
  y  <- rbinom(n, 1, 0.5); n0 <- sum(y == 0); n1 <- n - n0
  X  <- matrix(0, n, p)
  ns <- p - p_relevant

  put <- function(f0, f1) {          # f0/f1: functions of the sample size
    for (j in seq_len(p_relevant)) { X[y == 0, j] <<- f0(n0); X[y == 1, j] <<- f1(n1) }
    if (ns > 0) for (j in (p_relevant + 1):p) {   # noise: class-0 law in both classes
      X[y == 0, j] <<- f0(n0); X[y == 1, j] <<- f0(n1)
    }
  }

  if (scenario == "C1: Normal") {
    put(function(m) rnorm(m, 0, 1), function(m) rnorm(m, shift, 1))

  } else if (scenario == "C2: SkewNormal") {
    put(function(m) sn::rsn(m, 0, 1.5, 5), function(m) sn::rsn(m, shift, 1.5, 5))

  } else if (scenario == "C3: SkewT") {
    put(function(m) sn::rst(m, 0, 1.5, 3, 3), function(m) sn::rst(m, shift, 1.5, 3, 3))

  } else if (scenario == "C4: MirrorTPSC") {
    # mean-matched: class 1 is the mirror image of class 0, so the classes share every even-order moment and only the SIGN of the skewness separates them.
    put(function(m) -.pl_rtpsc0(m), function(m) .pl_rtpsc0(m))
    if (ns > 0) for (j in (p_relevant + 1):p) X[, j] <- -.pl_rtpsc0(n)

  } else if (scenario == "C5: Bimodal") {
    # mean-matched two-component mixture, the unimodal families are misspecified.
    pw <- 0.62; m1 <- -1.3; s <- 0.55; m2 <- -pw * m1 / (1 - pw)
    rmix <- function(m, mir) {
      v <- ifelse(runif(m) < pw, m1, m2) + s * rnorm(m)
      if (mir) -v else v
    }
    put(function(m) rmix(m, FALSE), function(m) rmix(m, TRUE))
    if (ns > 0) for (j in (p_relevant + 1):p) X[, j] <- rmix(n, FALSE)

  } else if (scenario == "C6: Contam") {
    # 10% of each class is scale-inflated tenfold.
    rc <- function(m, mu) { k <- round(0.1 * m)
      sample(c(rnorm(m - k, mu, 1), rnorm(k, mu, 10))) }
    put(function(m) rc(m, 0), function(m) rc(m, shift))

  } else stop("unknown continuous scenario: ", scenario)

  colnames(X) <- paste0("V", seq_len(p))
  list(X = X, y = as.factor(y), n_signal = p_relevant,
       signal = paste0("V", seq_len(p_relevant)),
       types = rep("tpsc", p))
}
