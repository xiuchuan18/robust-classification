# 14_ext_experiments.R
#
#   exp1  Appendix G.2  continuous settings C1-C6  *  seven working densities
#   exp2  Appendix G.3  discrete settings D1-D5    *  five working pmfs
#   exp3  Appendix G.3  mixed settings E1-E5       *  four variants
#

EXP2_MECHS <- c("E1: Poisson", "E2: NegBin", "E3: OrdinalLowCard",
                "E4: NominalNM", "E5: ContamCount")
EXP2_ORACLE <- c("E1: Poisson" = "poisson", "E2: NegBin" = "nbinom",
                 "E3: OrdinalLowCard" = "empirical",
                 "E4: NominalNM" = "empirical", "E5: ContamCount" = NA)

# Non-monotone nominal weights, symmetric about the midpoint of {0,...,5} and solved so that both the class mean (2.5) and the class variance (35/12) match
# the uniform pmf of class 0.
NOMINAL_Q_SIG  <- c(0.220, 0.007, 0.273, 0.273, 0.007, 0.220)
NOMINAL_Q_NULL <- rep(1/6, 6)

generate_exp2_data <- function(scenario, n, p, p_relevant) {
  y <- rbinom(n, 1, 0.5); n0 <- sum(y == 0); n1 <- n - n0
  X <- matrix(0, n, p)
  spec <- switch(scenario,
    "E1: Poisson" = list(
      g = function(m, inf) rpois(m, 5 * exp(if (inf) 0.30 else 0)),
      type = "poisson"),
    "E2: NegBin" = list(
      g = function(m, inf) rnbinom(m, size = 0.5, mu = 6 * exp(if (inf) 0.55 else 0)),
      type = "nbinom"),
    "E3: OrdinalLowCard" = list(
      g = function(m, inf) pmin(4L, rpois(m, 1.2 + if (inf) 1.0 else 0)),
      type = "empirical"),
    "E4: NominalNM" = list(
      # six unordered levels, mean- and variance-matched: the only mechanism in which the transform is not affine in x and the pmf is indispensable
      g = function(m, inf) sample(0:5, m, replace = TRUE,
                                  prob = if (inf) NOMINAL_Q_SIG else NOMINAL_Q_NULL),
      type = "empirical"),
    "E5: ContamCount" = list(
      # counts with 5% of the observations inflated tenfold: no pmf is correct
      g = function(m, inf) { v <- rpois(m, 4 * exp(if (inf) 0.40 else 0))
        k <- rbinom(m, 1, 0.05); v * (1 - k) + k * v * 10 },
      type = "poisson"),
    stop("unknown Exp-2 scenario: ", scenario))

  for (j in seq_len(p)) {
    X[y == 0, j] <- spec$g(n0, FALSE)
    X[y == 1, j] <- spec$g(n1, j <= p_relevant)
  }
  colnames(X) <- paste0("V", seq_len(p))
  list(X = X, y = as.factor(y), n_signal = p_relevant,
       signal = paste0("V", seq_len(p_relevant)), types = rep(spec$type, p))
}

## exp3
EXP3_MECHS <- c("X1: nom 0", "X2: nom 25", "X3: nom 50", "X4: nom 75",
                "X5: nom 100")
EXP3_NOMFRAC <- c("X1: nom 0" = 0, "X2: nom 25" = 0.25, "X3: nom 50" = 0.50,
                  "X4: nom 75" = 0.75, "X5: nom 100" = 1)

generate_exp3_data <- function(scenario, n, p, p_relevant, pi_disc = 0.5,
                               shift = 0.75) {
  nomf   <- unname(EXP3_NOMFRAC[scenario])
  n_disc <- round(pi_disc * p)
  y <- rbinom(n, 1, 0.5); n0 <- sum(y == 0); n1 <- n - n0
  X <- matrix(0, n, p); types <- character(p)

  k_disc_sig <- round(pi_disc * p_relevant)          # informative discrete cols
  sig_disc   <- seq_len(k_disc_sig)                  #   = 1 .. k_disc_sig
  noise_disc <- if (n_disc > k_disc_sig) (k_disc_sig + 1):n_disc else integer(0)
  sig <- c(sig_disc, n_disc + seq_len(p_relevant - k_disc_sig))


  is_nom <- rep(FALSE, p)
  if (k_disc_sig > 0)      is_nom[head(sig_disc,   round(nomf * k_disc_sig))]   <- TRUE
  if (length(noise_disc))  is_nom[head(noise_disc, round(nomf * length(noise_disc)))] <- TRUE

  for (j in seq_len(p)) {
    inf <- j %in% sig
    if (is_nom[j]) {                                    # unordered nominal
      types[j] <- "empirical"
      X[y == 0, j] <- sample(0:5, n0, TRUE, NOMINAL_Q_NULL)
      X[y == 1, j] <- sample(0:5, n1, TRUE,
                            if (inf) NOMINAL_Q_SIG else NOMINAL_Q_NULL)
    } else if (j <= n_disc) {                           # exponential-family count
      types[j] <- "poisson"
      X[y == 0, j] <- rpois(n0, 4)
      X[y == 1, j] <- rpois(n1, 4 * exp(if (inf) 0.45 else 0))
    } else {                                            # continuous, skew-heavy
      types[j] <- "tpsc"
      X[y == 0, j] <- sn::rst(n0, 0, 1.5, 3, 3)
      X[y == 1, j] <- sn::rst(n1, if (inf) shift else 0, 1.5, 3, 3)
    }
  }
  colnames(X) <- paste0("V", seq_len(p))
  list(X = X, y = as.factor(y), n_signal = length(sig),
       signal = paste0("V", sig), types = types, nominal_frac = nomf,
       n_nom_signal = sum(is_nom[sig]))
}


.fit_llr_binned <- function(x0, x1, K = 8L) {
  br <- unique(as.numeric(quantile(c(x0, x1), seq(0, 1, length.out = K + 1))))
  if (length(br) < 3L) return(.fit_llr_empirical(x0, x1))   # already coarse
  br[1] <- -Inf; br[length(br)] <- Inf
  f <- .fit_llr_empirical(cut(x0, br, labels = FALSE), cut(x1, br, labels = FALSE))
  function(x) f(cut(x, br, labels = FALSE))
}



EXP1_ARMS <- c("identity", "gaussian", "studentt", "snorm", "skewt",
               "tpsc", "mixture")
EXP2_ARMS <- c("identity", "poisson", "nbinom", "empirical", "tpsc")
EXP3_ARMS <- c("identity", "tpsc", "binned", "auto")

use_experiment <- function(which = c("exp1", "exp2", "exp3")) {
  which <- match.arg(which)
  ge <- globalenv()

  .DISC_FITTERS$binned <- .fit_llr_binned
  assign(".DISC_FITTERS", .DISC_FITTERS, envir = ge)

  mechs <- switch(which, exp1 = CONT_SCENARIOS, exp2 = EXP2_MECHS, exp3 = EXP3_MECHS)
  arms  <- switch(which, exp1 = EXP1_ARMS,      exp2 = EXP2_ARMS,  exp3 = EXP3_ARMS)

  assign("ALL_MECHANISMS", mechs, envir = ge)
  assign("ARM_ORDER", arms, envir = ge)
  assign("EXPERIMENT", which, envir = ge)
  assign("mech_block", function(mech) which, envir = ge)
  assign("arm_coverage", function(mech) arms, envir = ge)
  assign("generate_any", switch(which,
    exp1 = function(mech, n, p, k) generate_cont_data(mech, n, p, k),
    exp2 = function(mech, n, p, k) generate_exp2_data(mech, n, p, k),
    exp3 = function(mech, n, p, k) generate_exp3_data(mech, n, p, k)),
    envir = ge)
  assign("MECH_ORACLE", switch(which,
    exp1 = c(CONT_ORACLE, "C1: Normal" = "identity"),   # equal variances: the
    exp2 = EXP2_ORACLE,                                  # homoscedastic Gaussian
    exp3 = setNames(rep("auto", length(EXP3_MECHS)), EXP3_MECHS)),  # true oracle
    envir = ge)
  invisible(mechs)
}

## Manuscript labels
# Code setting id -> the label printed in the manuscript.
PAPER_MECH <- c(
  # exp1: already the manuscript's labels
  "C1: Normal" = "C1", "C2: SkewNormal" = "C2", "C3: SkewT" = "C3",
  "C4: MirrorTPSC" = "C4", "C5: Bimodal" = "C5", "C6: Contam" = "C6",
  # exp2: manuscript Table G.5
  "E1: Poisson" = "D1", "E2: NegBin" = "D2", "E3: OrdinalLowCard" = "D3",
  "E4: NominalNM" = "D4", "E5: ContamCount" = "D5",
  # exp3: manuscript Table G.6 (rho = nominal fraction of the discrete block)
  "X1: nom 0" = "E1 (rho=0)", "X2: nom 25" = "E2 (rho=0.25)",
  "X3: nom 50" = "E3 (rho=0.5)", "X4: nom 75" = "E4 (rho=0.75)",
  "X5: nom 100" = "E5 (rho=1)")

# Vector of code ids -> vector of manuscript labels (unknown ids pass through).
paper_mech <- function(x) {
  x <- as.character(x)
  lab <- unname(PAPER_MECH[x])
  ifelse(is.na(lab), x, lab)
}

paper_arm <- function(x) unname(ARM_LABEL[as.character(x)])
