# 04_datasets.R
# Real-data loaders return raw data. All supervised preprocessing (mode-diff screening, median imputation) is fit on the Training
# split only and then applied to the test split, matching the protocols of Fan & Lv (2008, screening on the n=training sample) and Xiong et al. (2025,
# "apply the MD filter to the training set").

## mode-diff (MD) filter, fit on training data (Xiong et al., 2025)
.kde_mode2 <- function(x) { d <- density(x, na.rm = TRUE); d$x[which.max(d$y)] }

# rank features by |mode(X_j | class1) - mode(X_j | class2)| using Train only.
md_scores <- function(X_train, y_train) {
  cl <- unique(y_train)
  s <- vapply(seq_len(ncol(X_train)), function(j)
    abs(.kde_mode2(X_train[y_train == cl[1], j]) - .kde_mode2(X_train[y_train == cl[2], j])), numeric(1))
  names(s) <- colnames(X_train)
  sort(s, decreasing = TRUE)
}

# Fit the screen on train and return the chosen column names (apply to both splits).
screen_md_fit <- function(X_train, y_train, keep = 200) {
  names(md_scores(X_train, y_train))[seq_len(min(keep, ncol(X_train)))]
}

## median imputation of 0-coded missings, fit on training data
# Returns a function that imputes a matrix using Train medians.
impute_zero_median_fit <- function(X_train, cols) {
  med <- sapply(cols, function(c) {
    v <- X_train[, c]; v[v == 0] <- NA; median(v, na.rm = TRUE)
  })
  function(X) {
    X <- as.data.frame(X)
    for (c in cols) { v <- X[[c]]; v[v == 0] <- med[[c]]; X[[c]] <- v }
    as.matrix(X)
  }
}

## raw data loaders
load_prostate_raw <- function() {
  data(prostate, package = "spls")
  X <- prostate$x; colnames(X) <- paste0("gene_", seq_len(ncol(X)))
  list(X = X, y = as.factor(prostate$y), screen = 200, impute = NULL)  # p=6033
}

load_pima_raw <- function() {
  data("PimaIndiansDiabetes", package = "mlbench")
  d <- get("PimaIndiansDiabetes")
  y <- as.factor(ifelse(d$diabetes == "pos", 1, 0))
  X <- as.matrix(d[, setdiff(colnames(d), "diabetes")])
  list(X = X, y = y, screen = NULL,
       impute = c("glucose", "pressure", "triceps", "insulin", "mass"))
}

load_htru2_raw <- function(path = here::here("data", "HTRU_2.csv")) {
  if (!is.null(path) && file.exists(path)) {
    d <- read.csv(path, header = FALSE)
  } else {
    tf <- tempfile(fileext = ".zip")
    utils::download.file(
      "https://archive.ics.uci.edu/static/public/372/htru2.zip", tf,
      quiet = TRUE, mode = "wb")
    d <- read.csv(unz(tf, "HTRU_2.csv"), header = FALSE); unlink(tf)
  }
  X <- as.matrix(d[, 1:8]); colnames(X) <- paste0("V", 1:8)
  list(X = X, y = as.factor(d[, 9]), screen = NULL, impute = NULL)
}

# Download HTRU2 once to a local CSV so every resample reads from disk instead of hitting the UCI mirror.
cache_htru2 <- function(path = here::here("data", "HTRU_2.csv")) {
  if (file.exists(path)) return(invisible(path))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tf <- tempfile(fileext = ".zip")
  utils::download.file(
    "https://archive.ics.uci.edu/static/public/372/htru2.zip", tf,
    quiet = TRUE, mode = "wb")
  d <- read.csv(unz(tf, "HTRU_2.csv"), header = FALSE); unlink(tf)
  write.table(d, path, sep = ",", row.names = FALSE, col.names = FALSE)
  invisible(path)
}

# Registry used by the real-data runner. Retinopathy is code-only (not in paper).
DATASET_LOADERS <- list(
  Prostate  = load_prostate_raw,
  Pima      = load_pima_raw,
  HTRU2     = load_htru2_raw
)

## assemble one train/test split with train-only preprocessing
make_split <- function(spec, train_idx) {
  X <- as.matrix(spec$X); y <- spec$y
  Xtr <- X[train_idx, , drop = FALSE]; ytr <- y[train_idx]
  Xte <- X[-train_idx, , drop = FALSE]; yte <- y[-train_idx]

  if (!is.null(spec$impute)) {
    imp <- impute_zero_median_fit(Xtr, spec$impute)
    Xtr <- imp(Xtr); Xte <- imp(Xte)
  }
  if (!is.null(spec$screen)) {
    keep <- screen_md_fit(Xtr, ytr, keep = spec$screen)
    Xtr <- Xtr[, keep, drop = FALSE]; Xte <- Xte[, keep, drop = FALSE]
  }
  list(X_train = Xtr, y_train = ytr, X_test = Xte, y_test = yte)
}
