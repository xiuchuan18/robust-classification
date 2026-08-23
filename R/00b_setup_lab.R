# 00b_setup_lab.R ------------------------------------------------------------
# Minimal setup for the TPSC-estimation experiments E1-E2.  Unlike 00_setup.R
# it needs *only* Rcpp + base/stats: E1-E2 concern stage-1 estimation only, so
# glmnet / e1071 / randomForest / lightgbm / mclust are never touched.
#
# Usage (from the project root, e.g. with the .Rproj open):
#
#     source("R/00b_setup_lab.R")
#
# The project root is the directory that contains BOTH src/tpsc.cpp and R/.

## ---- locate the project root ------------------------------------------------
.is_root <- function(d) file.exists(file.path(d, "src", "tpsc.cpp")) &&
                        dir.exists(file.path(d, "R"))

.find_root <- function() {
  cand <- getwd()
  for (i in 1:4) {                       # walk up in case wd is a subfolder (e.g. R/)
    if (.is_root(cand)) return(cand)
    parent <- dirname(cand)
    if (identical(parent, cand)) break
    cand <- parent
  }
  stop("Cannot find the project root from getwd() = '", getwd(), "'.\n",
       "  The root is the folder containing both src/tpsc.cpp and R/.\n",
       "  In RStudio, open the project's .Rproj file (or use\n",
       "  Session > Set Working Directory > To Project Directory) and retry.",
       call. = FALSE)
}

ROOT <- .find_root()

## ---- compile the C++ core ---------------------------------------------------
# Needs a working toolchain: Rtools on Windows, Xcode command line tools on
# macOS ("xcode-select --install"), r-base-dev on Linux.
suppressMessages(library(Rcpp))
Rcpp::sourceCpp(file.path(ROOT, "src", "tpsc.cpp"))

## ---- source the R modules ---------------------------------------------------
for (f in c("01_core_tpsc.R", "07_tpsc_lab.R",
            "08_exp_E1.R", "09_exp_E2.R")) {
  fp <- file.path(ROOT, "R", f)
  if (!file.exists(fp)) stop("missing module: ", fp, call. = FALSE)
  source(fp)
}
# 01_core_tpsc.R defines fit_srollr(), which calls cv.glmnet() at *run* time
# only; sourcing it without glmnet attached is therefore safe for E1-E2.

dir.create(file.path(ROOT, "results"), showWarnings = FALSE)

message("TPSC estimation lab ready.  root = ", ROOT,
        "\n  E1: run_E1()   E2: run_E2()")
