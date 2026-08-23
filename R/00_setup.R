# 00_setup.R

.pkgs <- c(
  "Rcpp",            # C++ core
  "glmnet",          # L1 logistic (L1-RL, S-RoLLR second stage)
  "e1071",           # SVM
  "randomForest",    # RF
  "lightgbm",        # LightGBM
  "mclust",          # GMM-NB
  "pROC",            # AUC
  "MASS", "Matrix",  # data generation
  "dplyr", "tidyr",  # summaries
  "ggplot2", "gridExtra", "RColorBrewer",
  "foreach", "doParallel",
  "here",            # project-root resolution used by the .Rmd entry points
  "spls",            # prostate cancer data
  "mlbench"          # Pima Indians diabetes data
)
.missing <- setdiff(.pkgs, rownames(installed.packages()))
if (length(.missing)) {
  message("Installing missing packages: ", paste(.missing, collapse = ", "))
  install.packages(.missing, repos = "https://cloud.r-project.org")
}
invisible(lapply(.pkgs, function(p) suppressMessages(library(p, character.only = TRUE))))

## compile the C++ core
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA)
if (is.na(.here) || !nzchar(.here)) .here <- getwd() else .here <- dirname(.here)
.cpp <- file.path(.here, "src", "tpsc.cpp")
if (!file.exists(.cpp)) .cpp <- file.path(getwd(), "src", "tpsc.cpp")
stopifnot(file.exists(.cpp))
Rcpp::sourceCpp(.cpp)

## source the R modules
.rdir <- dirname(.cpp); .rdir <- file.path(dirname(.rdir), "R")
for (f in c("01_core_tpsc.R", "02_competitors.R", "03_data_gen.R",
            "04_datasets.R", "05_metrics.R", "06_fit_all.R")) {
  fp <- file.path(.rdir, f); if (file.exists(fp)) source(fp)
}
message("RoLLR core ready (Rcpp compiled, modules loaded).")
