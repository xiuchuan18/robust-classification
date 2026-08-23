# 00c_setup_ext.R -------------------------------------------------------------
# Single entry point for the Section 6 extension study (the working-density /
# plug-in experiments of Section 6, reported in Appendix G, Tables G.3-G.6).
#
# Loads everything the study needs, in the right order, from an EXPLICIT project
# directory.  Both the master session and every PSOCK worker call ext_boot(dir),
# so the two environments are guaranteed identical.
#
# Usage (from the project root, e.g. with the .Rproj open):
#
#     source("R/00c_setup_ext.R"); ext_boot()
#     source("R/14_ext_experiments.R")   # loaded by ext_boot() already
#     use_experiment("exp1")             # or "exp2" / "exp3"
#
# The project root is the directory that contains BOTH src/tpsc.cpp and R/.
# Relation to the other two setup files:
#   R/00_setup.R      Sections 4-5   (all eight competing classifiers)
#   R/00b_setup_lab.R Supplement E1-E2 (stage-1 estimation lab, Rcpp only)
#   R/00c_setup_ext.R Section 6      (this file; stage-1 plug-in families)

## ---- extra packages used only by Section 6 ----------------------------------
# sn      : skew-normal / skew-t MLE (Azzalini), the `snorm` and `skewt` arms
# tibble  : column_to_rownames() in the analysis scripts
.ext_pkgs <- c("sn", "tibble")
.ext_missing <- setdiff(.ext_pkgs, rownames(installed.packages()))
if (length(.ext_missing)) {
  message("Installing missing packages: ", paste(.ext_missing, collapse = ", "))
  install.packages(.ext_missing, repos = "https://cloud.r-project.org")
}

ext_boot <- function(dir = NULL) {
  if (is.null(dir)) dir <- .ext_find_root()
  dir <- normalizePath(dir, winslash = "/", mustWork = TRUE)

  # 00_setup.R compiles src/tpsc.cpp and sources R/01-R/06; the extension needs
  # 01_core_tpsc.R (fit_tpsc_mle, d_tpsc) and 05_metrics.R (classification and
  # selection metrics) from it, and reuses the identical Rcpp core.
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  suppressMessages(source(file.path(dir, "R", "00_setup.R")))
  suppressMessages(library(sn))

  for (f in c("10_ext_densities.R", "11_ext_mech_continuous.R",
              "12_ext_pmf_discrete.R", "13_ext_harness.R",
              "14_ext_experiments.R")) {
    fp <- file.path(dir, "R", f)
    if (!file.exists(fp)) stop("missing module: ", fp, call. = FALSE)
    source(fp)
  }

  assign("ROOT", dir, envir = globalenv())
  dir.create(file.path(dir, "results"), showWarnings = FALSE)
  invisible(TRUE)
}

## ---- locate the project root ------------------------------------------------
.ext_is_root <- function(d) file.exists(file.path(d, "src", "tpsc.cpp")) &&
                            dir.exists(file.path(d, "R"))

.ext_find_root <- function() {
  cand <- Sys.getenv("EXT_WD", getwd())
  for (i in 1:4) {                     # walk up in case wd is a subfolder (R/)
    if (.ext_is_root(cand)) return(cand)
    parent <- dirname(cand)
    if (identical(parent, cand)) break
    cand <- parent
  }
  stop("Cannot find the project root from '", Sys.getenv("EXT_WD", getwd()),
       "'.\n  The root is the folder containing both src/tpsc.cpp and R/.\n",
       "  In RStudio, open the project's .Rproj file (or use Session > Set\n",
       "  Working Directory > To Project Directory) and retry.", call. = FALSE)
}

# Where Section 6 writes and reads its result objects.
ext_res <- function(f) {
  root <- if (exists("ROOT", inherits = TRUE)) get("ROOT", inherits = TRUE) else getwd()
  dir.create(file.path(root, "results"), showWarnings = FALSE)
  file.path(root, "results", f)
}
