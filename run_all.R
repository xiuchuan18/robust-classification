# run_all.R
# One-shot re-run of the whole TPSC estimation lab: E1 and E2.  Works both ways:
#
#   RStudio Console (from the project root, .Rproj open):
#       source("run_all.R")
#
#   Shell:
#       Rscript run_all.R

t_start <- Sys.time()
.say <- function(...) cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

source("R/00b_setup_lab.R")

OUTDIR <- file.path(ROOT, "results")
dir.create(OUTDIR, showWarnings = FALSE)

.block <- function(name, outfile, expr) {
  fp <- file.path(OUTDIR, outfile)
  if (file.exists(fp)) {
    .say("%s: %s already exists -- skipping (delete it to redo)", name, outfile)
    return(invisible(NULL))
  }
  .say("%s: starting", name)
  t0  <- Sys.time()
  res <- force(expr)
  .say("%s: done in %.1f min", name,
       as.numeric(difftime(Sys.time(), t0, units = "mins")))
  invisible(res)
}

## 1. E1: estimators (M0-M2) + parameter box
E1 <- .block("E1", "E1a_summary.csv", run_E1())

## 2. E2: finite-sample properties of the TPSC MLE
E2 <- .block("E2", "E2a_summary.csv", run_E2())

.say("ALL DONE in %.1f min.  Files in %s:",
     as.numeric(difftime(Sys.time(), t_start, units = "mins")), basename(OUTDIR))
print(list.files(OUTDIR))
cat("\nSend these for analysis:\n",
    paste(" -", c("E1a_summary.csv", "E1b_summary.csv", "E2a_summary.csv",
                  "E2b_summary.csv"),
          collapse = "\n"), "\n")
