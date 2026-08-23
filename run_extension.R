# run_extension.R
# Section 6 ("Extension of S-RoLLR"): runs one of the three extension
# experiments and saves its raw per-(setting, arm, replicate) table.
#
#   EXP = "exp1"   continuous settings C1-C6, seven working densities
#                  -> manuscript Tables G.3, G.4     -> results/Sec6_continuous_raw.rds
#   EXP = "exp2"   discrete settings D1-D5, five working pmfs
#                  -> manuscript Table G.5           -> results/Sec6_discrete_raw.rds
#   EXP = "exp3"   mixed settings E1-E5, four variants
#                  -> manuscript Table G.6           -> results/Sec6_mixed_raw.rds
#
# Manuscript design: n = 800, p = 100, |I| = 20, 70% training split, 300 Monte
# Carlo replicates (Appendix G.1).  Those are the defaults below.
#
# Usage (from the project root):
#
#   Shell (PowerShell):
#       $env:EXP="exp1"; Rscript run_extension.R
#       foreach ($e in "exp1","exp2","exp3") { $env:EXP=$e; Rscript run_extension.R }
#
#   Shell (bash):
#       EXP=exp1 Rscript run_extension.R
#
#   RStudio console:
#       Sys.setenv(EXP = "exp1"); source("run_extension.R")

PROJ <- normalizePath(Sys.getenv("EXT_WD", getwd()), winslash = "/")
source(file.path(PROJ, "R", "00c_setup_ext.R")); ext_boot(PROJ)

EXP    <- Sys.getenv("EXP", "exp1")
N_REPS <- as.integer(Sys.getenv("N_REPS", "300"))
NN     <- as.integer(Sys.getenv("NN",   "800"))
PP     <- as.integer(Sys.getenv("PP",   "100"))
PREL   <- as.integer(Sys.getenv("PREL", "20"))
NCORE  <- as.integer(Sys.getenv("NCORE", as.character(max(1, detectCores() - 2))))

DEFAULT_OUT <- c(exp1 = "Sec6_continuous_raw.rds", exp2 = "Sec6_discrete_raw.rds",
                 exp3 = "Sec6_mixed_raw.rds")
OUT <- Sys.getenv("OUT", unname(DEFAULT_OUT[EXP]))
OUT <- if (basename(OUT) == OUT) ext_res(OUT) else OUT

use_experiment(EXP)
TGRIDS <- build_tgrids(PP, PREL)

jobs <- expand.grid(mech_index = seq_along(ALL_MECHANISMS), rep = seq_len(N_REPS))
jobs$mechanism <- ALL_MECHANISMS[jobs$mech_index]
# distinct seed space per experiment so the three never share a data stream
SEED_BASE <- SEED_BASE + 100000000L * match(EXP, c("exp1", "exp2", "exp3"))

cat(sprintf("%s: %d settings x %d reps = %d jobs | arms: %s\n", EXP,
            length(ALL_MECHANISMS), N_REPS, nrow(jobs),
            paste(ARM_ORDER, collapse = ", ")))

cl <- makeCluster(NCORE); registerDoParallel(cl)
clusterExport(cl, c("PROJ", "NN", "PP", "PREL", "TGRIDS", "jobs", "EXP",
                    "SEED_BASE"), envir = environment())
cpp_exports <- c("d_tpsc_cpp", "tpsc_nll_cpp", "tpsc_transform_cpp",
                 "tpsc_rowloglik_cpp")

t0 <- proc.time()[["elapsed"]]
results <- foreach(k = seq_len(nrow(jobs)), .combine = "rbind",
                   .packages = c("Rcpp", "glmnet", "pROC", "sn", "MASS"),
                   .noexport = cpp_exports) %dopar% {
  if (!exists("EXT_READY", envir = globalenv())) {
    source(file.path(PROJ, "R", "00c_setup_ext.R")); ext_boot(PROJ)
    use_experiment(EXP)
    assign("SEED_BASE", SEED_BASE, envir = globalenv())
    assign("EXT_READY", TRUE, envir = globalenv())
  }
  j <- jobs[k, ]
  tryCatch(run_plugin_job(j$mechanism, j$rep, NN, PP, PREL,
                          TGRIDS[[j$mechanism]], j$mech_index),
           error = function(e) { message(sprintf("FAILED %s rep %d: %s",
                                  j$mechanism, j$rep, conditionMessage(e))); NULL })
}
stopCluster(cl); registerDoSEQ()

expected <- length(ALL_MECHANISMS) * length(ARM_ORDER) * N_REPS
cat(sprintf("elapsed %.1f min | rows %d / %d | arm failures %d\n",
            (proc.time()[["elapsed"]] - t0) / 60, nrow(results), expected,
            sum(!results$ok)))
if (any(!results$ok)) print(head(unique(results[!results$ok,
                                  c("mechanism","arm","error")]), 10))
attr(results, "config") <- list(experiment = EXP, n = NN, p = PP,
                                p_relevant = PREL, n_reps = N_REPS)
saveRDS(results, OUT); cat("saved:", OUT, "\n")
