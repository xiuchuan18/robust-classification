# Robust Classification with Feature Selection Using Non-Gaussian Data

Replication package for the manuscript

> Xiuchuan Liu and Xianzheng Huang. *Robust Classification with Feature Selection Using Non-Gaussian Data.*

This repository contains all code for modelling and numerical experiments in the main article and its supplement.

Everything is R, with one Rcpp translation unit (`src/tpsc.cpp`) implementing the two-piece scale-Cauchy (TPSC) density, its likelihood, and the
log-likelihood-ratio transformation. **Run every script with this folder as the working directory** — open `robust-classification.Rproj` in RStudio, or `cd`
here before calling `Rscript`. Each entry point locates the project root as the directory containing both `src/tpsc.cpp` and `R/`.

---

## 1. Paper to code map

The main article carries Figures 1-3 and Table 1. Everything else is in the supplement, numbered by appendix (`E.1`, `F.2`, `G.4`, ...).

### Main article

| Manuscript | What it is | Entry point | Result |
|---|---|---|---|
| **Section 4, Figure 1** | AUC boxplots, eight classifiers over twelve simulation Cases | `run_simulation.Rmd`, then `visualize_results.R` | `results/simulation_results.rds`, `docs/run_simulation_visualization.pdf` (page 1) |
| **Section 4, Table 1** | Size of the selected set, TPR, FPR and selection F1 for S-RoLLR vs. L1-RL, averaged over 300 runs | `run_simulation.Rmd` | `results/simulation_results.rds` |
| **Section 5, Figure 2** | AUC boxplots on prostate, Pima and HTRU2 | `run_real_data.Rmd` | written on run: `results/real_data_results.rds` |
| **Section 5.2, Figure 3** | Pima diagnostics: BMI densities, the BMI transformation, two waterfall plots | `run_real_data.Rmd` | written on run |
| **Section 6** | Narrative only. The experiments behind it are reported in Appendix G | see the Appendix G rows below | |

### Supplement, Appendix E: Additional Simulation Results

| Manuscript | What it is | Entry point | Result |
|---|---|---|---|
| **Table E.1** | The twelve simulation Cases (generative models of the informative features) | implemented in `R/03_data_gen.R` | -- |
| **Table E.2** | Median computing time across 300 replicates, twelve Cases | `run_simulation.Rmd` | `results/simulation_results.rds` |
| **Figure E.1** | Selection frequencies of the top-15 genes, S-RoLLR vs. L1-RL (prostate) | `run_variable_selection.Rmd` | `results/variable_selection_freq.rds` |
| **Figure E.2** | Why gene 2215 is picked only by S-RoLLR | `gene2215_analysis.Rmd` | written on run: `figure_gene2215.pdf` |

### Supplement, Appendix F: Numerical implementation details

Stage 1 only: how the TPSC parameter vector is estimated and how the estimator behaves. All four experiments run from `run_all.R` and are reported by
`run_tpsc_estimation.Rmd`.

| Manuscript | What it is | Lab code | Result |
|---|---|---|---|
| **Table F.1** | The four TPSC configurations S1-S4 used throughout Appendix F | `E_CONFIGS` in `R/07_tpsc_lab.R` | -- |
| **Table F.2** | Three initialisation/optimisation strategies over 500 replicates | E1a, `R/08_exp_E1.R` | `results/E1a_summary.csv`, `results/E1a_raw.rds` |
| **Table F.3** | Sensitivity of the fitted distribution to the parameter box | E1b, `R/08_exp_E1.R` | `results/E1b_summary.csv`, `results/E1b_raw.rds` |
| **Table F.4** | Empirical convergence rates from log-log regressions | E2a, `R/09_exp_E2.R` | `results/E2a_parameter_rate.csv`, `results/E2a_kl_rate.csv`, `results/E2a_summary.csv` |
| **Table F.5** | The delta-oracle comparison at n = 1000 | E2b, `R/09_exp_E2.R` | `results/E2b_oracle.csv`, `results/E2b_summary.csv` |
| **Table F.6** | Number of cross-validation folds used to select lambda | fold rule in `R/06_fit_all.R` | -- |

The manuscript renames the three strategies of Table F.2, while the code keeps the lab ids:

| code | manuscript |
|---|---|
| `M0: 4 deterministic starts (shipped)` | **Ours** |
| `M1: 50 LHS random starts` | **Greedy** |
| `M2: Nelder-Mead` | **Simplex** |

### Supplement, Appendix G: The S-RoLLR(D) family

The working-distribution study summarized in Section 6 of the main article.

| Manuscript | What it is | Entry point | Result |
|---|---|---|---|
| **Table G.1** | Generative models of the extension experiments | `R/11_ext_mech_continuous.R` (C1-C6), `R/14_ext_experiments.R` (D1-D5, E1-E5) | -- |
| **Table G.2** | The distribution-estimation procedures | `R/10_ext_densities.R` (continuous), `R/12_ext_pmf_discrete.R` (discrete) | -- |
| **Table G.3** | Transformation sampling variability, three competitive continuous variants | `run_extension.R` with `EXP=exp1`, then `analyze_extension.R` | `results/Sec6_continuous_raw.rds` |
| **Table G.4** | Continuous settings C1-C6 over seven working densities | `run_extension.R` with `EXP=exp1`, then `analyze_extension.R` | `results/Sec6_continuous_raw.rds` |
| **Table G.5** | Discrete settings D1-D5 over five working pmfs | `run_extension.R` with `EXP=exp2`, then `analyze_extension.R` | `results/Sec6_discrete_raw.rds` |
| **Table G.6** | Mixed settings E1-E5 over four variants | `run_extension.R` with `EXP=exp3`, then `analyze_extension.R` | `results/Sec6_mixed_raw.rds` |
| **Appendix G.2**, the sample-size paragraph | Sweep on C3 behind the 0.061 to 5.6e-5 statement | `run_extension_smalln.R`, then `analyze_extension_smalln.R` | `results/Sec6_smalln_raw.rds` |
| **Table G.4**, the "Median time" row | Reproducibility audit of the timings | `audit_extension_timing.R` | `results/logs/Sec6_timing_audit.txt` |

## 2. Source tree

```
R/
  00_setup.R                packages + Rcpp compile + modules 01-06   (Sections 4-5)
  00b_setup_lab.R           boot for the stage-1 estimation lab       (Appendix F)
  00c_setup_ext.R           boot for the extension study, ext_boot()  (Appendix G)
  01_core_tpsc.R            TPSC MLE, the LLR transform, RoLLR / S-RoLLR
  02_competitors.R          MOKE, GMM-NB, L1-RL, SVM, RF, LightGBM
  03_data_gen.R             the twelve simulation Cases of Section 4
  04_datasets.R             real-data loaders and train-only screening/imputation
  05_metrics.R              classification, selection and timing metrics
  06_fit_all.R              one train/test split across all eight classifiers
  07_tpsc_lab.R             shared machinery for the estimation experiments
  08_exp_E1.R               E1: estimators M0-M2 and the parameter box
  09_exp_E2.R               E2: finite-sample behaviour of the TPSC MLE
  10_ext_densities.R        continuous working densities + the plug-in S-RoLLR API
  11_ext_mech_continuous.R  continuous generative settings C1-C6
  12_ext_pmf_discrete.R     discrete working pmfs + the cardinality-adaptive rule
  13_ext_harness.R          Appendix G harness: seed contract, one replicate, T-grid
  14_ext_experiments.R      setting/arm registry for exp1-exp3 + paper labels
src/
  tpsc.cpp                  Rcpp core shared by everything above
data/                       empty by design; see data/README.md
results/                    Result objects (see the tables in Section 1)
results/logs/               console transcripts of the archived Appendix G runs
docs/                       rendered outputs (simulation notebook, figure PDF, session info)
```

The three setup files are deliberately separate. `00_setup.R` pulls in every
competing classifier and is what Sections 4-5 need. `00b_setup_lab.R` needs only
Rcpp and base R, because the estimation lab concerns stage 1 alone.
`00c_setup_ext.R` adds `sn` and the five extension modules, and is the single
loader that both the master session and every parallel worker call, so the two
environments cannot drift apart.

## 3. Reproducing Appendix G (the extension study)

```bash
Rscript -e "source('R/00c_setup_ext.R'); ext_boot()"
```

Run the three experiments (each is independent; `EXP` selects which):

```bash
EXP=exp1 Rscript run_extension.R
```

```bash
EXP=exp2 Rscript run_extension.R
```

```bash
EXP=exp3 Rscript run_extension.R
```

```bash
Rscript run_extension_smalln.R
```

On PowerShell, set the variable first (`$env:EXP="exp1"`) and then call `Rscript run_extension.R`. Defaults are the manuscript's design: n = 800,
p = 100, 20 informative features, 70% training split, 300 Monte Carlo replicates. Everything is overridable through environment variables (`N_REPS`,
`NN`, `PP`, `PREL`, `NCORE`, `OUT`, and `NGRID`/`TEST_N`/`MECHS` for the sweep), so a two-minute smoke test is

```bash
N_REPS=2 PP=20 PREL=6 NN=200 NCORE=2 OUT=smoke.rds EXP=exp1 Rscript run_extension.R
```

Then produce the tables:

```bash
Rscript analyze_extension.R
```

```bash
Rscript analyze_extension_smalln.R
```

```bash
Rscript audit_extension_timing.R
```

`analyze_extension.R` prints Tables G.4, G.5, G.6 and G.3 in that order, followed by the directional test behind the Appendix G.3 claim that the adaptive rule's
advantage grows with the nominal fraction.

Cost of the archived runs on 14 cores: exp1 about 5.7 h, exp2 about 32 min, exp3 about 36 min, sweep about 1 h. The analysis scripts read the stored `.rds`
files and finish in seconds.

### Reproducibility contract

One job is one (setting, replicate); its seed is `SEED_BASE + 1e6 * setting_index + replicate`, so seeds are unique, far apart,
and independent of how `foreach` schedules jobs across workers. Re-running a single job in isolation reproduces it. The seed is consumed by the data
generation and the train/test split only; every arm then sees byte-identical training data, an identical split, and an identical vector of cross-validation
fold labels. The comparison is therefore fully paired, which is what makes the small per-replicate differences in Appendix G estimable. No stage-1 fitter
consumes randomness, so the arm order cannot affect any result. The T-evaluation grids are precomputed once in the master under `GRID_SEED` and
exported to the workers, so no worker ever calls `set.seed()` itself. Details in `R/13_ext_harness.R`.

### Setting labels

The continuous setting ids in the code (`C1: Normal` through `C6: Contam`) are the manuscript's own. The discrete and mixed ids were fixed before the appendix
was written, and are left unchanged so that the archived objects in `results/` stay byte-consistent with the run logs in `results/logs/`. `PAPER_MECH` in
`R/14_ext_experiments.R` maps them one-to-one onto the manuscript's labels, and every analysis script prints the manuscript's labels:

| code id | manuscript | code id | manuscript |
|---|---|---|---|
| `E1: Poisson` | D1 | `X1: nom 0` | E1 (rho = 0) |
| `E2: NegBin` | D2 | `X2: nom 25` | E2 (rho = 0.25) |
| `E3: OrdinalLowCard` | D3 | `X3: nom 50` | E3 (rho = 0.5) |
| `E4: NominalNM` | D4 | `X4: nom 75` | E4 (rho = 0.75) |
| `E5: ContamCount` | D5 | `X5: nom 100` | E5 (rho = 1) |

Arm ids map to the manuscript's variant names through `ARM_LABEL` in `R/13_ext_harness.R`: `identity` to L1-RL, `snorm` to skew-normal, `mixture` to
mixture (BIC), `empirical` to empirical pmf, `auto` to adaptive, and so on.

> **Note on the archived logs.** The transcripts in `results/logs/` were
> > captured before the manuscript's tables were renumbered into appendices, so
> their headers read `TABLE 4/5/6/7` where the paper now reads `G.3/G.4/G.5/G.6`.
> The numbers in the bodies are unchanged. Re-running `analyze_extension.R`
> prints the current numbering.

## 4. Reproducing Sections 4-5 and Appendices E-F

```bash
Rscript -e "rmarkdown::render('run_simulation.Rmd')"
```

`run_simulation.Rmd` reuses `results/simulation_results.rds` if it is present;
set `FORCE_RERUN <- TRUE` in the chunk to recompute. Then

```bash
Rscript visualize_results.R
```

writes `docs/run_simulation_visualization.pdf`: 7 pages, in order AUC, accuracy,
classification F1, selection TPR, selection FPR, selection F1, computing time.
Page 1 is Figure 1 of the main article.

`run_real_data.Rmd` downloads HTRU2 once through `cache_htru2()` and always
recomputes. `run_variable_selection.Rmd` runs 500 resamples of the prostate data
and writes `results/variable_selection_freq.rds`, which `gene2215_analysis.Rmd`
then reads for Figure E.2.

For Appendix F's estimation lab:

```bash
Rscript run_all.R
```

`run_all.R` skips any block whose summary CSV already exists in `results/`;
delete the CSV to recompute. `run_tpsc_estimation.Rmd` renders the resulting
tables.

## 5. Requirements

R >= 4.3 with a working C++ toolchain (Rtools on Windows, Xcode command line tools on macOS, `r-base-dev` on Linux). The setup scripts install anything
missing from CRAN on first use:

`Rcpp`, `glmnet`, `e1071`, `randomForest`, `lightgbm`, `mclust`, `pROC`, `MASS`, `Matrix`, `dplyr`, `tidyr`, `tibble`, `ggplot2`, `gridExtra`, `RColorBrewer`,
`foreach`, `doParallel`, `sn`, `here`, `spls` (prostate data), `mlbench` (Pima data), `rmarkdown`.

The Results were produced with R 4.4.1 on Windows 11 with 14 worker processes; the full session record is in `docs/session_info.txt`.

## 6. Data availability

No dataset is redistributed in this repository. Prostate and Pima ship with the CRAN packages `spls` and `mlbench`; HTRU2 is fetched once from the UCI Machine
Learning Repository into the git-ignored `data/` directory. See `data/README.md`.

## 7. License and citation

Code is released under the MIT License (see `LICENSE`). If you use it, please
cite the manuscript; `CITATION.cff` carries the machine-readable metadata.
