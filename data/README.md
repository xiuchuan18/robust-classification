# `data/`

This directory is intentionally empty in the repository. None of the three
datasets used in the manuscript is redistributed here:

| dataset | Section | how it is obtained |
|---|---|---|
| Prostate cancer (`p = 6033`) | 5.1 | `data(prostate, package = "spls")` |
| Pima Indians diabetes | 5.2 | `data("PimaIndiansDiabetes", package = "mlbench")` |
| HTRU2 pulsar | 5.3 | downloaded once from the UCI repository by `cache_htru2()` in `R/04_datasets.R`, which writes `data/HTRU_2.csv` |

`data/HTRU_2.csv` is git-ignored, so the first run of `run_real_data.Rmd`
fetches it and every later run reads it from disk. If the UCI mirror is
unavailable, place a copy of `HTRU_2.csv` (the raw comma-separated file, no
header, nine columns) in this directory and everything downstream works
unchanged.

HTRU2 is distributed by the UCI Machine Learning Repository under CC BY 4.0;
please cite Lyon et al. (2016) if you use it.
