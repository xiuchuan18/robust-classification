# `data/`

This directory is intentionally empty in the repository. None of the three
datasets used in the manuscript is redistributed here:

| dataset | Section | how it is obtained |
|---|---|---|
| Prostate cancer (`p = 6033`) | 5.1 | `data(prostate, package = "spls")` |
| Pima Indians diabetes | 5.2 | `data("PimaIndiansDiabetes", package = "mlbench")` |
| HTRU2 pulsar | 5.3 | downloaded once from the UCI repository by `cache_htru2()` in `R/04_datasets.R`, which writes `data/HTRU_2.csv` |
