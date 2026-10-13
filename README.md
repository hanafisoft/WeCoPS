# WeCoPS — Weighting CPCA for Common and Proper Structures

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23337145.svg)](https://doi.org/10.5281/zenodo.23337145)

Does a block weighting lead to the same interpretation as the default weighting in consensus
PCA (CPCA)? This repository contains two implementations with the same functions and results:

* `R/WeCoPS/` — R package `WeCoPS` (base R only);
* `python/` — Python package `wecops` (numpy, scipy; matplotlib optional).

For a fixed number `R` of components, both read the retained subspace in two bases (principal
components; common and proper components obtained by concentration, rotation towards a
data-driven target and orientation inside each labelled space), measure how the principal
components mix common and proper components, apply the equivalence rule (same motif, and every
labelled space of the default remains present in the subspace of the other weighting,
`kappa_L >= 0.9`), return the canonical factorisation `X_k ~ B D_k V_k'` on both bases, and draw
the two-panel diagnostic sheet.

## Installation

R:

```r
# install.packages("remotes")
remotes::install_github("hanafisoft/WeCoPS", subdir = "R/WeCoPS")
```

Python:

```bash
pip install "git+https://github.com/hanafisoft/WeCoPS.git#subdirectory=python"
```

## Quick start

```r
library(WeCoPS)
b <- load_butters()$blocks
r <- compare_weighting(b, R = 4, weights = "nvar", block_names = c("VIS", "NIR", "MIR"))
r                 # verdict, motifs, presence of each labelled space, mixing, threshold
plot_fiche(r)     # diagnostic sheet
```

```python
import wecops as wc
b = wc.load_butters()["blocks"]
r = wc.compare_weighting(b, 4, "nvar", block_names=["VIS", "NIR", "MIR"])
print(r.summary())
```

## Reproducing the article

* `R/WeCoPS/inst/scripts/reproduce_article.R` and `python/examples/reproduce_article.py`
  reproduce the 28 verdicts (19 equivalent, 10 false alarms of the principal-component reading).
* `R/WeCoPS/inst/scripts/verify_article.R` checks under R the results of the article against the
  Python reference values (47 checks: mixing table, 28 verdicts, 196 cases of the choice of R,
  sensitivity to the thresholds).
* `R/check_WeCoPS.R` checks the R package against the Python reference (58 checks).

## Data

The packages ship the datasets of the article: wines (sensory profiles, four blocks), butters and
margarines (VIS, NIR and MIR spectra of 21 fat products; Jaillais et al., 2007, distributed with
the agreement of B. Jaillais) and four simulated datasets.

## Author, licence, citation

Author and maintainer: Mohamed Hanafi (StatSC, Oniris-INRAE) <mohammed44hanafi@gmail.com>.
Licence: MIT. Cite as: Hanafi M. (2026). WeCoPS: Weighting CPCA for Common and Proper Structures, version 0.2.1. Zenodo. https://doi.org/10.5281/zenodo.23337145 (see also `CITATION.cff`).
