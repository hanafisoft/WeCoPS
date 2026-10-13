# WeCoPS (R)

**Weighting CPCA for Common and Proper Structures.** Does a block weighting lead to the same
interpretation as the default weighting in consensus PCA?

## Installation

```r
# install.packages("remotes")
remotes::install_github("hanafisoft/WeCoPS", subdir = "R/WeCoPS")
```

The package depends only on base R (graphics, grDevices, stats, utils).

## Quick start

```r
library(WeCoPS)
d <- load_butters()                     # VIS, NIR and MIR spectra of 21 fat products
res <- compare_weighting(d$blocks, R = 4, weights = "nvar", block_names = d$block_names)
res                                     # verdict, motifs, presence, reweighting threshold
plot(res)                               # diagnostic sheet

a <- wecops(d$blocks, R = 4, block_names = d$block_names)   # one analysis, two bases
a$motif_pc; a$motif                     # principal / common and proper components
a$mixing$mu_star                        # mixing of each principal component
a$factorization$ccp[[1]]$D              # D_1 of the canonical factorisation
```

Your own blocks: `prep_blocks(list(X1, X2, ...))` (column centring, unit total inertia).
Weightings: `"rv"`, `"mfa"`, `"statis"`, `"nvar"`, or a vector of K positive weights.

Reproduce the article: `source(system.file("scripts", "reproduce_article.R", package = "WeCoPS"))`.

Author and maintainer: Mohamed Hanafi <mohammed44hanafi@gmail.com>. Licence: MIT.
