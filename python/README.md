# WeCoPS (Python)

**Weighting CPCA for Common and Proper Structures.**

In exploratory multiblock analysis, blocks are usually scaled to unit total inertia, but other
weightings are common (multiple factor analysis, RV normalisation, STATIS compromise, number of
variables). WeCoPS answers the user's question: *I have my results with the default weighting;
would they have changed with another weighting?*

For a fixed number `R` of components, WeCoPS

* reads the retained subspace of the consensus PCA (CPCA) in two bases: the **principal
  components** and the **common and proper components** (concentration, data-driven target at the
  salience threshold `alpha`, rotation towards the target, then orientation inside each labelled
  space by the multiple co-inertia analysis of the blocks of the label, which makes the basis
  unique);
* measures how much each principal component **mixes** common and proper components
  (indices `mu` and `mu_star`);
* applies the **equivalence rule**: a weighting is equivalent to the default if it gives the same
  motif of common and proper components and every labelled space of the default (space spanned by
  the components sharing a label) remains **present** in its subspace (`kappa >= 0.9`);
* returns the four-part decomposition of the inertia of each block, the reweighting threshold
  `lambda_R / lambda_{R+1}`, and the **canonical factorisation** `X_k ~ B D_k V_k'` on both bases;
* draws the two-panel diagnostic sheet.

## Installation

```bash
pip install "git+https://github.com/hanafisoft/WeCoPS.git#subdirectory=python"
pip install matplotlib   # optional, for plot_fiche
```

## Quick start

```python
import wecops as wc

data = wc.load_butters()                       # VIS, NIR and MIR spectra of 21 fat products
blocks, names = data["blocks"], data["block_names"]

res = wc.compare_weighting(blocks, R=4, weights="nvar", block_names=names)
print(res.summary())
wc.plot_fiche(res, file="fiche.png")

one = wc.wecops(blocks, R=4, block_names=names)  # one analysis, two bases
one.motif_pc, one.motif                          # labels of principal / common and proper components
one.mixing["mu_star"]                            # mixing of each principal component
one.factorization["ccp"][0]["D"]                 # D_1 of the canonical factorisation (common and proper basis)
```

Your own blocks must describe the same individuals in the same row order; preprocess them with
`wc.prep_blocks(list_of_arrays)` (column centring and unit total inertia). Weightings:
`"rv"`, `"mfa"`, `"statis"`, `"nvar"`, or any vector of K positive weights relative to the default.

## Main functions

| function | purpose |
|---|---|
| `prep_blocks`, `snv` | preprocessing |
| `block_weights` | usual weightings relative to the default |
| `cpca`, `saliences` | weighted consensus PCA, relative saliences |
| `ccp` | common and proper components of a subspace |
| `mixing` | mixing indices of principal components |
| `presence`, `four_parts`, `reweighting_bound` | ingredients of the rule |
| `canonical_factorization` | `D_k` and `V_k` on a given basis |
| `wecops` | one analysis read in both bases |
| `compare_weighting` | equivalence rule, weighting versus default |
| `plot_fiche` | diagnostic sheet |
| `load_wines`, `load_butters`, `load_simulated` | datasets of the article |

## Reproducing the article

`python examples/reproduce_article.py` reproduces the 28 verdicts (four weightings, seven cases).

## Reference

Hanafi M. When do two weightings give the same interpretation in common and proper structures
between several blocks of data? (submitted).

## Author and maintainer

Mohamed Hanafi <mohammed44hanafi@gmail.com>. Licence: MIT.
