"""WeCoPS: Weighting CPCA for Common and Proper Structures.

Does a block weighting lead to the same interpretation as the default weighting in consensus
PCA? WeCoPS reads the retained subspace in two bases (principal components; common and proper
components), measures how principal components mix common and proper components, and applies the
equivalence rule (same motif, presence of every default component above a threshold).
"""
from .core import prep_blocks, snv, block_weights, cpca, saliences, WEIGHTINGS
from .rotation import ccp, concentration, target_rotation, block_labels, label_names, motif
from .compare import (mixing, presence, four_parts, canonical_factorization, reweighting_bound,
                      compare_weighting, WeCoPSResult)
from .analysis import wecops
from .plot import plot_fiche
from .datasets import load_wines, load_butters, load_simulated

__version__ = "0.2.1"
