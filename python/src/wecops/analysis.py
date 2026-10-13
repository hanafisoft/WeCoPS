"""One weighted CPCA read in two bases: principal components and common and proper components."""
from __future__ import annotations

from .core import _as_blocks, block_weights, cpca, saliences
from .compare import WeCoPSResult, canonical_factorization, mixing, four_parts
from .rotation import ccp, label_names, block_labels, motif

__all__ = ["wecops"]


def wecops(blocks, R, weights="default", alpha=0.05, nstart=10, seed=0, block_names=None):
    """Weighted CPCA with its principal components and its common and proper components.

    Parameters
    ----------
    blocks : list of (n, J_k) arrays, preprocessed with :func:`wecops.prep_blocks`.
    R : int
        Number of retained components.
    weights : str or array-like
        ``"default"``, ``"rv"``, ``"mfa"``, ``"statis"``, ``"nvar"`` or K positive weights.
    alpha : float
        Salience threshold for labels.

    Returns
    -------
    WeCoPSResult with
        ``U`` (principal components), ``V`` (common and proper components), their saliences and
        labels, the motif, the mixing of the principal components (``gamma``, ``mu``, ``mu_star``),
        the four-part decomposition of each block, and the canonical factorisations
        ``X_k ~ B D_k V_k'`` on both bases.
    """
    blocks = _as_blocks(blocks)
    w = block_weights(blocks, weights)
    a = cpca(blocks, w, R)
    U = a["U"]
    c = ccp(blocks, U, alpha, nstart=nstart, seed=seed, weights=w)
    V, T = c["V"], c["T"]
    S_pc = saliences(U, blocks)
    T_pc = block_labels(S_pc, alpha)
    S = saliences(V, blocks)
    return WeCoPSResult(
        weights=w, R=R, alpha=alpha, eigenvalues=a["explained"],
        U=U, saliences_pc=S_pc, labels_pc=label_names(T_pc, block_names), motif_pc=motif(T_pc, block_names),
        V=V, saliences=S, T=T, labels=label_names(T, block_names), motif=motif(T, block_names),
        mixing=mixing(U, V, T), parts=four_parts(S, T),
        factorization=dict(ccp=canonical_factorization(blocks, V), pc=canonical_factorization(blocks, U)),
        block_names=block_names or [f"block {k + 1}" for k in range(len(blocks))],
    )
