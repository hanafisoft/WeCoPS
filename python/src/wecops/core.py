"""Block preprocessing, block weightings and weighted consensus PCA (CPCA)."""
from __future__ import annotations

import numpy as np

__all__ = ["prep_blocks", "snv", "block_weights", "cpca", "saliences", "WEIGHTINGS"]

WEIGHTINGS = ("default", "rv", "mfa", "statis", "nvar")


def _as_blocks(blocks):
    out = [np.asarray(X, dtype=float) for X in blocks]
    if len(out) < 2:
        raise ValueError("at least two blocks are required")
    n = out[0].shape[0]
    for k, X in enumerate(out):
        if X.ndim != 2:
            raise ValueError(f"block {k + 1} is not a matrix")
        if X.shape[0] != n:
            raise ValueError("all blocks must have the same number of rows (individuals)")
        if not np.all(np.isfinite(X)):
            raise ValueError(f"block {k + 1} contains missing or infinite values")
    return out


def snv(X):
    """Standard normal variate: centre and scale each row (spectrum) to unit standard deviation."""
    X = np.asarray(X, dtype=float)
    return (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, ddof=1, keepdims=True)


def prep_blocks(blocks, center=True, unit_inertia=True):
    """Centre each block by column and scale it to unit total inertia (the default weighting).

    Parameters
    ----------
    blocks : sequence of (n, J_k) arrays
        Blocks describing the same n individuals, in the same row order.
    center : bool
        Centre each variable (column).
    unit_inertia : bool
        Divide each block by its Frobenius norm, so that ``trace(X_k X_k') = 1``.

    Returns
    -------
    list of arrays
    """
    out = []
    for X in _as_blocks(blocks):
        if center:
            X = X - X.mean(axis=0)
        if unit_inertia:
            nrm = np.linalg.norm(X)
            if nrm == 0:
                raise ValueError("a block has zero inertia")
            X = X / nrm
        out.append(X)
    return out


def _statis(blocks):
    W = [X @ X.T for X in blocks]
    nW = np.array([np.linalg.norm(Wk) for Wk in W])
    Wt = [Wk / nk for Wk, nk in zip(W, nW)]
    S = np.array([[np.sum(a * b) for b in Wt] for a in Wt])  # RV coefficients between blocks
    _, vec = np.linalg.eigh(S)
    beta = np.abs(vec[:, -1])
    return beta / nW, beta / beta.sum(), S


def block_weights(blocks, method="default"):
    """Weights of usual weightings, expressed relative to the default weighting.

    The blocks are assumed to be preprocessed with :func:`prep_blocks` (unit total inertia).

    Parameters
    ----------
    method : {"default", "rv", "mfa", "statis", "nvar"} or array-like
        ``"default"``: all weights equal to 1 (unit inertia);
        ``"mfa"``: inverse of the first eigenvalue of each block (multiple factor analysis);
        ``"rv"``: inverse of ``||X_k X_k'||_F`` (RV normalisation of STATIS);
        ``"statis"``: RV normalisation times the STATIS compromise weights;
        ``"nvar"``: inverse of the number of variables.
        An array of K positive numbers is returned unchanged.

    Returns
    -------
    ndarray of shape (K,)
    """
    blocks = _as_blocks(blocks)
    if not isinstance(method, str):
        w = np.asarray(method, dtype=float)
        if w.shape != (len(blocks),) or np.any(w <= 0):
            raise ValueError("weights must be K positive numbers")
        return w
    m = method.lower()
    if m == "default":
        return np.ones(len(blocks))
    if m == "mfa":
        return np.array([1.0 / np.linalg.svd(X, compute_uv=False)[0] ** 2 for X in blocks])
    if m == "rv":
        return np.array([1.0 / np.linalg.norm(X @ X.T) for X in blocks])
    if m == "statis":
        return _statis(blocks)[0]
    if m == "nvar":
        return np.array([1.0 / X.shape[1] for X in blocks])
    raise ValueError(f"unknown weighting '{method}'; use one of {WEIGHTINGS} or a vector")


def cpca(blocks, weights=None, R=2):
    """Weighted consensus PCA.

    The principal components are the R dominant eigenvectors of ``M(w) = sum_k w_k X_k X_k'``,
    that is the principal components of the concatenated table ``[sqrt(w_1) X_1 | ... ]``.
    Only the ratios between weights matter.

    Returns
    -------
    dict with
        ``U``: (n, R) orthonormal principal components;
        ``eigenvalues``: all eigenvalues of ``M(w)`` (weights normalised to sum 1), decreasing;
        ``explained``: eigenvalues as a share of their sum;
        ``weights``: the weights used.
    """
    blocks = _as_blocks(blocks)
    K = len(blocks)
    w = np.ones(K) if weights is None else block_weights(blocks, weights)
    wn = w / w.sum()
    M = sum(a * X @ X.T for a, X in zip(wn, blocks))
    lam, vec = np.linalg.eigh(M)
    lam, vec = lam[::-1], vec[:, ::-1]
    lam = np.clip(lam, 0, None)
    if not 1 <= R < len(lam):
        raise ValueError("R must be between 1 and n - 1")
    U = vec[:, :R].copy()
    return dict(U=U, eigenvalues=lam, explained=lam / lam.sum(), weights=w)


def saliences(U, blocks):
    """Relative saliences ``s_kh = u_h' X_k X_k' u_h``: share of the inertia of block k carried by
    component h (blocks of unit inertia). Returns a (K, R) array."""
    U = np.atleast_2d(np.asarray(U, dtype=float))
    return np.array([np.sum((X.T @ U) ** 2, axis=0) for X in blocks])
