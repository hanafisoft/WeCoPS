"""Common and proper components: concentration (INDORT criterion restricted to the retained
subspace), data-driven target at threshold alpha, and rotation towards the target.

All rotations are carried out in the retained subspace: the components are ``V = U Q`` with
``Q`` orthogonal (R x R), and the saliences are ``s_kh = q_h' (U' W_k U) q_h``.
"""
from __future__ import annotations

import numpy as np

__all__ = ["concentration", "target_rotation", "orient_groups", "label_groups", "ccp", "block_labels", "motif", "label_names"]



def _projected(U, blocks):
    return [(X.T @ U).T @ (X.T @ U) for X in blocks]  # U' X_k X_k' U


def _sal(Q, Wt):
    return np.array([np.einsum("ih,ij,jh->h", Q, Wk, Q) for Wk in Wt])


def concentration(Wt, nstart=10, maxit=1000, tol=1e-12, seed=0):
    """Maximise the concentration ``C(Q) = sum_k sum_h s_kh^2`` over orthogonal Q.

    Polar (monotone) algorithm, started from the identity and ``nstart - 1`` random rotations;
    the best solution is kept. ``Wt`` is the list of projected blocks ``U' W_k U``.
    """
    R = Wt[0].shape[0]
    rng = np.random.default_rng(seed)
    best_f, best_Q = -np.inf, None
    for st in range(nstart):
        Q = np.eye(R) if st == 0 else np.linalg.qr(rng.standard_normal((R, R)))[0]
        for _ in range(maxit):
            S = _sal(Q, Wt)
            G = sum((Wk @ Q) * S[k] for k, Wk in enumerate(Wt))
            A, _, Bt = np.linalg.svd(G)
            Qn = A @ Bt
            if np.abs(Qn - Q).max() < tol:
                Q = Qn
                break
            Q = Qn
        f = np.sum(_sal(Q, Wt) ** 2)
        if f > best_f + 1e-12:
            best_f, best_Q = f, Q
    return best_Q


def target_rotation(Q, Wt, T, sweeps=100, tol=1e-12):
    """Rotate towards a binary target T (K x R): minimise the inertia of each block on the
    components it should not carry, ``Phi = sum_h sum_{k: T_kh = 0} s_kh``, by sweeps of plane
    rotations. In the plane (q_i, q_j), Phi is ``const + B cos 2t + C sin 2t``; the optimal angle
    is ``t = atan2(-C, -B) / 2`` (no rotation when the criterion is flat in that plane)."""
    Q = Q.copy()
    R = Q.shape[1]
    for _ in range(sweeps):
        moved = 0.0
        for i in range(R):
            for j in range(i + 1, R):
                qi, qj = Q[:, i].copy(), Q[:, j].copy()
                B = C = 0.0
                for k in range(len(Wt)):
                    zi, zj = T[k, i] == 0, T[k, j] == 0
                    if not (zi or zj) or (zi and zj):
                        continue  # both or none excluded: Phi does not depend on the angle
                    a, b, d = qi @ Wt[k] @ qi, qj @ Wt[k] @ qj, qi @ Wt[k] @ qj
                    sg = 1.0 if zi else -1.0
                    B += sg * (a - b) / 2
                    C += sg * d
                if B * B + C * C < tol * tol:
                    continue
                t = 0.5 * np.arctan2(-C, -B)
                c, s = np.cos(t), np.sin(t)
                Q[:, i], Q[:, j] = c * qi + s * qj, -s * qi + c * qj
                moved = max(moved, abs(t))
        if moved < 1e-10:
            break
    return Q


def block_labels(S, alpha=0.05, previous=None, band=None):
    """Binary labels (K x R): block k carries component h if ``s_kh >= alpha``.

    With ``previous`` and ``band = (lo, hi)``, a membership changes only if the salience leaves
    the band (hysteresis): set to 1 when ``s >= hi``, to 0 when ``s < lo``, unchanged otherwise.
    """
    S = np.asarray(S)
    if previous is None or band is None:
        return (S >= alpha).astype(int)
    lo, hi = band
    T = np.array(previous, dtype=int).copy()
    T[S >= hi] = 1
    T[S < lo] = 0
    return T


def label_names(T, block_names=None):
    """Label of each component: 'common 1-2', 'proper 3' or 'empty'."""
    K, R = T.shape
    names = block_names or [str(k + 1) for k in range(K)]
    out = []
    for h in range(R):
        b = [names[k] for k in range(K) if T[k, h]]
        out.append("common " + "-".join(b) if len(b) > 1 else ("proper " + b[0] if b else "empty"))
    return out


def motif(T, block_names=None):
    """Motif of a solution: the sorted list of its component labels."""
    return sorted(label_names(T, block_names))


def label_groups(T):
    """Groups of components sharing a label: list of (label tuple, column indices)."""
    labs = [tuple(int(x) for x in T[:, j]) for j in range(T.shape[1])]
    out = []
    for L in sorted(set(labs), key=lambda L: labs.index(L)):
        out.append((L, [j for j, l in enumerate(labs) if l == L]))
    return out


def orient_groups(Q, Wt, T, U, weights=None):
    """Orientation of the components inside each labelled space.

    The rotation towards the target determines, for each label L, the space spanned by the
    components of that label, but not their orientation inside it. In each group of two or more
    components, the components are replaced by the axes of the multiple co-inertia analysis
    (Chessel & Hanafi, 1996) of the blocks of the label restricted to that space: eigenvectors of
    ``sum_{k in L} w_k Q_L' (U' W_k U) Q_L``, in decreasing order of eigenvalues. Saliences of
    each group, labels and motif are unchanged. Signs are fixed so that the largest entry (in
    absolute value) of each component ``v_h = U q_h`` is positive.
    """
    Q = Q.copy()
    K = len(Wt)
    w = np.ones(K) if weights is None else np.asarray(weights, dtype=float)
    for L, idx in label_groups(T):
        if len(idx) < 2:
            continue
        ks = [k for k in range(K) if L[k] == 1] or list(range(K))
        QL = Q[:, idx]
        M = sum(w[k] * QL.T @ Wt[k] @ QL for k in ks)
        val, vec = np.linalg.eigh((M + M.T) / 2)
        Q[:, idx] = QL @ vec[:, ::-1]
    V = U @ Q
    sgn = np.sign(V[np.argmax(np.abs(V), axis=0), np.arange(V.shape[1])])
    sgn[sgn == 0] = 1
    return Q * sgn


def ccp(blocks, U, alpha=0.05, nstart=10, seed=0, maxit=20, sweeps=100, weights=None):
    """Common and proper components of the subspace spanned by U.

    1. Concentration (INDORT criterion restricted to the subspace).
    2. Target: ``T = (s_kh >= alpha)``.
    3. Rotation towards the target; steps 2-3 are repeated until the target is stable.
       This determines the space spanned by each group of components sharing a label.
    4. Orientation inside each labelled space (:func:`orient_groups`), which makes the basis
       unique (up to signs) when the eigenvalues of each group are distinct. The labels returned
       are those of the oriented components (threshold alpha on their saliences).

    Returns
    -------
    dict with ``V`` (n x R), ``Q`` (V = U Q), ``T`` (K x R labels), ``S`` (K x R saliences).
    """
    Wt = _projected(U, blocks)
    Q = concentration(Wt, nstart=nstart, seed=seed)
    T = block_labels(_sal(Q, Wt), alpha)
    for _ in range(maxit):
        Q = target_rotation(Q, Wt, T, sweeps=sweeps)
        T2 = block_labels(_sal(Q, Wt), alpha)
        if np.array_equal(T2, T):
            break
        T = T2
    Q = orient_groups(Q, Wt, T, U, weights)
    S = _sal(Q, Wt)
    T = block_labels(S, alpha)  # labels of the final (oriented) components
    return dict(V=U @ Q, Q=Q, T=T, S=S)


def ccp_from(blocks, U, Q0, T_prev, band, maxit=20, sweeps=100, weights=None):
    """Common and proper components of U starting from rotation Q0, with labels compared to
    T_prev through a tolerance band (used to compare a weighting with the default)."""
    Wt = _projected(U, blocks)
    Q = Q0.copy()
    T = block_labels(_sal(Q, Wt), previous=T_prev, band=band)
    for _ in range(maxit):
        Q = target_rotation(Q, Wt, T, sweeps=sweeps)
        T2 = block_labels(_sal(Q, Wt), previous=T_prev, band=band)
        if np.array_equal(T2, T):
            break
        T = T2
    Q = orient_groups(Q, Wt, T, U, weights)
    S = _sal(Q, Wt)
    T = block_labels(S, previous=T_prev, band=band)  # labels of the final (oriented) components
    return dict(V=U @ Q, Q=Q, T=T, S=S)
