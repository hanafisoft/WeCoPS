"""Mixing of principal components, presence, four-part decomposition, canonical factorisation
and the equivalence rule between a weighting and the default weighting."""
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment

from .core import _as_blocks, block_weights, cpca, saliences
from .rotation import ccp, ccp_from, block_labels, label_names, label_groups, motif

__all__ = ["mixing", "presence", "group_presence", "four_parts", "canonical_factorization", "reweighting_bound",
           "compare_weighting", "WeCoPSResult"]


def mixing(U, V, T):
    """Mixing of the principal components U with respect to the common and proper components V.

    ``U = V Gamma`` with ``Gamma = V'U``; ``gamma_jh^2`` is the share of u_h carried by v_j.

    Returns
    -------
    dict with ``gamma`` (R x R), ``mu`` (``1 - max_j gamma_jh^2``) and ``mu_star`` (same, after
    pooling the components v_j that share a label: mixing between different labels only).
    """
    G = V.T @ U
    G2 = G ** 2
    mu = 1 - G2.max(axis=0)
    labs = [tuple(T[:, j]) for j in range(T.shape[1])]
    groups = sorted(set(labs))
    pooled = np.array([[G2[[j for j, l in enumerate(labs) if l == L], h].sum() for h in range(U.shape[1])]
                       for L in groups])
    mu_star = 1 - pooled.max(axis=0)
    return dict(gamma=G, mu=np.clip(mu, 0, None), mu_star=np.clip(mu_star, 0, None))


def presence(U_alt, V0):
    """Presence ``kappa_j = ||U_alt' v_j||^2`` of each component v_j of V0 in the subspace of U_alt."""
    return np.sum((U_alt.T @ V0) ** 2, axis=0)


def group_presence(U_alt, V0, T0, block_names=None):
    """Presence of each labelled space of V0 in the subspace of U_alt.

    For a group of components sharing a label, only the space E_L they span is determined by the
    rotation towards the target; its presence is the squared cosine of the largest principal angle
    between E_L and the subspace of U_alt, i.e. the smallest eigenvalue of ``V_L' U_alt U_alt' V_L``.
    For a group of one component it equals ``kappa_j``.

    Returns
    -------
    list of dicts with ``label`` (name), ``components`` (indices) and ``kappa``.
    """
    out = []
    P = U_alt.T @ V0
    for L, idx in label_groups(T0):
        G = P[:, idx].T @ P[:, idx]
        out.append(dict(label=label_names(T0[:, idx[:1]], block_names)[0], components=idx,
                        kappa=float(np.linalg.eigvalsh((G + G.T) / 2).min())))
    return out


def four_parts(S, T):
    """Decomposition of the inertia of each block (unit inertia) into four parts:
    common components, its own proper components, proper components of other blocks (and empty
    components), and the residual outside the retained subspace. Returns a (K, 4) array."""
    K, R = S.shape
    P = np.zeros((K, 4))
    nb = T.sum(axis=0)
    for k in range(K):
        for h in range(R):
            if nb[h] >= 2:
                P[k, 0] += S[k, h]
            elif T[k, h]:
                P[k, 1] += S[k, h]
            else:
                P[k, 2] += S[k, h]
        P[k, 3] = 1 - S[k].sum()
    return P


def canonical_factorization(blocks, B):
    """Canonical factorisation ``X_k ~ B D_k V_k'`` of each block on the basis B (n x R).

    ``d_kh = ||X_k' b_h||`` (so that ``d_kh^2`` is the salience s_kh) and
    ``V_k = X_k' B D_k^{-1}`` (unit-norm loadings).

    Returns
    -------
    list of dicts with ``D`` (R x R diagonal) and ``V`` (J_k x R), one per block.
    """
    out = []
    for X in blocks:
        L = X.T @ B
        d = np.linalg.norm(L, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            Vk = np.where(d > 0, L / d, 0.0)
        out.append(dict(D=np.diag(d), V=Vk))
    return out


def reweighting_bound(blocks, R, weights, V0=None):
    """Reweighting threshold (sufficient condition for equivalence).

    With ``A = max w / min w`` and the eigenvalues of ``M(1) = sum_k X_k X_k'``: if
    ``A < lambda_R / lambda_{R+1}``, no reweighting can change the retained subspace much, and each
    component v_j of V0 satisfies ``kappa_j >= (rho_j / A - lambda_{R+1}) / (lambda_1 - lambda_{R+1})``,
    with ``rho_j = v_j' M(1) v_j``.
    """
    blocks = _as_blocks(blocks)
    w = block_weights(blocks, weights)
    M1 = sum(X @ X.T for X in blocks)
    lam = np.clip(np.linalg.eigvalsh(M1)[::-1], 0, None)
    ratio = lam[R - 1] / lam[R] if lam[R] > 0 else np.inf
    A = w.max() / w.min()
    out = dict(A=A, ratio=ratio, guaranteed=bool(A < ratio))
    if V0 is not None:
        rho = np.einsum("ij,ik,kj->j", V0, M1, V0)
        out["kappa_lower_bound"] = (rho / A - lam[R]) / (lam[0] - lam[R]) if out["guaranteed"] else None
    return out


def _match(U0, U):
    C = np.abs(U0.T @ U)
    r, c = linear_sum_assignment(-C)
    Um = U[:, c] * np.sign(np.sum(U0 * U[:, c], axis=0))
    return Um, C[r, c] ** 2


class WeCoPSResult(dict):
    """Result of :func:`compare_weighting` (a dict with attribute access and a readable print)."""

    __getattr__ = dict.__getitem__

    def __repr__(self):  # pragma: no cover - display only
        return self.summary()

    def summary(self):
        lines = [
            f"WeCoPS comparison: weighting '{self['weighting_name']}' versus default (R = {self['R']}, alpha = {self['alpha']:g})",
            f"  verdict: {'EQUIVALENT' if self['equivalent'] else 'INTERPRETATION CHANGED'}"
            f" (same motif: {self['same_motif']}, min presence: {min(self['kappa']):.3f}, threshold {self['kappa_min']:g})",
            f"  default motif:     {', '.join(self['motif_default'])}",
            f"  alternative motif: {', '.join(self['motif_alt'])}",
            "  presence of each labelled space of the default: " + ", ".join(
                f"{g['label']}" + (f" ({len(g['components'])} comp.)" if len(g['components']) > 1 else "")
                + f": {g['kappa']:.2f}" for g in self['groups']),
            f"  weight ratio A = {self['bound']['A']:.3g}; lambda_R/lambda_R+1 = {self['bound']['ratio']:.3g}"
            + ("  -> below the reweighting threshold: no component can leave the subspace entirely" if self['bound']['guaranteed'] else ""),
            f"  principal-component reading: {'same' if self['pc_same'] else 'changed'}"
            f" (max mixing mu* default {self['mixing_default']['mu_star'].max():.2f}, alternative {self['mixing_alt']['mu_star'].max():.2f})",
        ]
        if self['equivalent'] and not self['pc_same']:
            lines.append("  note: reading the principal components one by one would wrongly signal a change")
        return "\n".join(lines)


def compare_weighting(blocks, R, weights, alpha=0.05, kappa_min=0.9, band=0.02, nstart=10, seed=0,
                      block_names=None):
    """Does a weighting lead to the same interpretation as the default weighting?

    Parameters
    ----------
    blocks : list of (n, J_k) arrays, preprocessed with :func:`wecops.prep_blocks`.
    R : int
        Number of retained components (fixed).
    weights : str or array-like
        Alternative weighting: ``"rv"``, ``"mfa"``, ``"statis"``, ``"nvar"`` or K positive weights
        relative to the default weighting.
    alpha : float
        Salience threshold for labels.
    kappa_min : float
        Presence threshold of the equivalence rule.
    band : float
        Half-width of the tolerance band around alpha used to compare labels (``0`` = strict).

    Returns
    -------
    WeCoPSResult
        Verdict, motifs, presence, four-part decompositions, mixing of principal components,
        reweighting threshold, and canonical factorisations (D_k, V_k) of both analyses on both bases.
    """
    blocks = _as_blocks(blocks)
    K = len(blocks)
    w = block_weights(blocks, weights)
    name = weights if isinstance(weights, str) else "user"
    lo, hi = (max(alpha - band, alpha / 2), alpha + band) if band > 0 else (alpha, alpha)

    a0 = cpca(blocks, None, R)
    U0 = a0["U"]
    c0 = ccp(blocks, U0, alpha, nstart=nstart, seed=seed)
    V0, T0 = c0["V"], c0["T"]

    ag = cpca(blocks, w, R)
    Ug = ag["U"]
    P_, _, Pt = np.linalg.svd(Ug.T @ V0)
    cg = ccp_from(blocks, Ug, P_ @ Pt, T0, (lo, hi), weights=w)
    Vg, Tg = cg["V"], cg["T"]

    kap_c = presence(Ug, V0)
    groups = group_presence(Ug, V0, T0, block_names)
    kap = np.array([g["kappa"] for g in groups])
    m0, mg = motif(T0, block_names), motif(Tg, block_names)
    same = m0 == mg
    equivalent = bool(same and kap.min() >= kappa_min)

    # usual reading of the principal components, one by one
    Um, cos2 = _match(U0, Ug)
    S_pc0, S_pcg = saliences(U0, blocks), saliences(Um, blocks)
    pc_same = bool(cos2.min() >= kappa_min and np.array_equal(S_pc0 >= alpha, S_pcg >= alpha))

    S0, Sg = saliences(V0, blocks), saliences(Vg, blocks)
    parts0, partsg = four_parts(S0, T0), four_parts(Sg, Tg)

    return WeCoPSResult(
        weighting_name=name, weights=w, R=R, alpha=alpha, kappa_min=kappa_min, band=(lo, hi),
        equivalent=equivalent, same_motif=bool(same), kappa=kap, groups=groups, kappa_components=kap_c,
        motif_default=m0, motif_alt=mg,
        labels_default=label_names(T0, block_names), labels_alt=label_names(Tg, block_names),
        T_default=T0, T_alt=Tg, V_default=V0, V_alt=Vg, U_default=U0, U_alt=Ug,
        saliences_default=S0, saliences_alt=Sg,
        parts_default=parts0, parts_alt=partsg, parts_change=partsg - parts0,
        eigenvalues_default=a0["explained"], eigenvalues_alt=ag["explained"],
        mixing_default=mixing(U0, V0, T0), mixing_alt=mixing(Ug, Vg, Tg),
        pc_same=pc_same, pc_cos2=cos2,
        bound=reweighting_bound(blocks, R, w, V0),
        factorization=dict(
            default_ccp=canonical_factorization(blocks, V0), default_pc=canonical_factorization(blocks, U0),
            alt_ccp=canonical_factorization(blocks, Vg), alt_pc=canonical_factorization(blocks, Ug)),
        block_names=block_names or [f"block {k + 1}" for k in range(K)],
    )
