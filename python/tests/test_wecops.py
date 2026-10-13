"""Tests of WeCoPS against the results of the article."""
import numpy as np
import pytest

import wecops as wc


def test_preprocessing_unit_inertia():
    blocks = wc.prep_blocks([np.arange(12.).reshape(4, 3) ** 2, np.random.default_rng(0).normal(size=(4, 5))])
    for X in blocks:
        assert np.allclose(X.mean(0), 0)
        assert np.isclose(np.linalg.norm(X), 1)


def test_saliences_sum_to_one_on_full_basis():
    X = wc.load_simulated(1)["blocks"]
    a = wc.cpca(X, None, 5)
    full = np.linalg.eigh(sum(x @ x.T for x in X))[1]
    assert np.allclose(wc.saliences(full, X).sum(1), 1)
    assert np.all(wc.saliences(a["U"], X).sum(1) <= 1 + 1e-12)


def test_weights_scale_invariance():
    X = wc.load_wines()["blocks"]
    U1 = wc.cpca(X, [1, 2, 3, 4], 3)["U"]
    U2 = wc.cpca(X, [10, 20, 30, 40], 3)["U"]
    assert np.allclose(np.abs(U1.T @ U2), np.eye(3), atol=1e-10)


def test_simulated_1_two_bases():
    a = wc.wecops(wc.load_simulated(1)["blocks"], 5)
    # principal components mix the common structure of blocks 1-2 with the proper structure of block 3
    assert a.motif_pc == ["common 1-2", "common 1-2-3", "common 1-2-3", "proper 1", "proper 2"]
    # common and proper components recover the true motif
    assert a.motif == ["common 1-2", "common 1-2", "proper 1", "proper 2", "proper 3"]
    assert a.mixing["mu_star"].max() == pytest.approx(0.16, abs=0.01)
    # canonical factorisation: d_kh^2 equals the salience
    for k, f in enumerate(a.factorization["ccp"]):
        assert np.allclose(np.diag(f["D"]) ** 2, a.saliences[k])


def test_coincidence_when_no_coupling():
    # blocks with orthogonal supports: no mixing between different labels
    rng = np.random.default_rng(1)
    Z = rng.normal(size=(20, 3))
    Z = np.linalg.qr(Z - Z.mean(0))[0]  # centred, orthonormal structures
    X1 = np.outer(Z[:, 0], rng.normal(size=6)) * 3 + np.outer(Z[:, 1], rng.normal(size=6))
    X2 = np.outer(Z[:, 2], rng.normal(size=5))
    X = wc.prep_blocks([X1, X2])
    a = wc.wecops(X, 3)
    assert a.mixing["mu_star"].max() < 1e-8


@pytest.mark.parametrize("case,R,expected", [
    (1, 5, {"rv": (True, False), "mfa": (True, False), "statis": (True, False), "nvar": (True, True)}),
    (3, 5, {"rv": (False, None), "mfa": (False, None), "statis": (True, None), "nvar": (True, None)}),
])
def test_verdicts_simulated(case, R, expected):
    X = wc.load_simulated(case)["blocks"]
    for w, (eq, pc_same) in expected.items():
        r = wc.compare_weighting(X, R, w)
        assert r.equivalent == eq, w
        if pc_same is not None:
            assert r.pc_same == pc_same, w


def test_verdicts_butters():
    b = wc.load_butters()
    r = wc.compare_weighting(b["blocks"], 4, "nvar", block_names=b["block_names"])
    assert not r.equivalent
    assert r.kappa.min() == pytest.approx(0.37, abs=0.01)
    assert wc.compare_weighting(b["blocks"], 4, "mfa").equivalent


def test_reweighting_threshold():
    X = wc.load_simulated(1)["blocks"]
    r = wc.compare_weighting(X, 5, "mfa")
    assert r.bound["guaranteed"]
    assert np.all(r.kappa_components >= r.bound["kappa_lower_bound"] - 1e-10)


def test_labels_after_orientation():
    b = wc.load_wines()["blocks"]
    r = wc.compare_weighting(b, 3, "rv", alpha=0.01)
    assert r.equivalent
    assert np.array_equal(r.T_default, (r.saliences_default >= 0.01).astype(int))
