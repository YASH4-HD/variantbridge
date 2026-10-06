import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from scipy import stats

from variantbridge import association as A


@pytest.fixture
def data():
    r = np.random.default_rng(1)
    n, m = 250, 8
    G = r.binomial(2, 0.3, (n, m)).astype(float)
    G[3, 2] = np.nan
    G[10:14, 5] = np.nan
    C = r.normal(size=(n, 3))
    y = 0.4 * np.nan_to_num(G[:, 0]) + C @ [0.5, -0.2, 0.3] + r.normal(size=n)
    return G, C, y


def test_linear_matches_statsmodels_including_missing(data):
    G, C, y = data
    res = A.scan_linear(G, y, C)
    for j in range(G.shape[1]):
        ok = np.isfinite(G[:, j])
        f = sm.OLS(y[ok], sm.add_constant(np.column_stack([C[ok], G[ok, j]]))).fit()
        assert res.effect[j] == pytest.approx(f.params[-1], rel=1e-9)
        assert res.se[j] == pytest.approx(f.bse[-1], rel=1e-9)
        assert res.p[j] == pytest.approx(f.pvalues[-1], rel=1e-7)
        lo, hi = np.asarray(f.conf_int())[-1]
        assert res.ci_low[j] == pytest.approx(lo, rel=1e-8) and res.ci_high[j] == pytest.approx(hi, rel=1e-8)
        assert res.n[j] == ok.sum()


def test_vectorised_scan_equals_single_fit(data):
    G, C, y = data
    res = A.scan_linear(G, y, C)
    one = A.fit_ols(y, G[:, 4], C, variant="v4")
    for k in ("effect", "se", "p", "ci_low", "ci_high", "df", "n"):
        assert res.loc[4, k] == pytest.approx(one[k], rel=1e-10)


def test_logistic_matches_statsmodels(data):
    G, C, y = data
    yb = (y > np.median(y)).astype(float)
    res = A.scan_logistic(G, yb, C)
    for j in range(G.shape[1]):
        ok = np.isfinite(G[:, j])
        f = sm.Logit(yb[ok], sm.add_constant(np.column_stack([C[ok], G[ok, j]]))).fit(disp=0, tol=1e-12)
        assert res.effect[j] == pytest.approx(f.params[-1], rel=1e-6)
        assert res.se[j] == pytest.approx(f.bse[-1], rel=1e-6)
        assert res.p[j] == pytest.approx(f.pvalues[-1], rel=1e-5)
        assert res.odds_ratio[j] == pytest.approx(np.exp(f.params[-1]), rel=1e-6)


def test_p_value_is_covariate_adjusted_not_univariate_correlation():
    """Regression test for the v0.1 flaw: p came from an unadjusted correlation while the effect
    came from an adjusted model. With a confounder driving both x and y the two differ a lot."""
    r = np.random.default_rng(5)
    n = 400
    conf = r.normal(size=n)
    x = conf + 0.3 * r.normal(size=n)
    y = 2.0 * conf + r.normal(size=n)  # no direct x effect
    res = A.scan_linear(x[:, None], y, conf[:, None])
    _, p_unadjusted = stats.pearsonr(x, y)
    assert p_unadjusted < 1e-10          # naive v0.1-style statistic is wildly significant
    assert res.p[0] > 0.01               # adjusted model correctly shows no association
    f = sm.OLS(y, sm.add_constant(np.column_stack([conf, x]))).fit()
    assert res.p[0] == pytest.approx(f.pvalues[-1], rel=1e-7)


def test_consistency_checker_passes_and_detects_tampering(data):
    G, C, y = data
    res = A.scan_linear(G, y, C)
    assert A.check_statistical_consistency(res) == []
    bad = res.copy()
    bad.loc[0, "p"] = 0.123456                       # p no longer derived from effect/se/df
    assert any("p-value" in s for s in A.check_statistical_consistency(bad))
    bad = res.copy()
    bad.loc[1, "ci_high"] += 0.5
    assert any("CI" in s for s in A.check_statistical_consistency(bad))
    with pytest.raises(AssertionError):
        A.assert_statistical_consistency(bad)
    yb = (y > np.median(y)).astype(float)
    rl = A.scan_logistic(G, yb, C)
    assert A.check_statistical_consistency(rl) == []


def test_degenerate_and_missing_handling():
    r = np.random.default_rng(0)
    G = np.column_stack([np.ones(50), r.binomial(2, .4, 50)]).astype(float)
    y = r.normal(size=50)
    res = A.scan_linear(G, y)
    assert res.status[0] == "degenerate" and np.isnan(res.p[0])
    assert res.status[1] == "ok"
    with pytest.raises(ValueError):
        A.scan_linear(G, np.where(np.arange(50) == 3, np.nan, y))


def test_logistic_separation_is_reported_not_hidden():
    n = 60
    x = np.r_[np.zeros(30), np.ones(30)]
    y = x.copy()                                      # perfect separation
    r = A.fit_logistic(y, x, variant="sep")
    assert r["status"] in ("separation", "no_convergence") and np.isnan(r["p"])


def test_dispatcher_auto_selects_model(data):
    G, C, y = data
    assert A.scan(G, y, C).model.iloc[0] == "OLS"
    assert A.scan(G, (y > 0).astype(float), C).model.iloc[0] == "logistic-Wald"
