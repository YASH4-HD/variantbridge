import numpy as np
import pandas as pd
import pytest

from variantbridge.amr import cohorts, phenotypes as P, prediction as PR, structure


def test_lineage_blocked_folds_hold_out_whole_lineages_deterministically():
    groups = np.repeat([f"ST{i}" for i in range(12)], [30, 25, 20, 18, 15, 12, 10, 9, 8, 6, 5, 4])
    f1 = PR.lineage_blocked_folds(groups, 5); f2 = PR.lineage_blocked_folds(groups, 5)
    assert (f1 == f2).all()
    for g in set(groups):
        assert len(set(f1[groups == g])) == 1, "a lineage was split across folds"
    assert set(f1) == set(range(5))


def test_too_few_lineages_refuses():
    with pytest.raises(ValueError, match="lineage groups"):
        PR.lineage_blocked_folds(["a", "b", "a"], 5)


def test_random_folds_do_share_lineages():
    groups = np.repeat([f"ST{i}" for i in range(8)], 25)
    f = PR.random_folds(len(groups), 5, seed=0)
    assert any(len(set(f[groups == g])) > 1 for g in set(groups))


def test_interval_mae_is_censoring_aware():
    lo = np.array([-np.inf, 2.0, 0.0]); hi = np.array([-2.0, np.inf, 0.0])
    assert PR.interval_mae(np.array([-5.0, 5.0, 0.0]), lo, hi) == 0.0      # inside the censored/exact interval
    assert PR.interval_mae(np.array([0.0, 0.0, 1.0]), lo, hi) == pytest.approx((2 + 2 + 1) / 3)


def test_ridge_tobit_recovers_sign_and_handles_censoring():
    rng = np.random.default_rng(0)
    n = 800; X = rng.integers(0, 2, (n, 6)).astype(float)
    y = -1 + 2.5 * X[:, 0] + rng.normal(0, 1, n)
    lo = np.where(y <= -2, -np.inf, np.where(y >= 2, 2.0, np.round(y))); hi = np.where(y <= -2, -2.0, np.where(y >= 2, np.inf, np.round(y)))
    m = PR.fit_ridge_tobit(X, lo, hi, 1.0)
    assert m["converged"] and m["coef"][0] > 1.5 and np.abs(m["coef"][1:]).max() < 0.5


def _lineage_only_scenario(seed=0, n_lin=24, per=30, n_feat=150):
    """Phenotype is a property of the LINEAGE only (random per lineage); features are clonal lineage signatures.
    A random split can memorise signature -> phenotype; held-out lineages carry no such information."""
    rng = np.random.default_rng(seed)
    st = np.repeat([f"ST{i}" for i in range(n_lin)], per)
    high = {f"ST{i}": (i % 2 == 0) for i in range(n_lin)}
    base = {f"ST{i}": rng.random(n_feat) < 0.5 for i in range(n_lin)}      # unrelated to `high`
    X = np.array([np.where(rng.random(n_feat) < 0.03, ~base[s], base[s]) for s in st], float)
    y = np.where([high[s] for s in st], 3.0, -3.0) + rng.normal(0, 0.3, len(st))
    ms = []
    for v in y:
        r = int(np.round(v)); ms.append("<=0.25" if r <= -2 else (">4" if r > 2 else f"{2.0**r:g}"))
    ph = P.parse_mic_table(pd.DataFrame({"measurement": ms}))
    ph["cohort_tier"] = "discovery"
    return X, ph, st


def test_lineage_blocked_cv_is_lower_than_random_split_leakage_lesson(cfg):
    X, ph, st = _lineage_only_scenario()
    ids = [f"f{i}" for i in range(X.shape[1])]
    blocked = PR.nested_lineage_blocked_cv(X, ph, st, ids, cfg)
    rnd = PR.random_split_leakage_demo(X, ph, st, ids, cfg)
    gap = PR.leakage_gap(blocked, rnd)
    assert rnd["cv_strategy"] == "random_demo_leakage_only" and blocked["cv_strategy"] == "nested_lineage_blocked"
    assert rnd["mean_enet_extreme_auc"] > 0.95 and blocked["mean_enet_extreme_auc"] < 0.75 and gap["gap"] > 0.2
    assert (blocked["per_fold"]["lineage_overlap_train_test"] == 0).all()           # no lineage shared train/test
    assert (rnd["per_fold"]["lineage_overlap_train_test"] > 0).any()
    assert "never be reported as performance" in gap["note"]


def test_blocked_cv_reports_honest_signal_when_signal_is_isolate_level(syn, disc, cfg):
    blocked = PR.nested_lineage_blocked_cv(disc["G"], disc["pheno"], disc["st"].to_numpy(), disc["ids"], cfg)
    assert blocked["mean_enet_extreme_auc"] > 0.85 and np.isfinite(blocked["mean_interval_mae"])
    assert blocked["importance_warning"].startswith("Feature importance reflects predictive utility")
    top = blocked["importance"].head(5)["feature_id"].tolist()
    assert "gyrA_S83L" in top
    assert (blocked["per_fold"]["lineage_overlap_train_test"] == 0).all()


def test_screening_uses_training_data_only(disc, cfg):
    """Selected feature indices come from the training fold: changing TEST labels cannot change them."""
    X, ph = disc["G"], disc["pheno"]
    tr = np.arange(len(ph)) < 350
    a = PR._screen(X[tr], ph["ordinal_code"].to_numpy()[tr], 20)
    ph2 = ph.copy(); ph2.loc[~tr, "ordinal_code"] = 0
    b = PR._screen(X[tr], ph2["ordinal_code"].to_numpy()[tr], 20)
    assert (a == b).all()


def test_discovery_to_replication_is_primary_and_tier_guarded(syn, disc, cfg):
    k = syn["kept"]; rp = k[k["cohort_tier"] == "replication"].reset_index(drop=True)
    Xr = syn["features"].loc[rp["genome_id"]].to_numpy()
    out = PR.discovery_to_replication(disc["G"], disc["pheno"], disc["st"].to_numpy(), disc["ids"], Xr, rp, disc["ids"], cfg)
    assert out["cv_strategy"] == "discovery_to_replication" and out["enet_extreme_auc"] > 0.8 and "UNVERIFIED comparability" in out["caveat"]
    with pytest.raises(cohorts.CohortSeparationError):
        PR.discovery_to_replication(Xr, rp, np.array(["x"] * len(rp)), disc["ids"], disc["G"], disc["pheno"], disc["ids"], cfg)   # swapped


def test_exploratory_data_cannot_produce_model_metrics(syn, cfg):
    k = syn["kept"]; ex = k[k["cohort_tier"] == "exploratory"].reset_index(drop=True)
    X = syn["features"].loc[ex["genome_id"]].to_numpy()
    with pytest.raises(cohorts.CohortSeparationError):
        PR.nested_lineage_blocked_cv(X, ex, syn["st"].loc[ex["genome_id"]].to_numpy(), list(syn["features"].columns), cfg)
    with pytest.raises(cohorts.CohortSeparationError):
        PR.random_split_leakage_demo(X, ex, syn["st"].loc[ex["genome_id"]].to_numpy(), list(syn["features"].columns), cfg)
    d = PR.describe_exploratory(ex)
    assert "descriptive only" in d["use"] and "auc" not in " ".join(d).lower()


def test_structure_helpers_deterministic(syn, cfg):
    st = syn["st"].loc[syn["kept"]["genome_id"]].reset_index(drop=True)
    g = structure.lineage_groups(st, 10, 5)
    assert g.nunique() <= 6 and set(g) - set(st) <= {"other_small"}
    L1, n1 = structure.lineage_design(st, cfg); L2, n2 = structure.lineage_design(st, cfg)
    assert (L1 == L2).all() and n1 == n2
    X = syn["features"].to_numpy()[:50]
    D = structure.jaccard_distance(X); assert np.allclose(D, D.T) and np.allclose(np.diag(D), 0)
    p1, p2 = structure.pcoa(D, 3), structure.pcoa(D, 3); assert np.allclose(p1, p2)
