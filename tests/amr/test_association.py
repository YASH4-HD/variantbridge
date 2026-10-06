import inspect

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from scipy import stats

from variantbridge.amr import association as A
from variantbridge.amr import cohorts, phenotypes as P
from variantbridge.amr.config import BreakpointUnverified


def _numeric_check():
    rng = np.random.default_rng(0)
    n = 60
    mu = rng.normal(size=n)
    lo = rng.normal(size=n); hi = lo.copy()
    kind = np.zeros(n, int); kind[:20] = -1; kind[20:40] = 1
    lo[:20] = -np.inf; hi[:20] = mu[:20] + rng.normal(size=20) * .5
    hi[20:40] = np.inf; lo[20:40] = mu[20:40] + rng.normal(size=20) * .5
    return mu, lo, hi, kind


def test_analytic_derivatives_match_numeric():
    mu, lo, hi, kind = _numeric_check()
    s, e = 0.3, 1e-6
    t = A._terms(mu, s, lo, hi, kind)
    f = lambda m, ss, i: A._terms(m, ss, lo, hi, kind)[i]
    assert np.abs((f(mu + e, s, 0) - f(mu - e, s, 0)) / (2 * e) - t[1]).max() < 1e-6
    assert np.abs((f(mu, s + e, 0) - f(mu, s - e, 0)) / (2 * e) - t[2]).max() < 1e-6
    assert np.abs((f(mu + e, s, 1) - f(mu - e, s, 1)) / (2 * e) - t[3]).max() < 1e-6
    assert np.abs((f(mu, s + e, 1) - f(mu, s - e, 1)) / (2 * e) - t[4]).max() < 1e-6
    assert np.abs((f(mu, s + e, 2) - f(mu, s - e, 2)) / (2 * e) - t[5]).max() < 1e-6


def test_uncensored_fit_equals_ols():
    rng = np.random.default_rng(1)
    n = 400
    g = rng.integers(0, 2, n).astype(float); c = rng.normal(size=n)
    y = 1.0 + 0.8 * g + 0.5 * c + rng.normal(0, 1, n)
    X = np.column_stack([np.ones(n), c, g])
    f = A.fit_censored(y, y, X)
    ols = sm.OLS(y, X).fit()
    assert f["status"] == "ok"
    assert np.allclose(f["beta"], ols.params, atol=1e-6)
    se = np.sqrt(np.diag(f["cov"])[:-1])
    assert np.allclose(se, ols.bse, rtol=0.03)                       # MLE vs df-corrected OLS SE
    assert np.exp(f["log_sigma"]) == pytest.approx(np.sqrt(ols.ssr / n), rel=1e-6)


def test_censored_model_beats_naive_linear_on_boundary_values():
    """The user-critical failure mode: naive OLS on '<=0.25' / '>4' treated as numbers attenuates the effect."""
    rng = np.random.default_rng(2)
    n = 3000
    g = rng.integers(0, 2, n).astype(float)
    latent = -1.0 + 3.0 * g + rng.normal(0, 1.2, n)
    lo = np.where(latent <= -2, -np.inf, np.where(latent >= 2, 2.0, np.round(latent)))
    hi = np.where(latent <= -2, -2.0, np.where(latent >= 2, np.inf, np.round(latent)))
    X = np.column_stack([np.ones(n), g])
    tobit = A.fit_censored(lo, hi, X)["beta"][1]
    y_naive = np.where(np.isfinite(lo) & np.isfinite(hi), lo, np.where(np.isinf(lo), hi, lo))
    naive = sm.OLS(y_naive, X).fit().params[1]
    assert abs(tobit - 3.0) < 0.25 and naive < 2.4 and abs(tobit - 3.0) < abs(naive - 3.0) / 2


def test_both_sided_censoring_rejected():
    with pytest.raises(ValueError):
        A.fit_censored(np.array([-np.inf, 0.0]), np.array([np.inf, 0.0]), np.ones((2, 1)))


def test_planted_determinants_recovered_with_correct_effect(assoc):
    r = assoc.set_index("feature_id")
    assert assoc.iloc[0]["feature_id"] == "gyrA_S83L"
    assert abs(r.loc["gyrA_S83L", "effect"] - 3.0) < 0.6
    assert abs(r.loc["parC_S80I", "effect"] - 1.6) < 0.8
    for f in ("gyrA_S83L", "gyrA_D87N", "parC_S80I"):
        assert r.loc[f, "q_value"] < 0.05 and r.loc[f, "effect"] > 0
    assert (r["status"] == "ok").all()


def test_effect_se_ci_p_from_same_model(assoc):
    assert A.check_amr_consistency(assoc) == []
    bad = assoc.copy(); bad.loc[0, "p_value"] = 0.5
    assert any("p inconsistent" in m for m in A.check_amr_consistency(bad))
    bad = assoc.copy(); bad.loc[1, "ci_high"] += 0.3
    assert any("CI inconsistent" in m for m in A.check_amr_consistency(bad))
    bad = assoc.copy(); bad.loc[2, "se"] *= 2
    assert A.check_amr_consistency(bad)
    # LRT is a diagnostic column that is NOT the reported p
    assert "p_lrt_diagnostic" in assoc and not np.allclose(assoc["p_value"], assoc["p_lrt_diagnostic"])


def test_bh_and_bonferroni_over_tested_features(assoc):
    from variantbridge.stats_utils import benjamini_hochberg
    assert np.allclose(assoc["q_value"], benjamini_hochberg(assoc["p_value"].to_numpy()))
    assert np.allclose(assoc["bonferroni_p"], np.minimum(1, assoc["p_value"] * len(assoc)))
    assert assoc.attrs["correction_method"] == "bh_fdr" and assoc.attrs["n_features_tested"] == len(assoc)


def test_lineage_adjustment_changes_confounded_signal_and_unadjusted_reported(assoc):
    assert assoc["lineage_adjusted"].all() and (assoc["structure_control"] == "lineage_fixed_effects").all()
    assert assoc["unadjusted_p"].notna().all() and not np.allclose(assoc["unadjusted_p"], assoc["p_value"])


def test_run_requires_registered_answer_key_hash(disc, cfg):
    for bad in ("", None, "abc"):
        with pytest.raises(ValueError, match="answer_key_sha256"):
            A.run_association(disc["G"][:, :5], disc["ids"][:5], "gene_pa", disc["pheno"], disc["L"], cfg, bad)


def test_discovery_only(syn, disc, cfg, answer_key):
    k = syn["kept"]
    for tier in ("replication", "exploratory"):
        sub = k[k["cohort_tier"] == tier].reset_index(drop=True)
        G = syn["features"].loc[sub["genome_id"]].to_numpy()[:, :5]
        with pytest.raises(cohorts.CohortSeparationError):
            A.run_association(G, disc["ids"][:5], "gene_pa", sub, None, cfg, answer_key[1])
    with pytest.raises(cohorts.CohortSeparationError):
        A.run_association(syn["features"].loc[k["genome_id"]].to_numpy()[:, :5], disc["ids"][:5], "gene_pa", k, None, cfg, answer_key[1])


def test_association_is_unbiased_no_annotation_inputs():
    """Addendum item 2/3: the association API has no way to receive annotation / answer-key / pathway information."""
    params = set(inspect.signature(A.run_association).parameters)
    assert not {p for p in params if any(w in p.lower() for w in ("annot", "pathway", "efflux", "database", "prior", "weight", "gene_set"))}
    mod_src = inspect.getsource(A)
    import re
    assert not re.search(r"^\s*(from\s+\.annotation|from\s+\.evidence|import\s+.*annotation)", mod_src, re.M), "association must not import the annotation/evidence layers"
    assert "label_answer_key" not in mod_src and "load_answer_key" not in mod_src and "answer_key.yaml" not in inspect.getsource(A.run_association)


def test_results_invariant_to_feature_naming_and_annotation(disc, cfg, answer_key, assoc):
    """Renaming features to hide any biological name must not change a single statistic."""
    anon = [f"F{i}" for i in range(len(disc["ids"]))]
    r2 = A.run_association(disc["G"], anon, "gene_pa", disc["pheno"], disc["L"], cfg, answer_key[1])
    a = assoc.set_index(pd.Index(assoc["feature_id"])).loc[disc["ids"]]
    b = r2.set_index("feature_id").loc[anon]
    assert np.allclose(a["effect"].to_numpy(), b["effect"].to_numpy()) and np.allclose(a["p_value"].to_numpy(), b["p_value"].to_numpy())


def test_frequency_filter_recorded_not_silent(disc, cfg, answer_key):
    G = disc["G"][:, :6].copy()
    G[:, 0] = 0; G[:3, 0] = 1                      # 3 isolates < min_isolates_with_feature (5)
    G[:, 1] = 1                                     # constant
    G[0, 2] = np.nan                                # missing -> excluded, never imputed
    r = A.run_association(G, disc["ids"][:6], "gene_pa", disc["pheno"], disc["L"], cfg, answer_key[1])
    reasons = dict(r.attrs["filtered"])
    assert reasons[disc["ids"][0]] == "frequency_filter" and reasons[disc["ids"][1]] == "constant" and reasons[disc["ids"][2]] == "missing_values"
    assert r.attrs["n_features_filtered"] == 3 and len(r) == 3


def test_replication_uses_same_model_and_flags_untestable(syn, disc, assoc, cfg):
    k = syn["kept"]; rp = k[k["cohort_tier"] == "replication"].reset_index(drop=True)
    G = syn["features"].loc[rp["genome_id"]].to_numpy()
    from variantbridge.amr import structure
    L, _ = structure.lineage_design(syn["st"].loc[rp["genome_id"]].reset_index(drop=True), cfg)
    out = A.replicate_features(assoc, G, disc["ids"], rp, L, cfg)
    top = out.set_index("feature_id")
    assert top.loc["gyrA_S83L", "same_direction"] and top.loc["gyrA_S83L", "replicated"]
    # drop one hit's column from the replication feature set => not_testable
    keep = [i for i, f in enumerate(disc["ids"]) if f != "gyrA_D87N"]
    out2 = A.replicate_features(assoc, G[:, keep], [disc["ids"][i] for i in keep], rp, L, cfg).set_index("feature_id")
    assert out2.loc["gyrA_D87N", "status"] == "not_testable"
    with pytest.raises(cohorts.CohortSeparationError):
        A.replicate_features(assoc, disc["G"], disc["ids"], disc["pheno"], disc["L"], cfg)      # discovery data is not replication


def test_binary_breakpoint_model_not_run_when_unverified(disc, cfg):
    with pytest.raises(BreakpointUnverified):
        A.run_binary_verified_breakpoint(disc["G"][:, :5], disc["ids"][:5], disc["pheno"], disc["L"], cfg)


def test_binary_breakpoint_model_runs_only_with_verified_test_config(disc, verified_cfg):
    r = A.run_binary_verified_breakpoint(disc["G"][:, 100:106], disc["ids"][100:106], disc["pheno"], disc["L"], verified_cfg)
    assert len(r) == 6 and (r["censoring_handling"] == "binary_verified_breakpoint").all()


# ------------------------------------------------------------------ pyseer wrapper (format fixture + fake runner only)
def test_pyseer_wrapper_not_run_without_binary(disc, cfg, tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", "")
    assert not A.pyseer_available()
    with pytest.raises(A.PyseerUnavailable):
        A.run_pyseer("gene_pa", "x.Rtab", disc["pheno"], cfg, tmp_path, "continuous_censoring_naive_SENSITIVITY_ONLY", acknowledge_censoring_naive=True)


def test_pyseer_modes_guard_censoring_and_breakpoints(disc, cfg, tmp_path):
    with pytest.raises(BreakpointUnverified):
        A.write_pyseer_phenotype(disc["pheno"], cfg, "binary_verified_breakpoint", tmp_path / "p.tsv")
    with pytest.raises(ValueError, match="acknowledge_censoring_naive"):
        A.write_pyseer_phenotype(disc["pheno"], cfg, "continuous_censoring_naive_SENSITIVITY_ONLY", tmp_path / "p.tsv")
    m = A.write_pyseer_phenotype(disc["pheno"], cfg, "continuous_censoring_naive_SENSITIVITY_ONLY", tmp_path / "p.tsv", acknowledge_censoring_naive=True)
    assert m["censoring_handling"] == "naive_boundary_substitution_SENSITIVITY_ONLY"
    with pytest.raises(cohorts.CohortSeparationError):
        A.write_pyseer_phenotype(disc["pheno"].assign(cohort_tier="replication"), cfg, "continuous_censoring_naive_SENSITIVITY_ONLY", tmp_path / "p.tsv", True)


def test_pyseer_binary_mode_with_test_breakpoints(disc, verified_cfg, tmp_path):
    m = A.write_pyseer_phenotype(disc["pheno"], verified_cfg, "binary_verified_breakpoint", tmp_path / "p.tsv")
    t = pd.read_csv(tmp_path / "p.tsv", sep="\t")
    assert set(t["phenotype"]) == {0, 1} and m["dropped"] > 0 and m["censoring_handling"] == "binary_verified_breakpoint"


def test_pyseer_parse_and_run_with_fake_runner(disc, cfg, tmp_path):
    fixture = pd.DataFrame({"variant": ["gene_a", "gene_b", "gene_c"], "af": [.3, .2, .1], "filter-pvalue": [1e-5, .2, .5], "lrt-pvalue": [1e-6, .3, .6],
                            "beta": [1.2, .1, -.2], "beta-std-err": [.2, .1, .3], "notes": [np.nan, np.nan, "bad-chisq"]})
    def fake(cmd, stdout, stderr, text):
        fixture.to_csv(stdout.name, sep="\t", index=False)
        class R: returncode = 0; stderr = ""
        assert cmd[0] == "pyseer" and "--pres" in cmd and "--continuous" in cmd
        return R()
    r = A.run_pyseer("gene_pa", "x.Rtab", disc["pheno"], cfg, tmp_path, "continuous_censoring_naive_SENSITIVITY_ONLY", acknowledge_censoring_naive=True, runner=fake)
    assert (r["censoring_handling"] == "naive_boundary_substitution_SENSITIVITY_ONLY").all() and r.iloc[0]["feature_id"] == "gene_a"
    assert r.set_index("feature_id").loc["gene_c", "status"].startswith("pyseer_flag") and np.isnan(r.set_index("feature_id").loc["gene_c", "q_value"])
    zc = stats.norm.ppf(0.975); a = r.set_index("feature_id").loc["gene_a"]
    assert a["p_value"] == pytest.approx(2 * stats.norm.sf(1.2 / .2)) and a["ci_low"] == pytest.approx(1.2 - zc * .2)       # p is Wald from beta/SE, consistent with CI
    with pytest.raises(KeyError, match="UNVERIFIED"):
        A.parse_pyseer_output(pd.DataFrame({"variant": ["a"]}), "gene_pa", 10, "x", "none")
