import json

import numpy as np
import pandas as pd
import pytest

from variantbridge import artifacts, association, demo, evidence, io as vio, manifest as mf, prs, reporting, capability, stats_utils as su


def test_stats_utils():
    p = np.array([0.5, 0.01, 0.2, np.nan])
    q = su.benjamini_hochberg(p)
    assert np.isnan(q[3]) and q[1] == pytest.approx(0.03) and (q[:3] >= p[:3]).all()
    assert su.genomic_inflation_lambda(np.random.default_rng(0).uniform(size=20000)) == pytest.approx(1.0, abs=0.03)
    assert su.empirical_pvalue(5.0, [1, 2, 3]) == pytest.approx(1 / 4)
    assert su.empirical_pvalue(0.0, [1, 2, 3]) == 1.0
    e, o = su.qq_points(np.array([0.5, 0.1, 0.01]))
    assert len(e) == len(o) == 3 and (np.diff(o) <= 0).all() and (np.diff(e) <= 0).all()


def test_manifest_requires_all_fields_and_roundtrips(tmp_path):
    full = {k: "x" for k in mf.REQUIRED_FIELDS}
    full.pop("software_versions"); full.pop("analysis_date")
    m = mf.build_manifest(**full)
    assert "numpy" in m["software_versions"] and "T" in m["analysis_date"]
    incomplete = dict(full); incomplete.pop("citation")
    with pytest.raises(ValueError, match="citation"):
        mf.build_manifest(**incomplete)
    mf.write_manifest(tmp_path / "m.json", m)
    assert mf.read_manifest(tmp_path / "m.json")["dataset_name"] == "x" and mf.validate_manifest(m) == []


def test_artifact_loader_checks_manifest_and_checksum(tmp_path):
    assert artifacts.load_table_artifact("nope.csv.gz", tmp_path)[0] is None
    f = tmp_path / "t.csv.gz"; pd.DataFrame({"a": [1, 2]}).to_csv(f, index=False)
    df, m, msg = artifacts.load_table_artifact("t.csv.gz", tmp_path)
    assert df is None and "no manifest" in msg
    full = {k: "x" for k in mf.REQUIRED_FIELDS}
    mf.write_manifest(tmp_path / "t.manifest.json", mf.build_manifest(**{**full, "output_sha256": mf.file_sha256(f)}))
    df, m, msg = artifacts.load_table_artifact("t.csv.gz", tmp_path)
    assert msg == "" and len(df) == 2
    pd.DataFrame({"a": [9]}).to_csv(f, index=False)
    df, m, msg = artifacts.load_table_artifact("t.csv.gz", tmp_path)
    assert df is None and "checksum" in msg


def test_prs_leakage_guard_and_evaluation():
    Gd, yd, Gt, yt, beta = demo.simulate_discovery_target(seed=3)
    prs.assert_disjoint_cohorts(Gd.index, Gt.index)
    with pytest.raises(ValueError, match="leakage"):
        prs.assert_disjoint_cohorts(list(Gd.index) + ["T00001"], Gt.index)
    disc = association.scan(Gd.to_numpy(), yd.to_numpy(), None, list(Gd.columns))
    w = prs.select_weights(disc, 1e-4)
    score = prs.polygenic_score(Gt[w["variant"]].to_numpy(), w["effect"].to_numpy())
    ev = prs.evaluate_polygenic_score(score, yt.to_numpy(), n_permutations=200, n_bootstrap=100)
    assert ev["metric"] == "incremental R2" and ev["observed"] > 0.05 and ev["empirical_p_one_sided"] < 0.01
    assert "controlled-access" not in prs.DISEASE_PRS_LIMITATION and "independent" in prs.DISEASE_PRS_LIMITATION


def test_evidence_card_never_claims_causation():
    ge, ph, me = demo.human_demo()
    f = ge.drop(columns="sample_id").fillna(ge.drop(columns="sample_id").mean())
    res = association.scan(f.to_numpy(), ph["phenotype"].to_numpy(), None, list(f.columns))
    card = evidence.make_card(res.sort_values("p").iloc[0], "phenotype", "none", capability.SIMULATED)
    txt = evidence.render_text(card)
    assert "CLAIMS NOT JUSTIFIED" in txt and "causality" in txt and "SIMULATED" in txt
    assert "causes" not in card.claim_justified.lower() and "cause" not in card.claim_justified.lower().replace("because", "")
    assert card.ci_low < card.effect < card.ci_high


def test_load_csv_inputs_validation():
    g = pd.DataFrame({"sample_id": ["a", "b", "c"], "v1": [0, 1, np.nan], "v2": [2, 1, 0]})
    p = pd.DataFrame({"sample_id": ["a", "b", "c"], "phenotype": [0, 1, 1]})
    m, X, y, feats = vio.load_csv_inputs(g, p, "phenotype")
    assert X["v1"].isna().sum() == 1                      # NOT silently imputed
    with pytest.raises(ValueError, match="duplicate"):
        vio.load_csv_inputs(pd.concat([g, g.iloc[:1]]), p, "phenotype")
    with pytest.raises(ValueError, match="no overlapping"):
        vio.load_csv_inputs(g.assign(sample_id=["x", "y", "z"]), p, "phenotype")
    with pytest.raises(ValueError, match="not found"):
        vio.load_csv_inputs(g, p, "missing_col")


def test_genotype_block_roundtrip_and_validation(tmp_path):
    blk = vio.GenotypeBlock("1", np.array([1, 5, 9]), np.array(["a", "b", "c"]), np.array(["G"] * 3), np.array(["A"] * 3),
                            np.arange(6, dtype=np.float32).reshape(2, 3), np.array(["s1", "s2"]))
    blk.save(tmp_path / "b.npz")
    b2 = vio.GenotypeBlock.load(tmp_path / "b.npz")
    assert (b2.dosage == blk.dosage).all() and b2.chrom == "1"
    assert b2.subset_samples(["s2"]).dosage.tolist() == [[3.0, 4.0, 5.0]]
    with pytest.raises(KeyError):
        blk.subset_samples(["zzz"])
    with pytest.raises(ValueError):
        vio.GenotypeBlock("1", np.array([9, 1, 5]), np.array(["a", "b", "c"]), np.array(["G"] * 3), np.array(["A"] * 3),
                          np.zeros((2, 3), np.float32), np.array(["s1", "s2"]))


def test_report_shows_not_run_by_default(tmp_path):
    txt = reporting.methods_report(tmp_path)
    assert "NOT RUN" in txt and "PASSED" not in txt.replace("NOT RUN until executed", "")
    assert "Genetic Expression Score" in txt and "PRS" in txt   # PRS appears only as 'Disease-trait PRS: not possible'
    assert "Genetic Expression Score (eQTLGen -> GEUVADIS)" in txt


@pytest.mark.parametrize("script", sorted(__import__("pathlib").Path(__file__).resolve().parents[1].glob("scripts/0*.py")))
def test_every_script_help_runs(script):
    """--help must work for every offline script (catches argparse '%' help-string bugs)."""
    import subprocess, sys
    r = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr[-500:]
