import json

import numpy as np
import pandas as pd
import pytest

from variantbridge.amr import data_io as D, validation as V
from variantbridge.amr import association as A


# ------------------------------------------------------------------ status registry
def test_missing_status_is_not_run(tmp_path):
    s = V.load_status(tmp_path)
    assert set(s) == set(V.AMR_VALIDATION_IDS) and all(v["status"] == "NOT RUN" for v in s.values())


def test_record_requires_fields_and_scope(tmp_path):
    for bad in ({"status": "PASSED"}, {"status": "PASSED", "executed_at": "x"}, {"status": "PASSED", "executed_at": "x", "data_scope": "maybe"}):
        with pytest.raises(ValueError):
            V.record_status("A1", bad, tmp_path)
    with pytest.raises(KeyError):
        V.record_status("Z9", {"status": "PASSED", "executed_at": "x", "data_scope": "real"}, tmp_path)


def test_synthetic_record_never_reads_as_validation(tmp_path):
    V.record_status("A1", {"status": "PASSED", "executed_at": V.now(), "data_scope": "synthetic"}, tmp_path)
    r = V.load_status(tmp_path)["A1"]
    assert r["status"] == "NOT RUN" and "SYNTHETIC" in r["note"]
    V.record_status("A1", {"status": "PASSED", "executed_at": V.now(), "data_scope": "real"}, tmp_path)
    assert V.load_status(tmp_path)["A1"]["status"] == "PASSED"


def test_shipped_status_file_is_honest():
    """The repository ships ONLY the executed literature check (A6 UNVERIFIED); nothing is pre-filled as passed."""
    from pathlib import Path
    root = Path(__file__).resolve().parents[2] / "artifacts" / "amr"
    s = V.load_status(root)
    assert s["A6"]["status"] == "UNVERIFIED" and s["A6"]["data_scope"] == "real"
    assert all(s[k]["status"] == "NOT RUN" for k in ("A1", "A2", "A3", "A4", "A5"))
    assert not any(v["status"] == "PASSED" for v in s.values())


# ------------------------------------------------------------------ answer-key validation logic
def _assoc(rows):
    d = pd.DataFrame(rows, columns=["feature_id", "effect", "p_value"])
    d["status"] = "ok"; d["se"] = 0.1; d["z"] = d["effect"] / 0.1
    from variantbridge.stats_utils import benjamini_hochberg
    d["q_value"] = benjamini_hochberg(d["p_value"].to_numpy())
    return d


def test_validation_passes_fails_and_incomplete(cfg, answer_key):
    key = answer_key[0]
    good = _assoc([("gyrA_S83L", 3.0, 1e-30), ("parC_S80I", 1.5, 1e-8)] + [(f"g{i}", 0.01, 0.5 + i * 1e-3) for i in range(60)])
    inc = V.validate_answer_key(good, key, cfg)
    assert inc["status"] == V.INCOMPLETE and inc["criteria"]["C3_lineage_blocked_auc"]["status"] == "NOT RUN"
    ok = V.validate_answer_key(good, key, cfg, cv_blocked={"mean_enet_extreme_auc": 0.9}, cv_random={"mean_enet_extreme_auc": 0.99})
    assert ok["status"] == V.PASSED and ok["criteria"]["C3_lineage_blocked_auc"]["leakage_gap"] == pytest.approx(0.09)
    low = V.validate_answer_key(good, key, cfg, cv_blocked={"mean_enet_extreme_auc": 0.6})
    assert low["status"] == V.FAILED and low["failure_handling"]
    wrong_dir = _assoc([("gyrA_S83L", -3.0, 1e-30)] + [(f"g{i}", 0.01, 0.5) for i in range(60)])
    assert V.validate_answer_key(wrong_dir, key, cfg, cv_blocked={"mean_enet_extreme_auc": 0.9})["status"] == V.FAILED
    missed = _assoc([(f"g{i}", 0.5, 1e-6 * (i + 1)) for i in range(25)] + [("gyrA_S83L", 3.0, 1e-2)] + [(f"h{i}", 0, 0.9) for i in range(30)])
    r = V.validate_answer_key(missed, key, cfg, cv_blocked={"mean_enet_extreme_auc": 0.9})
    assert r["criteria"]["C1_qrdr_recovered_top_hits"]["passed"] is False and r["status"] == V.FAILED


def test_secondary_determinants_reported_not_scored(cfg, answer_key):
    a = _assoc([("gyrA_S83L", 3.0, 1e-30), ("acrR_Q15*", 0.3, 0.6), ("qnrS1", 1.0, 1e-5)] + [(f"g{i}", 0.0, 0.7) for i in range(40)])
    r = V.validate_answer_key(a, answer_key[0], cfg, cv_blocked={"mean_enet_extreme_auc": 0.9})
    sec = {s["feature_id"]: s for s in r["secondary_reported_not_scored"]}
    assert sec["acrR_Q15*"]["significant"] is False and sec["qnrS1"]["significant"] is True     # a null is reported, never hidden
    assert r["status"] == V.PASSED


def test_linkage_clusters_count_perfectly_linked_once(cfg):
    rng = np.random.default_rng(0)
    a = rng.integers(0, 2, 200).astype(float)
    G = np.column_stack([a, a, 1 - a, rng.integers(0, 2, 200)])
    cl = V.linkage_clusters(G, ["x", "y", "z", "w"], 0.95)
    assert cl["x"] == cl["y"] == cl["z"] != cl["w"]


def test_synthetic_end_to_end_answer_key_recovery(assoc, syn, disc, cfg, answer_key):
    """SYNTHETIC: shows the software recovers PLANTED determinants. Not biological validation."""
    r = V.validate_answer_key(assoc, answer_key[0], cfg, syn["feature_meta"], disc["G"], disc["ids"], {"mean_enet_extreme_auc": 0.95})
    assert r["criteria"]["C1_qrdr_recovered_top_hits"]["passed"] and r["criteria"]["C2_direction_positive"]["passed"]


# ------------------------------------------------------------------ data_io (no network)
def _page(n, start=0):
    return json.dumps([{"genome_id": f"562.{i}", "measurement": "<=0.25", "pmid": "38052776", "antibiotic": "ciprofloxacin"} for i in range(start, start + n)]).encode()


def test_no_cache_no_network_is_not_run(cfg, tmp_path):
    with pytest.raises(D.DataNotAvailable):
        D.fetch_genome_amr(cfg, tmp_path)
    with pytest.raises(D.DataNotAvailable):
        D.fetch_genome_amr(cfg, tmp_path, http_get=lambda u: _page(1), allow_network=False)


def test_paged_fetch_is_resumable_and_logged(cfg, tmp_path, monkeypatch):
    monkeypatch.setattr(D, "PAGE", 3)
    calls = []
    def get(url):
        calls.append(url); return _page(3 if len(calls) == 1 else 2, 0 if len(calls) == 1 else 3)
    pages = D.fetch_genome_amr(cfg, tmp_path, get, allow_network=True)
    assert len(pages) == 2 and len(calls) == 2 and "limit(3,0)" in calls[0] and "limit(3,3)" in calls[1]
    again = D.fetch_genome_amr(cfg, tmp_path, get, allow_network=True)
    assert len(again) == 2 and len(calls) == 2                       # nothing re-downloaded
    log = json.loads((tmp_path / "DOWNLOAD_LOG.json").read_text())
    assert set(log) == {p.name for p in pages} and all(len(v["sha256"]) == 64 for v in log.values())
    rec = D.load_cached_records(pages); assert len(rec) == 5


def test_field_alias_failure_is_loud(tmp_path):
    p = tmp_path / "genome_amr_page_0000.json"; p.write_text(json.dumps([{"foo": 1}]))
    with pytest.raises(KeyError, match="UNVERIFIED"):
        D.load_cached_records([p])


def test_fetch_phenotypes_adds_tier_and_censoring(cfg):
    rec = pd.DataFrame({"genome_id": ["a", "b", "c"], "pmid": ["38052776", "34485958", "999"], "measurement": ["<=0.25", ">=4", "1"], "laboratory_typing_method": "MIC"})
    out = D.fetch_phenotypes(cfg, records=rec)
    assert out["cohort_tier"].tolist() == ["discovery", "replication", "exploratory"] and out["censoring"].tolist() == ["left", "right", "exact"]


def test_metadata_fetch_resumable(tmp_path):
    n = []
    def get(url):
        n.append(url); return json.dumps([{"genome_id": "a", "mlst": "ST10"}]).encode()
    D.fetch_genome_metadata(["a", "b"], tmp_path, get, True, batch=1)
    assert len(n) == 2
    D.fetch_genome_metadata(["a", "b"], tmp_path, get, True, batch=1)
    assert len(n) == 2
    with pytest.raises(D.DataNotAvailable):
        D.fetch_genome_metadata(["a", "z"], tmp_path / "other", None, False)


def test_contigs_to_fasta_and_validation():
    fa = D.contigs_json_to_fasta([{"sequence_id": "c1", "sequence": "acgt" * 30}], "562.1")
    assert fa.startswith(">562.1|c1\n") and max(len(l) for l in fa.splitlines()[1:]) <= 80 and fa.replace("\n", "").endswith("ACGT")
    with pytest.raises(ValueError):
        D.contigs_json_to_fasta([{"sequence": "ACGJ"}], "g")
