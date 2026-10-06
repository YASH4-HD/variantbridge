"""Offline-script plumbing on SYNTHETIC inputs (scripts amr_00..06). Proves wiring, manifests, scope propagation and
guards - NOT biological validation."""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
S = ROOT / "scripts"


def run(script, *args, check=True):
    r = subprocess.run([sys.executable, str(S / script), *map(str, args)], capture_output=True, text=True, timeout=900, cwd=ROOT)
    if check:
        assert r.returncode == 0, f"{script}\n{r.stdout[-800:]}\n{r.stderr[-1500:]}"
    return r


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    d = tmp_path_factory.mktemp("amr_chain")
    run("amr_00_make_synthetic_inputs.py", "--out", d / "in")
    art = d / "art"
    run("amr_02_build_phenotypes.py", "--records-csv", d / "in/records.csv", "--meta-csv", d / "in/meta.csv", "--data-scope", "synthetic", "--artifact-dir", art)
    f = d / "in/features.npz"
    p = art / "amr_phenotypes.csv.gz"
    run("amr_03_run_association.py", "--phenotypes", p, "--features", f, "--artifact-dir", art)
    run("amr_04_run_prediction.py", "--phenotypes", p, "--features", f, "--artifact-dir", art)
    run("amr_05_validate_answer_key.py", "--features", f, "--artifact-dir", art)
    run("amr_06_build_evidence_cards.py", "--features", f, "--artifact-dir", art, "--top", 8)
    return d


def test_artifacts_have_manifests_and_synthetic_scope(chain):
    art = chain / "art"
    for name in ("amr_phenotypes", "amr_association_discovery", "amr_replication", "amr_prediction_importance"):
        m = json.loads((art / f"{name}.manifest.json").read_text())
        assert m["data_scope"] == "synthetic" and len(m["output_sha256"]) == 64 and m["breakpoint_status"]["discovery"] == "UNVERIFIED"
        assert m["qc_thresholds"] and all("provenance" in t for t in m["qc_thresholds"])
    m = json.loads((art / "amr_association_discovery.manifest.json").read_text())
    assert len(m["answer_key_sha256"]) == 64 and m["consistency_violations"] == []


def test_synthetic_run_records_no_validation_status(chain):
    from variantbridge.amr import validation as V
    s = V.load_status(chain / "art")
    assert all(v["status"] == "NOT RUN" for v in s.values())
    assert not (chain / "art" / "AMR_VALIDATION_STATUS.json").exists()


def test_phenotype_artifact_keeps_censoring_explicit(chain):
    ph = pd.read_csv(chain / "art/amr_phenotypes.csv.gz")
    for c in ("mic_raw", "mic_operator", "mic_boundary_mg_l", "censoring", "log2_mic_exact", "log2_boundary", "log2_lower", "log2_upper", "ordinal_code", "cohort_tier"):
        assert c in ph
    cens = ph[ph["censoring"] != "exact"]
    assert len(cens) > 0 and cens["log2_mic_exact"].isna().all() and cens["log2_boundary"].notna().all()
    assert set(ph["cohort_tier"]) == {"discovery", "replication", "exploratory"} and not ph["genome_id"].duplicated().any()
    assert "resistant_phenotype" in ph and ph["resistant_phenotype"].isna().all()       # reported calls are not used


def test_association_artifact_discovery_only_and_answer_key_not_in_columns(chain):
    a = pd.read_csv(chain / "art/amr_association_discovery.csv.gz")
    m = json.loads((chain / "art/amr_association_discovery.manifest.json").read_text())
    ph = pd.read_csv(chain / "art/amr_phenotypes.csv.gz")
    assert m["inclusion_exclusion"]["n_discovery"] == int((ph["cohort_tier"] == "discovery").sum())
    assert not {"in_answer_key", "answer_key_role"} & set(a.columns)           # labels are added only after the run
    assert a.iloc[0]["feature_id"] == "gyrA_S83L" and a["status"].eq("ok").all()


def test_evidence_cards_written_validated_and_synthetic(chain):
    from variantbridge.amr import evidence as E
    doc = json.loads((chain / "art/amr_evidence_cards.json").read_text())
    assert doc["data_scope"] == "synthetic" and doc["rejected"] == [] and len(doc["cards"]) >= 8
    for c in doc["cards"]:
        E.validate_card(c)
        assert c["claim_currently_justified"] == [E.SYNTHETIC_CLAIM] and c["causal_status"]["status"] == "not_established"


def test_answer_key_validation_json_scope(chain):
    v = json.loads((chain / "art/amr_answer_key_validation.json").read_text())
    assert v["data_scope"] == "synthetic" and v["status"] in ("PASSED", "FAILED", "INCOMPLETE") and len(v["answer_key_sha256"]) == 64


def test_scope_mismatch_refused(chain, tmp_path):
    from variantbridge.amr.features_io import load_feature_matrix, save_feature_matrix
    X, samples, fids, meta, m = load_feature_matrix(chain / "in/features.npz")
    save_feature_matrix(tmp_path / "real_labelled.npz", X, samples, fids, meta, {k: v for k, v in m.items() if k not in ("data_scope", "output_sha256", "n_samples", "n_features", "feature_meta_file", "feature_meta_sha256")}, data_scope="real")
    r = run("amr_03_run_association.py", "--phenotypes", chain / "art/amr_phenotypes.csv.gz", "--features", tmp_path / "real_labelled.npz",
            "--artifact-dir", tmp_path / "x_art", check=False)
    assert r.returncode != 0 and "data_scope mismatch" in (r.stdout + r.stderr)


def test_tampered_phenotype_artifact_refused(chain, tmp_path):
    import shutil
    art = tmp_path / "art"; shutil.copytree(chain / "art", art)
    with open(art / "amr_phenotypes.csv.gz", "ab") as fh:
        fh.write(b"x")
    r = run("amr_03_run_association.py", "--phenotypes", art / "amr_phenotypes.csv.gz", "--features", chain / "in/features.npz", "--artifact-dir", art, check=False)
    assert r.returncode != 0 and "checksum mismatch" in (r.stdout + r.stderr)


def test_answer_key_changed_after_association_is_refused(chain, tmp_path):
    import yaml
    from variantbridge.amr import config as C
    k = yaml.safe_load((C.CONFIG_DIR / "answer_key.yaml").read_text()); k["registered_on"] = "2026-10-06"
    kp = tmp_path / "changed_key.yaml"; kp.write_text(yaml.safe_dump(k))
    r = run("amr_05_validate_answer_key.py", "--answer-key", kp, "--features", chain / "in/features.npz", "--artifact-dir", chain / "art", check=False)
    assert r.returncode != 0 and "hash mismatch" in (r.stdout + r.stderr)


def test_synthetic_records_require_explicit_scope(chain):
    r = run("amr_02_build_phenotypes.py", "--records-csv", chain / "in/records.csv", "--meta-csv", chain / "in/meta.csv", "--artifact-dir", chain / "art2", check=False)
    assert r.returncode != 0 and "--data-scope is required" in (r.stdout + r.stderr)


def test_fetch_script_dry_run_and_blocked_network_is_not_run(tmp_path):
    r = run("amr_01_fetch_bvbrc.py", "--dry-run")
    assert "genome_amr" in r.stdout and "limit(" in r.stdout
    r2 = run("amr_01_fetch_bvbrc.py", "--raw-dir", tmp_path / "raw", check=False)
    assert r2.returncode == 2 and "NOT RUN" in r2.stderr
    r3 = run("amr_02_build_phenotypes.py", "--raw-dir", tmp_path / "raw", "--meta-csv", tmp_path / "none.csv", check=False)
    assert r3.returncode != 0 and "NOT RUN" in (r3.stdout + r3.stderr)


# ------------------------------------------------------------------ Streamlit viewer
def _app(monkeypatch, art=None):
    st_testing = pytest.importorskip("streamlit.testing.v1")
    if art is not None:
        monkeypatch.setenv("VARIANTBRIDGE_AMR_ARTIFACTS", str(art))
    else:
        monkeypatch.delenv("VARIANTBRIDGE_AMR_ARTIFACTS", raising=False)
    at = st_testing.AppTest.from_file(str(ROOT / "app.py"), default_timeout=300).run()
    at.radio[0].set_value("EvoResist-AI - Pathogen/AMR").run()
    at.radio[1].set_value("EvoResist v0.2 artifacts (offline results)").run()
    return at


def test_viewer_without_results_says_not_run_and_shows_only_a6(monkeypatch, tmp_path):
    at = _app(monkeypatch)       # the shipped repository artifacts: only the executed A6 literature check
    assert not at.exception, [e.value for e in at.exception]
    assert any("NOT RUN" in w.value for w in at.warning)


def test_viewer_with_synthetic_artifacts_labels_scope_and_separates_dimensions(monkeypatch, chain):
    at = _app(monkeypatch, chain / "art")
    assert not at.exception, [e.value for e in at.exception]
    assert any("SYNTHETIC" in e.value for e in at.error)
    md = " ".join(m.value for m in at.markdown)
    for h in ("1 · Association evidence", "2 · Predictive contribution", "3 · Biological plausibility", "4 · Causal status"):
        assert h in md
    assert "composite" not in md.lower()


def test_existing_cohortomics_and_demo_branches_still_render(monkeypatch):
    st_testing = pytest.importorskip("streamlit.testing.v1")
    at = st_testing.AppTest.from_file(str(ROOT / "app.py"), default_timeout=300).run()
    assert not at.exception, [e.value for e in at.exception]
    at.radio[0].set_value("CohortOmics - Human cohorts").run()
    at.radio[1].set_value("Validated artifacts (real data)").run()
    assert not at.exception, [e.value for e in at.exception]
