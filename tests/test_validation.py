import json

import numpy as np
import pandas as pd
import pytest

from variantbridge import capability, validation as V


def test_v5_comparator_passes_identical_and_fails_perturbed():
    df = pd.DataFrame({"variant": ["a", "b", "c"], "effect": [0.1, -0.5, 2.0], "se": [0.05, 0.1, 0.4], "p": [0.04, 1e-6, 0.3]})
    assert V.v5_compare_to_plink2(df, df)["status"] == capability.PASSED
    bad = df.copy(); bad.loc[1, "effect"] *= 1.01
    assert V.v5_compare_to_plink2(df, bad)["status"] == capability.FAILED
    bad = df.copy(); bad.loc[0, "p"] = 0.2
    assert V.v5_compare_to_plink2(df, bad)["status"] == capability.FAILED
    assert V.v5_compare_to_plink2(df, df.iloc[:2])["status"] == capability.FAILED      # unmatched variant fails


def test_v3_null_calibration():
    assert V.v3_null_calibration([0.98, 1.02, 1.01])["status"] == capability.PASSED
    assert V.v3_null_calibration([1.4, 1.5])["status"] == capability.FAILED


def test_v2_sex_check():
    assert V.v2_sex_check({"n_compared": 100, "n_concordant": 100, "concordance": 1.0})["status"] == capability.PASSED
    assert V.v2_sex_check({"n_compared": 0, "n_concordant": 0, "concordance": float("nan")})["status"] == capability.FAILED


def synth_v4(same_dir=1.0, n=300, seed=0):
    r = np.random.default_rng(seed)
    genes = [f"ENSG{i:05d}" for i in range(n)]
    snps = [f"rs{i}" for i in range(n)]
    rho = r.choice([-1, 1], n) * r.uniform(0.3, 0.7, n)
    flip = r.random(n) > same_dir
    eff = np.where(flip, -np.sign(rho), np.sign(rho)) * r.uniform(0.2, 0.8, n)
    pv = 10 ** -r.uniform(3, 12, n)
    key = pd.DataFrame({"gene_id": genes, "snp_id": snps, "pvalue": pv, "rho": rho})
    mine_pairs = pd.DataFrame({"gene_id": genes, "variant": snps, "effect": eff, "p": pv * r.uniform(0.5, 2, n)})
    mine_best = pd.DataFrame({"gene_id": genes, "best_variant": snps, "effect": eff, "p_nominal": pv, "q_bh": 1e-3})
    return mine_best, mine_pairs, key


def test_v4_logic_on_synthetic_tables():
    best, pairs, key = synth_v4(1.0)
    ok = V.v4_answer_key_replication("EUR", best, pairs, key, key)
    assert ok["status"] == capability.PASSED and ok["same_direction_rate"] == 1.0 and ok["spearman_neglog10p"] > 0.9
    best, pairs, key = synth_v4(0.7)
    assert V.v4_answer_key_replication("EUR", best, pairs, key, key)["status"] == capability.FAILED
    best, pairs, key = synth_v4(1.0, n=40)
    assert V.v4_answer_key_replication("EUR", best, pairs, key, key)["status"] == capability.FAILED   # too few pairs for a verdict


def test_v4_flags_global_orientation_flip_and_does_not_autoflip():
    best, pairs, key = synth_v4(1.0)
    pairs = pairs.assign(effect=-pairs["effect"])
    res = V.v4_answer_key_replication("EUR", best, pairs, key, key)
    assert res["same_direction_rate"] == 0.0 and res["status"] == capability.FAILED
    assert any("orientation" in f for f in res["flags"])


def test_status_defaults_to_not_run_and_records_only_with_execution(tmp_path):
    st = capability.load_status(tmp_path)
    assert all(v["status"] == capability.NOT_RUN for v in st.values()) and set(st) == set(capability.VALIDATION_IDS)
    with pytest.raises(ValueError):
        capability.record_status("V4", {"status": capability.PASSED}, tmp_path)           # no executed_at -> refused
    with pytest.raises(KeyError):
        capability.record_status("V9", {"status": "PASSED", "executed_at": "x"}, tmp_path)
    capability.record_status("V3", {"status": capability.PASSED, "executed_at": V.now()}, tmp_path)
    st = capability.load_status(tmp_path)
    assert st["V3"]["status"] == capability.PASSED and st["V4"]["status"] == capability.NOT_RUN


def test_shipped_status_file_contains_no_results():
    import pathlib
    p = pathlib.Path(__file__).resolve().parents[1] / "artifacts" / "VALIDATION_STATUS.json"
    if p.exists():
        data = json.loads(p.read_text())
        # only entries that were genuinely executed may exist, and each must carry an execution timestamp
        for k, v in data.items():
            assert "executed_at" in v, f"{k} recorded without execution timestamp"
