#!/usr/bin/env python
"""Evaluate the pre-registered answer-key criteria (A1, A3) from existing artifacts and record real-scope statuses.

    python scripts/amr_05_validate_answer_key.py --features data/work/amr_features.npz
Synthetic-scope inputs are evaluated and printed but written ONLY as synthetic (never shown as validation)."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

import _amr_common as K
from variantbridge.amr import validation as V
from variantbridge.amr.features_io import load_feature_matrix


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--answer-key", default=None)
    ap.add_argument("--features", required=True)
    ap.add_argument("--artifact-dir", default=str(K.DEFAULT_ART))
    a = ap.parse_args()
    cfg = K.C.load_config(a.config)
    key, sha = K.C.load_answer_key(a.answer_key)
    art = Path(a.artifact_dir)
    assoc, am = K.read_verified(art, "amr_association_discovery.csv.gz")
    if am.get("answer_key_sha256") != sha:
        sys.exit("answer_key.yaml changed after the association run (hash mismatch): re-run association; criteria must be fixed in advance")
    ph, _ = K.read_verified(art, "amr_phenotypes.csv.gz")
    X, samples, fids, fmeta, fm = load_feature_matrix(a.features)
    idx = {s: i for i, s in enumerate(samples)}
    d = ph[ph["cohort_tier"] == "discovery"]
    G = X[[idx[g] for g in d["genome_id"]]]
    cvb = cvr = None
    try:
        _, pm = K.read_verified(art, "amr_prediction_importance.csv.gz")
        res = pm["results"]
        cvb = {"mean_enet_extreme_auc": res["lineage_blocked"]["mean_enet_extreme_auc"]}
        cvr = {"mean_enet_extreme_auc": res["random_split_leakage_demo"]["mean_enet_extreme_auc"]}
    except SystemExit:
        print("prediction artifact absent: criterion C3 will be NOT RUN")
    out = V.validate_answer_key(assoc, key, cfg, fmeta, G, fids, cvb, cvr)
    scope = am["data_scope"]
    (art / "amr_answer_key_validation.json").write_text(json.dumps(dict(out, data_scope=scope, answer_key_sha256=sha, executed_at=V.now()), indent=2, default=str))
    print(f"answer-key validation: {out['status']} (data_scope={scope})")
    if scope == "real":
        V.record_status("A1", {"status": out["status"], "executed_at": V.now(), "data_scope": "real", "criteria": out["criteria"]}, art)
        if cvb is not None:
            c3 = out["criteria"]["C3_lineage_blocked_auc"]
            V.record_status("A3", {"status": V.PASSED if c3["passed"] else V.FAILED, "executed_at": V.now(), "data_scope": "real", "detail": c3}, art)
    else:
        print("SYNTHETIC scope: nothing recorded as validation.")


if __name__ == "__main__":
    main()
