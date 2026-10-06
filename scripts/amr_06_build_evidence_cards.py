#!/usr/bin/env python
"""Build schema-validated Evidence Cards (v1.1) for the top statistical candidates + answer-key features.

    python scripts/amr_06_build_evidence_cards.py --top 25
Ranking is PURE association (BH q, |z|); annotation only labels. Cards that violate the claim restrictions are rejected."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import _amr_common as K
from variantbridge.amr import annotation as AN, evidence as E, reporting, validation as V
from variantbridge.amr.config import breakpoint_block


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--answer-key", default=None)
    ap.add_argument("--features", required=True, help="feature npz (for feature metadata)")
    ap.add_argument("--artifact-dir", default=str(K.DEFAULT_ART))
    ap.add_argument("--top", type=int, default=25)
    a = ap.parse_args()
    cfg = K.C.load_config(a.config)
    key, sha = K.C.load_answer_key(a.answer_key)
    art = Path(a.artifact_dir)
    assoc, am = K.read_verified(art, "amr_association_discovery.csv.gz")
    assert am["answer_key_sha256"] == sha, "answer key changed after association"
    from variantbridge.amr.features_io import load_feature_matrix
    _, _, _, fmeta, _ = load_feature_matrix(a.features)
    lab = AN.label_answer_key(assoc[assoc["status"] == "ok"], fmeta, key)
    ranked = E.rank_candidates(assoc, lab[["feature_id", "in_answer_key", "answer_key_role", "answer_key_match"]])
    pick = pd.concat([ranked.head(a.top), ranked[ranked["in_answer_key"] == True]]).drop_duplicates("feature_id")  # noqa: E712
    try:
        rep, _ = K.read_verified(art, "amr_replication.csv.gz"); rep = rep.set_index("feature_id")
    except SystemExit:
        rep = None
    pred = {}
    try:
        imp, pm = K.read_verified(art, "amr_prediction_importance.csv.gz"); imp = imp.set_index("feature_id")
        auc = pm["results"]["lineage_blocked"]["mean_enet_extreme_auc"]
    except SystemExit:
        imp, auc = None, None
    scope = am["data_scope"]
    ctx = {"analysis_date": V.now()[:10], "pipeline_version": "EvoResist-AI 0.2.0", "dataset_name": cfg["dataset"]["name"],
           "dataset_citation": cfg["dataset"]["citation"], "random_seed": am.get("random_seed"), "cohort_tier": "discovery",
           "censoring_handling": "interval_censored_gaussian", "breakpoint_status": breakpoint_block(cfg, "discovery")["status"],
           "answer_key_sha256": sha, "data_scope": scope, "censoring_summary": am["model_specification"]["censoring"], "alpha_fdr": 0.05}
    cards, rejected = [], []
    fm = fmeta.set_index("feature_id")
    for i, (_, h) in enumerate(pick.iterrows(), 1):
        fid = h["feature_id"]
        ann = {"role": h.get("answer_key_role") if isinstance(h.get("answer_key_role"), str) else None, "in_database": False, "database": None,
               "is_answer_key_determinant": bool(h.get("in_answer_key"))}
        p = None
        if imp is not None and fid in imp.index:
            p = {"cv_strategy": "nested_lineage_blocked", "validation_scope": "lineage_blocked_cv_within_discovery", "outer_fold_auc": auc,
                 "coefficient": float(imp.loc[fid, "mean_enet_coefficient"]), "rank": int(imp.loc[fid, "rank_in_model"])}
        rp = None
        if rep is not None and fid in rep.index:
            r = rep.loc[fid]
            rp = {"cohort_pmid": "34485958", "effect": None if pd.isna(r.get("effect")) else float(r["effect"]), "standard_error": None if pd.isna(r.get("se")) else float(r["se"]),
                  "p_value": None if pd.isna(r.get("p_value")) else float(r["p_value"]), "same_direction": None if pd.isna(r.get("same_direction")) else bool(r["same_direction"]),
                  "replicated": None if pd.isna(r.get("replicated")) else bool(r["replicated"]), "status": str(r["status"])}
        meta = fm.loc[fid].to_dict() if fid in fm.index else {}
        try:
            cards.append(E.build_card(h, f"EC-{i:04d}", ctx, ann, p, rp, None, {k: (None if pd.isna(v) else v) for k, v in meta.items()}))
        except E.EvidenceCardError as e:
            rejected.append({"feature_id": fid, "problems": e.problems})
    (art / "amr_evidence_cards.json").write_text(json.dumps({"data_scope": scope, "cards": cards, "rejected": rejected}, indent=2, default=str))
    print(f"{len(cards)} evidence cards written, {len(rejected)} rejected by validation (data_scope={scope})")


if __name__ == "__main__":
    main()
