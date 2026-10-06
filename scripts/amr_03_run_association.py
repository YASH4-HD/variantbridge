#!/usr/bin/env python
"""Genome-wide, unbiased discovery association (censored MIC model) + replication of pre-specified hits.

    python scripts/amr_03_run_association.py --phenotypes artifacts/amr/amr_phenotypes.csv.gz --features data/work/amr_features.npz

Order enforced: the pre-registered answer_key.yaml is hashed FIRST and the hash is stored in the output manifest.
The answer key is not read by the association model. Only the discovery tier is analysed; replication hits are
re-tested with the same model. Real-data validation statuses are recorded only when all inputs are data_scope=real."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import _amr_common as K
from variantbridge.amr import association as A, reporting, structure, validation as V
from variantbridge.amr.features_io import load_feature_matrix
from variantbridge.amr.phenotypes import censoring_summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--answer-key", default=None)
    ap.add_argument("--phenotypes", default=str(K.DEFAULT_ART / "amr_phenotypes.csv.gz"))
    ap.add_argument("--features", required=True)
    ap.add_argument("--artifact-dir", default=str(K.DEFAULT_ART))
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    cfg = K.C.load_config(a.config)
    key, key_sha = K.C.load_answer_key(a.answer_key)          # 1. register the key BEFORE any association
    print(f"answer_key.yaml sha256 = {key_sha}")
    art = Path(a.artifact_dir)
    ph_all, ph_m = K.read_verified(Path(a.phenotypes).parent, Path(a.phenotypes).name)
    X, samples, fids, fmeta, f_m = load_feature_matrix(a.features)
    scope = ph_m["data_scope"]
    if f_m["data_scope"] != scope:
        sys.exit(f"data_scope mismatch: phenotypes={scope} features={f_m['data_scope']}; refusing to mix")
    idx = {s: i for i, s in enumerate(samples)}
    ftypes = fmeta.set_index("feature_id").loc[fids, "feature_type"].tolist()
    # ---- discovery
    d = ph_all[ph_all["cohort_tier"] == "discovery"].reset_index(drop=True)
    miss = d.loc[~d["genome_id"].isin(idx), "genome_id"]
    if len(miss):
        sys.exit(f"{len(miss)} discovery genomes lack features (e.g. {miss.iloc[0]}): upstream feature calling incomplete")
    Gd = X[[idx[g] for g in d["genome_id"]]]
    L, lnames = structure.lineage_design(d["mlst"], cfg)
    res = A.run_association(Gd, fids, ftypes, d, L, cfg, key_sha)
    bad = A.check_amr_consistency(res)
    m = reporting.amr_manifest(
        cfg, phenotype_definition="censored log2 MIC (interval-censored Gaussian)", inclusion_exclusion={"n_discovery": int(len(d)), "features_filtered": res.attrs["n_features_filtered"], "features_tested": res.attrs["n_features_tested"]},
        covariates=lnames, model_specification={"model": A.MODEL_NAME, "multiple_testing": "BH-FDR primary; Bonferroni column", "censoring": censoring_summary(d)},
        random_seed=a.seed, limitations=["breakpoint UNVERIFIED: no S/R model run", "Gaussian assumption on log2 MIC", "linked features not separated"],
        data_scope=scope, extra={"answer_key_sha256": key_sha, "consistency_violations": bad, "input_phenotypes_sha256": ph_m["output_sha256"], "input_features_sha256": f_m["output_sha256"]})
    K.write_table(res, art, "amr_association_discovery.csv.gz", m)
    print(f"discovery: {len(d)} isolates, {len(res)} features tested, {int((res['q_value'] < 0.05).sum())} with q<0.05; consistency violations: {len(bad)}")
    if scope == "real":
        V.record_status("A2", {"status": V.PASSED if not bad else V.FAILED, "executed_at": V.now(), "data_scope": "real", "violations": bad[:20]}, art)
    # ---- replication of pre-specified hits (same model)
    r = ph_all[ph_all["cohort_tier"] == "replication"].reset_index(drop=True)
    if len(r):
        Gr = X[[idx[g] for g in r["genome_id"] if g in idx]]
        r = r[r["genome_id"].isin(idx)].reset_index(drop=True)
        Lr, rn = structure.lineage_design(r["mlst"], cfg)
        rep = A.replicate_features(res, Gr, fids, r, Lr, cfg)
        mr = dict(m, model_specification={"model": A.MODEL_NAME, "role": "replication of pre-specified discovery hits", "n_replication": int(len(r)),
                                          "bonferroni_threshold_over_hits": rep.attrs.get("bonferroni_threshold_over_hits")})
        K.write_table(rep, art, "amr_replication.csv.gz", mr)
        n_rep = int(rep["replicated"].sum()) if "replicated" in rep else 0
        print(f"replication: {len(rep)} discovery hits tested, {n_rep} replicated")
        if scope == "real":
            conc = float(rep["same_direction"].mean()) if "same_direction" in rep and len(rep) else float("nan")
            V.record_status("A4", {"status": V.PASSED if (len(rep) and n_rep > 0) else V.FAILED, "executed_at": V.now(), "data_scope": "real",
                                   "n_hits_tested": int(len(rep)), "n_replicated": n_rep, "direction_concordance": conc,
                                   "note": "reported as observed; a failed replication is a result, not a pipeline error"}, art)


if __name__ == "__main__":
    main()
