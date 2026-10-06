#!/usr/bin/env python
"""Predictive baseline: lineage-blocked nested CV (primary within-discovery), random-split leakage demo (separate,
never a headline), discovery -> replication (primary external validation). Exploratory tier is not used here.

    python scripts/amr_04_run_prediction.py --phenotypes ... --features data/work/amr_features.npz
Targets are breakpoint-free: censored log2 MIC (ridge-Tobit, interval-aware MAE) and the extreme-class contrast (AUC)."""
import argparse
import sys
from pathlib import Path

import numpy as np

import _amr_common as K
from variantbridge.amr import prediction as PR, reporting
from variantbridge.amr.features_io import load_feature_matrix


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--phenotypes", default=str(K.DEFAULT_ART / "amr_phenotypes.csv.gz"))
    ap.add_argument("--features", required=True)
    ap.add_argument("--artifact-dir", default=str(K.DEFAULT_ART))
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    cfg = K.C.load_config(a.config)
    art = Path(a.artifact_dir)
    ph, ph_m = K.read_verified(Path(a.phenotypes).parent, Path(a.phenotypes).name)
    X, samples, fids, fmeta, f_m = load_feature_matrix(a.features)
    if f_m["data_scope"] != ph_m["data_scope"]:
        sys.exit("data_scope mismatch between phenotypes and features")
    idx = {s: i for i, s in enumerate(samples)}
    d = ph[(ph["cohort_tier"] == "discovery") & ph["genome_id"].isin(idx)].reset_index(drop=True)
    Xd = X[[idx[g] for g in d["genome_id"]]]
    st = d["mlst"].astype(object).where(d["mlst"].notna(), "unknown").astype(str).to_numpy()
    blocked = PR.nested_lineage_blocked_cv(Xd, d, st, fids, cfg)
    randdemo = PR.random_split_leakage_demo(Xd, d, st, fids, cfg, seed=a.seed)
    gap = PR.leakage_gap(blocked, randdemo)
    res = {"lineage_blocked": {k: blocked[k] for k in ("mean_enet_extreme_auc", "mean_tobit_extreme_auc", "mean_interval_mae", "cv_strategy")},
           "random_split_leakage_demo": {k: randdemo[k] for k in ("mean_enet_extreme_auc", "mean_tobit_extreme_auc", "mean_interval_mae")},
           "leakage_gap": gap, "discovery_to_replication": None}
    r = ph[(ph["cohort_tier"] == "replication") & ph["genome_id"].isin(idx)].reset_index(drop=True)
    if len(r):
        Xr = X[[idx[g] for g in r["genome_id"]]]
        res["discovery_to_replication"] = {k: v for k, v in PR.discovery_to_replication(Xd, d, st, fids, Xr, r, fids, cfg).items()
                                           if not k.startswith("_")}
    m = reporting.amr_manifest(cfg, phenotype_definition="breakpoint-free targets: censored log2 MIC; extreme-class contrast", inclusion_exclusion={"n_discovery": int(len(d))},
                               covariates="none (lineage held out by design)", model_specification={"models": "ridge-Tobit; elastic-net logistic", "cv": "nested lineage-blocked", "importance_warning": PR.IMPORTANCE_WARNING},
                               random_seed=a.seed, limitations=["extreme contrast is not an S/R call", "importance is predictive utility, not mechanism"],
                               data_scope=ph_m["data_scope"], extra={"results": res})
    K.write_table(blocked["importance"], art, "amr_prediction_importance.csv.gz", m)
    K.write_table(blocked["per_fold"], art, "amr_prediction_folds.csv.gz", m)
    print("lineage-blocked extreme AUC %.3f | random-split (LEAKAGE DEMO) %.3f | gap %.3f" % (gap["lineage_blocked_auc"], gap["random_split_auc"], gap["gap"]))


if __name__ == "__main__":
    main()
