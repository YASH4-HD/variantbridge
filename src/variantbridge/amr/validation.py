"""Pre-registered known-answer validation for EvoResist-AI (audit section 6) and status registry.

`validate_answer_key` reads the pre-registered key + config thresholds. It never changes a criterion
after seeing results. A software test passing on SYNTHETIC data is not a validation: statuses are
recorded with data_scope and only data_scope == 'real' records are shown as validation in the app."""
from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .annotation import label_answer_key
from .config import threshold

NOT_RUN, PASSED, FAILED, INCOMPLETE, UNVERIFIED = "NOT RUN", "PASSED", "FAILED", "INCOMPLETE", "UNVERIFIED"

AMR_VALIDATION_IDS = {
    "A1": "Answer-key recovery: gyrA/parC QRDR among top FDR-significant lineage-adjusted hits, positive direction",
    "A2": "Association output consistency (effect/SE/CI/p from one censored model) on the real discovery run",
    "A3": "Predictive baseline under lineage-blocked nested CV (extreme-contrast AUC >= pre-registered minimum) + leakage gap reported",
    "A4": "Discovery -> replication: direction and effect presence of discovery hits",
    "A5": "Real-data phenotype reconciliation with addendum 1.1 counts (informational)",
    "A6": "Discovery-cohort breakpoint/AST standard verified from source (gates any S/R analysis)",
}


def now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def status_file(root: str | Path = "artifacts/amr") -> Path:
    return Path(root) / "AMR_VALIDATION_STATUS.json"


def record_status(vid: str, record: dict, root: str | Path = "artifacts/amr") -> Path:
    if vid not in AMR_VALIDATION_IDS:
        raise KeyError(vid)
    for k in ("status", "executed_at", "data_scope"):
        if k not in record:
            raise ValueError(f"record must include {k!r}")
    if record["data_scope"] not in ("real", "synthetic"):
        raise ValueError("data_scope must be 'real' or 'synthetic'")
    p = status_file(root)
    p.parent.mkdir(parents=True, exist_ok=True)
    data = json.loads(p.read_text()) if p.exists() else {}
    data[vid] = record
    p.write_text(json.dumps(data, indent=2, default=str))
    return p


def load_status(root: str | Path = "artifacts/amr") -> dict:
    """Missing id => NOT RUN. A synthetic-scope record is NOT a validation and reads as NOT RUN (with a note)."""
    p = status_file(root)
    data = json.loads(p.read_text()) if p.exists() else {}
    out = {}
    for vid in AMR_VALIDATION_IDS:
        r = data.get(vid)
        if r is None:
            out[vid] = {"status": NOT_RUN}
        elif r.get("data_scope") != "real":
            out[vid] = {"status": NOT_RUN, "note": "only a SYNTHETIC-scope run is recorded; it is not a validation"}
        else:
            out[vid] = r
    return out


def linkage_clusters(G: np.ndarray, ids: list[str], r2_min: float) -> dict[str, int]:
    """Greedy clustering of features by squared Pearson correlation (input order = priority order)."""
    G = np.asarray(G, float)
    Z = G - G.mean(0)
    sd = Z.std(0)
    sd[sd == 0] = np.inf
    Z = Z / sd
    n = G.shape[0]
    cl, nxt = {}, 0
    for i, f in enumerate(ids):
        if f in cl:
            continue
        cl[f] = nxt
        r2 = (Z[:, i] @ Z[:, i + 1:] / n) ** 2 if i + 1 < len(ids) else np.array([])
        for jj, v in enumerate(r2):
            g = ids[i + 1 + jj]
            if g not in cl and v >= r2_min:
                cl[g] = nxt
        nxt += 1
    return cl


def validate_answer_key(assoc: pd.DataFrame, key: dict, cfg: dict, feature_meta: Optional[pd.DataFrame] = None,
                        G: Optional[np.ndarray] = None, feature_ids: Optional[list[str]] = None,
                        cv_blocked: Optional[dict] = None, cv_random: Optional[dict] = None) -> dict:
    """Evaluate the audit section-6 pass criteria with pre-registered thresholds. Returns a result dict; does not write files."""
    top_n = threshold(cfg, "answer_key.top_n_hits")
    r2 = threshold(cfg, "answer_key.linkage_cluster_r2")
    auc_min = threshold(cfg, "prediction.answer_key_auc_min")
    alpha = threshold(cfg, "association.alpha_fdr")
    lab = label_answer_key(assoc[assoc["status"] == "ok"], feature_meta, key)
    sig = lab[lab["q_value"] < alpha].sort_values("p_value").reset_index(drop=True)
    # clusters among significant features (so perfectly linked features count once)
    if G is not None and feature_ids is not None and len(sig):
        idx = {f: i for i, f in enumerate(feature_ids)}
        order = [f for f in sig["feature_id"] if f in idx]
        cl = linkage_clusters(G[:, [idx[f] for f in order]], order, r2)
    else:
        cl = {f: i for i, f in enumerate(sig["feature_id"])}
    sig["cluster"] = sig["feature_id"].map(cl)
    best = sig.drop_duplicates("cluster").reset_index(drop=True)
    best["cluster_rank"] = np.arange(1, len(best) + 1)
    prim = sig[sig["answer_key_role"] == "primary_QRDR"]
    prim_clusters = set(prim["cluster"])
    in_top = best[best["cluster"].isin(prim_clusters) & (best["cluster_rank"] <= top_n)]
    c1 = {"passed": bool(len(in_top) > 0), "n_significant_features": int(len(sig)), "n_significant_clusters": int(len(best)),
          "qrdr_features_significant": prim["feature_id"].tolist(), "qrdr_best_cluster_rank": int(best[best["cluster"].isin(prim_clusters)]["cluster_rank"].min()) if len(prim) else None,
          "top_n": top_n, "criterion": "at least one gyrA/parC QRDR cluster is FDR-significant and within the top-n clusters (lineage-adjusted model)"}
    c2 = {"passed": bool(len(prim) > 0 and (prim["effect"] > 0).all()), "criterion": "recovered QRDR features have positive effect",
          "directions": {r.feature_id: float(r.effect) for r in prim.itertuples()}}
    sec = lab[lab["answer_key_role"].isin(["secondary_PMQR", "secondary_efflux_regulator"])]
    secondary = [{"feature_id": r.feature_id, "role": r.answer_key_role, "effect": float(r.effect), "q_value": float(r.q_value),
                  "significant": bool(r.q_value < alpha)} for r in sec.itertuples()]
    if cv_blocked is None:
        c3 = {"passed": None, "status": NOT_RUN, "criterion": f"lineage-blocked extreme-contrast AUC >= {auc_min}"}
    else:
        auc = cv_blocked["mean_enet_extreme_auc"]
        c3 = {"passed": bool(np.isfinite(auc) and auc >= auc_min), "auc_lineage_blocked": float(auc), "auc_min": auc_min,
              "criterion": f"lineage-blocked extreme-contrast AUC >= {auc_min}",
              "leakage_gap": (None if cv_random is None else float(cv_random["mean_enet_extreme_auc"] - auc))}
    crit = {"C1_qrdr_recovered_top_hits": c1, "C2_direction_positive": c2, "C3_lineage_blocked_auc": c3}
    if c3["passed"] is None:
        overall = FAILED if not (c1["passed"] and c2["passed"]) else INCOMPLETE
    else:
        overall = PASSED if (c1["passed"] and c2["passed"] and c3["passed"]) else FAILED
    return {"status": overall, "criteria": crit, "secondary_reported_not_scored": secondary,
            "failure_handling": None if overall == PASSED else key.get("failure_handling")}
