"""Interpretable predictive baseline with lineage-blocked nested CV (audit step [9]).

Two breakpoint-FREE targets (no breakpoint is inferred):
  1. censored log2 MIC via L2-penalised interval-censored Gaussian regression ("ridge-Tobit");
     metric = interval-aware MAE (error is 0 when the prediction lies inside the observed interval).
  2. extreme-class contrast (phenotypes.extreme_contrast) via elastic-net logistic regression;
     metric = ROC-AUC. This is NOT a clinical S/R classifier.

Rules enforced:
  * the primary number is outer-fold, LINEAGE-BLOCKED (whole STs held out), feature screening and
    hyper-parameters chosen INSIDE training folds only;
  * random-split CV exists only as `random_split_leakage_demo` and is labelled leakage-only;
  * discovery -> replication is the primary external validation; exploratory-tier data are descriptive
    only (assert_tier rejects them here).
Feature importance is predictive utility, never mechanism."""
from __future__ import annotations

import warnings
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from .association import _kinds, _terms
from .cohorts import assert_tier
from .config import threshold
from .phenotypes import extreme_contrast

import sklearn

_SK_GE_18 = tuple(int(x) for x in sklearn.__version__.split(".")[:2]) >= (1, 8)
IMPORTANCE_WARNING = "Feature importance reflects predictive utility under this CV design, not causal contribution."


# ------------------------------------------------------------------------------------------ folds
def lineage_blocked_folds(groups: Sequence, n_folds: int) -> np.ndarray:
    """Deterministic fold id per sample; every lineage group lies entirely in ONE fold.
    Greedy bin-packing: largest groups first, each to the currently smallest fold."""
    g = pd.Series(list(groups)).astype(str)
    sizes = g.value_counts()
    if len(sizes) < n_folds:
        raise ValueError(f"need >= {n_folds} lineage groups for {n_folds} lineage-blocked folds, have {len(sizes)}")
    order = sorted(sizes.index, key=lambda k: (-sizes[k], k))
    load = np.zeros(n_folds, int)
    where = {}
    for k in order:
        f = int(np.argmin(load)); where[k] = f; load[f] += sizes[k]
    return g.map(where).to_numpy()


def random_folds(n: int, n_folds: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    f = np.arange(n) % n_folds
    rng.shuffle(f)
    return f


# ------------------------------------------------------------------------------------------ models
def fit_ridge_tobit(X: np.ndarray, lo: np.ndarray, hi: np.ndarray, l2: float) -> dict:
    n, p = X.shape
    kind = _kinds(lo, hi)
    A = np.column_stack([np.ones(n), X])
    y0 = np.where(kind == 0, lo, np.where(kind == -1, hi - 1.0, lo + 1.0))
    b0 = np.zeros(p + 1); b0[0] = y0.mean()
    s0 = np.log(max(y0.std(), 1e-2))
    th0 = np.concatenate([b0, [s0]])

    def f(th):
        b, s = th[:-1], th[-1]
        mu = A @ b
        ll, gm, gs, *_ = _terms(mu, s, lo, hi, kind)
        pen = 0.5 * l2 * float(b[1:] @ b[1:])
        g = np.concatenate([A.T @ gm, [gs.sum()]])
        g[1:-1] -= l2 * b[1:]
        return -(ll.sum()) + pen, -g

    r = minimize(f, th0, jac=True, method="L-BFGS-B", options={"maxiter": 500})
    return {"intercept": float(r.x[0]), "coef": r.x[1:-1].copy(), "log_sigma": float(r.x[-1]), "converged": bool(r.success)}


def predict_ridge_tobit(m: dict, X: np.ndarray) -> np.ndarray:
    return m["intercept"] + X @ m["coef"]


def interval_mae(pred: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> float:
    """Mean distance from the prediction to the observed interval (0 if inside). Censoring-aware."""
    d = np.where(pred < lo, lo - pred, np.where(pred > hi, pred - hi, 0.0))
    return float(np.mean(d))


def _screen(X: np.ndarray, code: np.ndarray, k: int) -> np.ndarray:
    """Top-k feature indices by |Spearman| with the ordinal MIC code (training data only)."""
    Xr = np.apply_along_axis(stats.rankdata, 0, X)
    cr = stats.rankdata(code)
    Xc = Xr - Xr.mean(0); cc = cr - cr.mean()
    den = np.sqrt((Xc ** 2).sum(0) * (cc ** 2).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(den > 0, (Xc * cc[:, None]).sum(0) / den, 0.0)
    return np.argsort(-np.abs(r), kind="stable")[: min(k, X.shape[1])]


def _fit_enet(X, y, C, l1):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        kw = {} if _SK_GE_18 else {"penalty": "elasticnet"}   # scikit-learn >=1.8 infers elastic-net from l1_ratio
        m = LogisticRegression(solver="saga", l1_ratio=l1, C=C, class_weight="balanced", max_iter=2000, random_state=0, **kw)
        return m.fit(X, y)


def _auc(y, s):
    return float(roc_auc_score(y, s)) if len(np.unique(y)) == 2 else float("nan")


def _fit_eval_fold(Xtr, ptr, Xte, pte, cfg, inner_groups, folds_inner):
    """Select hyper-parameters by inner CV on TRAIN only, refit on all TRAIN, evaluate on TEST."""
    top_k = threshold(cfg, "prediction.top_k_features")
    sel = _screen(Xtr, ptr["ordinal_code"].to_numpy(), top_k)
    A, B = Xtr[:, sel], Xte[:, sel]
    lo_tr, hi_tr = ptr["log2_lower"].to_numpy(float), ptr["log2_upper"].to_numpy(float)
    ex_tr = extreme_contrast(ptr).to_numpy()
    ex_te = extreme_contrast(pte).to_numpy()
    inner = (lineage_blocked_folds(inner_groups, folds_inner) if inner_groups is not None else random_folds(len(Xtr), folds_inner))
    # --- ridge-Tobit: choose l2 by inner interval-MAE
    best_l2, best_mae = None, np.inf
    for l2 in threshold(cfg, "prediction.ridge_tobit_l2_grid"):
        maes = []
        for f in range(folds_inner):
            tr, va = inner != f, inner == f
            if va.sum() == 0 or tr.sum() == 0:
                continue
            m = fit_ridge_tobit(A[tr], lo_tr[tr], hi_tr[tr], l2)
            maes.append(interval_mae(predict_ridge_tobit(m, A[va]), lo_tr[va], hi_tr[va]))
        if maes and np.mean(maes) < best_mae:
            best_l2, best_mae = l2, float(np.mean(maes))
    mt = fit_ridge_tobit(A, lo_tr, hi_tr, best_l2)
    mu = predict_ridge_tobit(mt, B)
    res = {"ridge_tobit_l2": best_l2, "interval_mae": interval_mae(mu, pte["log2_lower"].to_numpy(float), pte["log2_upper"].to_numpy(float)),
           "tobit_extreme_auc": _auc(ex_te[~np.isnan(ex_te)], mu[~np.isnan(ex_te)]) if (~np.isnan(ex_te)).any() else float("nan")}
    # --- elastic-net logistic on extreme contrast
    m_tr = ~np.isnan(ex_tr)
    res.update(enet_C=None, enet_extreme_auc=float("nan"), n_train_extreme=int(m_tr.sum()), n_test_extreme=int((~np.isnan(ex_te)).sum()))
    coefs = np.zeros(A.shape[1])
    if len(np.unique(ex_tr[m_tr])) == 2:
        best_C, best_a = None, -np.inf
        for C in threshold(cfg, "prediction.enet_C_grid"):
            aucs = []
            for f in range(folds_inner):
                tr, va = m_tr & (inner != f), m_tr & (inner == f)
                if len(np.unique(ex_tr[tr])) < 2 or len(np.unique(ex_tr[va])) < 2:
                    continue
                aucs.append(_auc(ex_tr[va], _fit_enet(A[tr], ex_tr[tr], C, threshold(cfg, "prediction.enet_l1_ratio")).decision_function(A[va])))
            if aucs and np.mean(aucs) > best_a:
                best_C, best_a = C, float(np.mean(aucs))
        if best_C is not None:
            me = _fit_enet(A[m_tr], ex_tr[m_tr], best_C, threshold(cfg, "prediction.enet_l1_ratio"))
            keep = ~np.isnan(ex_te)
            res["enet_C"] = best_C
            if keep.any():
                res["enet_extreme_auc"] = _auc(ex_te[keep], me.decision_function(B[keep]))
            coefs = me.coef_[0]
    res["_selected"] = sel
    res["_coef"] = coefs
    return res


def nested_lineage_blocked_cv(X: np.ndarray, pheno: pd.DataFrame, lineage: Sequence, feature_ids: Sequence[str], cfg: dict) -> dict:
    """Primary within-discovery predictive estimate (outer folds hold out whole STs)."""
    assert_tier(pheno, "discovery", "prediction.nested_lineage_blocked_cv")
    pheno = pheno.reset_index(drop=True)
    return _cv(np.asarray(X, float), pheno, np.asarray(lineage).astype(str), feature_ids, cfg, blocked=True)


def random_split_leakage_demo(X: np.ndarray, pheno: pd.DataFrame, lineage: Sequence, feature_ids: Sequence[str], cfg: dict, seed: int = 0) -> dict:
    """LEAKAGE DEMONSTRATION ONLY: random folds share lineages across train/test. Never a headline."""
    assert_tier(pheno, "discovery", "prediction.random_split_leakage_demo")
    out = _cv(np.asarray(X, float), pheno.reset_index(drop=True), np.asarray(lineage).astype(str), feature_ids, cfg, blocked=False, seed=seed)
    out["cv_strategy"] = "random_demo_leakage_only"
    return out


def _cv(X, pheno, lineage, feature_ids, cfg, blocked, seed=0):
    k_out, k_in = threshold(cfg, "prediction.outer_folds"), threshold(cfg, "prediction.inner_folds")
    outer = lineage_blocked_folds(lineage, k_out) if blocked else random_folds(len(X), k_out, seed)
    folds, sel_count = [], np.zeros(X.shape[1])
    coef_sum = np.zeros(X.shape[1])
    for f in range(k_out):
        te, tr = outer == f, outer != f
        if te.sum() == 0:
            continue
        ptr = pheno[tr].reset_index(drop=True)
        r = _fit_eval_fold(X[tr], ptr, X[te], pheno[te].reset_index(drop=True), cfg,
                           lineage[tr] if blocked else None, k_in)
        sel = r.pop("_selected"); coefs = r.pop("_coef")
        sel_count[sel] += 1
        if len(coefs) == len(sel):
            coef_sum[sel] += coefs
        r.update(fold=f, n_train=int(tr.sum()), n_test=int(te.sum()), n_test_lineages=int(len(set(lineage[te]))),
                 lineage_overlap_train_test=int(len(set(lineage[te]) & set(lineage[tr]))))
        folds.append(r)
    fd = pd.DataFrame(folds)
    imp = pd.DataFrame({"feature_id": list(feature_ids), "selected_in_n_folds": sel_count, "mean_enet_coefficient": coef_sum / max(1, len(fd))})
    imp = imp[imp["selected_in_n_folds"] > 0].sort_values(["selected_in_n_folds", "mean_enet_coefficient"], ascending=False).reset_index(drop=True)
    imp["rank_in_model"] = np.arange(1, len(imp) + 1)
    return {"cv_strategy": "nested_lineage_blocked" if blocked else "random_split", "per_fold": fd,
            "mean_enet_extreme_auc": float(np.nanmean(fd["enet_extreme_auc"])), "mean_tobit_extreme_auc": float(np.nanmean(fd["tobit_extreme_auc"])),
            "mean_interval_mae": float(np.nanmean(fd["interval_mae"])), "importance": imp, "importance_warning": IMPORTANCE_WARNING,
            "outer_folds": k_out, "inner_folds": k_in}


def leakage_gap(blocked: dict, random_demo: dict) -> dict:
    """Audit pass-criterion: report random-split minus lineage-blocked AUC as the leakage lesson."""
    return {"random_split_auc": random_demo["mean_enet_extreme_auc"], "lineage_blocked_auc": blocked["mean_enet_extreme_auc"],
            "gap": random_demo["mean_enet_extreme_auc"] - blocked["mean_enet_extreme_auc"],
            "note": "random-split value is a leakage demonstration only and must never be reported as performance"}


def discovery_to_replication(X_disc: np.ndarray, pheno_disc: pd.DataFrame, lineage_disc: Sequence, feats_disc: Sequence[str],
                             X_rep: np.ndarray, pheno_rep: pd.DataFrame, feats_rep: Sequence[str], cfg: dict) -> dict:
    """PRIMARY external validation: train (hyper-parameters by lineage-blocked inner CV) on discovery, test on replication.
    Only features present in both sets are used; the number dropped is reported. Extreme-contrast definitions are cohort-specific
    (discovery '>4' vs replication '>=4' are not identical classes) - reported as a comparability caveat."""
    assert_tier(pheno_disc, "discovery", "prediction.discovery_to_replication(train)")
    assert_tier(pheno_rep, "replication", "prediction.discovery_to_replication(test)")
    common = [f for f in feats_disc if f in set(feats_rep)]
    if not common:
        raise ValueError("no common features between discovery and replication")
    di = {f: i for i, f in enumerate(feats_disc)}; ri = {f: i for i, f in enumerate(feats_rep)}
    A = np.asarray(X_disc, float)[:, [di[f] for f in common]]
    B = np.asarray(X_rep, float)[:, [ri[f] for f in common]]
    r = _fit_eval_fold(A, pheno_disc.reset_index(drop=True), B, pheno_rep.reset_index(drop=True), cfg, np.asarray(lineage_disc).astype(str), threshold(cfg, "prediction.inner_folds"))
    r.pop("_selected"); r.pop("_coef")
    r.update(cv_strategy="discovery_to_replication", n_features_common=len(common), n_features_discovery_only=len(feats_disc) - len(common),
             caveat="extreme-contrast class definitions are cohort-specific; assay censoring bounds differ between cohorts (UNVERIFIED comparability)")
    return r


def describe_exploratory(pheno: pd.DataFrame) -> dict:
    """The ONLY thing exploratory-tier data may be used for here: description (no model metrics)."""
    from .phenotypes import censoring_summary
    return {"tier": "exploratory", "use": "descriptive only; never headline association or predictive metrics",
            "censoring": censoring_summary(pheno), "n": int(len(pheno))}
