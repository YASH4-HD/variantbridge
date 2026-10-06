"""Known-answer validation tests V1-V6 (audit section 6, with auditor amendments).

Functions here are pure computations returning result dicts. They never write files and never
decide that a validation "passed" on data that was not run. The offline scripts call them on
real data and then persist the outcome with ``capability.record_status``.

Pass/fail criteria are PRE-REGISTERED here (values below, each labelled [project design choice]
or [source]) BEFORE any real-data run. Changing a criterion after seeing results requires an
explicit CHANGELOG entry; do not tune them to make a validation pass.
"""
from __future__ import annotations

import datetime as _dt
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, silhouette_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neighbors import KNeighborsClassifier

from .capability import FAILED, PASSED
from .stats_utils import genomic_inflation_lambda

# ---- pre-registered criteria ------------------------------------------------------------
V1_MIN_BALANCED_ACCURACY = 0.90     # [project design choice] kNN CV balanced accuracy, leading PCs vs superpopulation
V1_N_PCS = 5                        # [project design choice]
V1_KNN_K = 15                       # [project design choice]
V1_CV_FOLDS = 5                     # [project design choice]
V3_LAMBDA_TOLERANCE = 0.10          # [project design choice] |median lambda - 1| <= 0.10 on permuted phenotypes
V4_MIN_SAME_DIRECTION = 0.90        # [project design choice stated in audit V4 as starting point]
V4_MIN_PAIRS = 100                  # [project design choice] minimum compared pairs for a meaningful verdict
V5_BETA_SE_RTOL = 1e-4              # [project design choice] PLINK2 prints ~6 significant digits
V5_NEGLOG10P_ATOL = 5e-3            # [project design choice] absolute difference in -log10(p) (for p > 1e-30)


def now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


# ------------------------------------------------------------------------------------ V1
def v1_population_structure(pcs: pd.DataFrame, labels: pd.Series, n_pcs: int = V1_N_PCS, k: int = V1_KNN_K,
                            cv: int = V1_CV_FOLDS, seed: int = 0,
                            min_balanced_accuracy: float = V1_MIN_BALANCED_ACCURACY) -> dict:
    """Do the leading PCs recover the established broad (super)population structure?

    Replaces the rigid "PC1 separates superpopulations" requirement (amendment 2). Quantified by:
    (a) per-PC variance explained by the population label (eta^2),
    (b) cross-validated kNN balanced accuracy of the label from the first ``n_pcs`` PCs,
    (c) silhouette of the labels in PC space, (d) confusion matrix.
    """
    d = pcs.set_index("sample_id")
    lab = labels.reindex(d.index)
    if lab.isna().any():
        raise ValueError("labels missing for some samples")
    cols = [f"PC{i + 1}" for i in range(min(n_pcs, sum(c.startswith("PC") for c in d.columns)))]
    X = d[cols].to_numpy()
    y = lab.astype(str).to_numpy()
    eta = {}
    for c in cols:
        v = d[c].to_numpy()
        grand = v.mean()
        ss_between = sum((y == g).sum() * (v[y == g].mean() - grand) ** 2 for g in np.unique(y))
        ss_total = ((v - grand) ** 2).sum()
        eta[c] = float(ss_between / ss_total) if ss_total > 0 else float("nan")
    classes = sorted(np.unique(y))
    min_class = min((y == c).sum() for c in classes)
    folds = int(min(cv, min_class))
    pred = cross_val_predict(KNeighborsClassifier(n_neighbors=min(k, max(1, min_class - 1))), X, y,
                             cv=StratifiedKFold(folds, shuffle=True, random_state=seed))
    bacc = float(balanced_accuracy_score(y, pred))
    sil = float(silhouette_score(X, y)) if len(classes) > 1 else float("nan")
    return {
        "test": "V1", "n_samples": int(len(y)), "classes": classes, "n_pcs_used": len(cols),
        "eta_squared_by_pc": eta, "knn_cv_balanced_accuracy": bacc, "silhouette": sil,
        "confusion_matrix": confusion_matrix(y, pred, labels=classes).tolist(),
        "criterion": f"kNN CV balanced accuracy >= {min_balanced_accuracy} using first {len(cols)} PCs",
        "status": PASSED if bacc >= min_balanced_accuracy else FAILED, "executed_at": now(),
    }


# ------------------------------------------------------------------------------------ V2
def v2_sex_check(concordance: dict, min_concordance: float = 0.99) -> dict:
    ok = concordance["n_compared"] > 0 and concordance["concordance"] >= min_concordance
    return {"test": "V2", **concordance, "criterion": f"concordance >= {min_concordance} [project design choice]",
            "status": PASSED if ok else FAILED, "executed_at": now()}


# ------------------------------------------------------------------------------------ V3
def v3_null_calibration(lambdas: list[float], tolerance: float = V3_LAMBDA_TOLERANCE) -> dict:
    lam = np.asarray(lambdas, dtype=float)
    med = float(np.nanmedian(lam))
    return {"test": "V3", "lambdas": lam.tolist(), "median_lambda": med,
            "criterion": f"|median lambda on permuted phenotypes - 1| <= {tolerance}",
            "status": PASSED if abs(med - 1) <= tolerance else FAILED, "executed_at": now(),
            "note": "Lambda on the REAL cis scan is expected to exceed 1 (true effects + LD); only the permuted-phenotype lambda is a calibration check."}


# ------------------------------------------------------------------------------------ V4
def v4_answer_key_replication(population: str, mine_best: pd.DataFrame, mine_pairs: pd.DataFrame,
                              key_best: pd.DataFrame, key_all: pd.DataFrame,
                              fdr: float = 0.05, min_same_direction: float = V4_MIN_SAME_DIRECTION,
                              min_pairs: int = V4_MIN_PAIRS) -> dict:
    """Compare re-derived cis-eQTL results with the published GEUVADIS answer key.

    key_best / key_all columns (normalised by the script): gene_id, snp_id, pvalue, rho.
    mine_best: gene_id, best_variant, effect, p_nominal, q_bh.  mine_pairs: gene_id, variant, effect, p.
    The orientation of the published rho (which allele is counted) is NOT assumed: a near-zero
    same-direction rate is reported as a possible global orientation difference and is never
    silently flipped.
    """
    key_best = key_best.drop_duplicates(["gene_id", "snp_id"])
    mp = mine_pairs.rename(columns={"variant": "snp_id"}).drop_duplicates(["gene_id", "snp_id"])
    # same direction among published best SNPs retained in my scan
    j = key_best.merge(mp, on=["gene_id", "snp_id"], how="inner", suffixes=("_key", "_mine"))
    n_dir = int(len(j))
    same = float((np.sign(j["rho"]) == np.sign(j["effect"])).mean()) if n_dir else float("nan")
    # rank correlation of -log10 p over shared gene-SNP pairs from the published 'all' file
    ja = key_all.drop_duplicates(["gene_id", "snp_id"]).merge(mp, on=["gene_id", "snp_id"], how="inner",
                                                              suffixes=("_key", "_mine"))
    if len(ja) > 2:
        rho_s = float(stats.spearmanr(-np.log10(ja["pvalue"].clip(lower=1e-300)),
                                      -np.log10(ja["p"].clip(lower=1e-300)))[0])
    else:
        rho_s = float("nan")
    # eGene overlap (descriptive)
    pub_genes = set(key_best["gene_id"])
    tested = pub_genes & set(mine_best["gene_id"])
    mine_sig = set(mine_best.loc[mine_best.get("q_bh", pd.Series(dtype=float)) < fdr, "gene_id"]) if "q_bh" in mine_best else set()
    ok = n_dir >= min_pairs and np.isfinite(same) and same >= min_same_direction
    flags = []
    if n_dir and np.isfinite(same) and same <= 0.10:
        flags.append("same-direction rate near 0: possible global allele-orientation convention difference between "
                     "published rho and counted allele; investigate, do not auto-flip")
    return {
        "test": "V4", "population": population,
        "n_published_best_pairs_compared": n_dir, "same_direction_rate": same,
        "n_shared_pairs_for_rank_correlation": int(len(ja)), "spearman_neglog10p": rho_s,
        "n_published_egenes": len(pub_genes), "n_published_egenes_tested_here": len(tested),
        "n_egenes_called_here_fdr": len(mine_sig), "n_overlap_egenes": len(tested & mine_sig),
        "published_egene_recall_among_tested": (len(tested & mine_sig) / len(tested)) if tested else float("nan"),
        "fdr_threshold": fdr, "flags": flags,
        "criterion": f"same-direction rate >= {min_same_direction} over >= {min_pairs} compared published best pairs",
        "status": PASSED if ok else FAILED, "executed_at": now(),
    }


# ------------------------------------------------------------------------------------ V5
def v5_compare_to_plink2(mine: pd.DataFrame, plink: pd.DataFrame, beta_se_rtol: float = V5_BETA_SE_RTOL,
                         neglog10p_atol: float = V5_NEGLOG10P_ATOL, model: str = "linear") -> dict:
    """Compare VariantBridge results with ``plink2 --glm`` on an identical subset/model.

    ``mine``: variant, effect, se, p.   ``plink``: variant, effect, se, p (effect = BETA or
    LOG(OR)). Both are matched on variant ID; any unmatched variant is reported and fails the test.
    """
    j = mine.merge(plink, on="variant", suffixes=("_vb", "_pl"), how="outer", indicator=True)
    unmatched = int((j["_merge"] != "both").sum())
    b = j[j["_merge"] == "both"]
    def rel(a, c):
        return np.abs(a - c) / np.maximum(np.abs(c), 1e-12)
    d_beta = rel(b["effect_vb"], b["effect_pl"])
    d_se = rel(b["se_vb"], b["se_pl"])
    lp_vb = -np.log10(b["p_vb"].clip(lower=1e-300))
    lp_pl = -np.log10(b["p_pl"].clip(lower=1e-300))
    mask = lp_pl < 30  # PLINK2 prints limited precision; very small p are compared only qualitatively
    d_p = np.abs(lp_vb - lp_pl)[mask]
    res = {
        "test": "V5", "model": model, "n_variants_matched": int(len(b)), "n_unmatched": unmatched,
        "max_rel_diff_effect": float(d_beta.max()) if len(b) else float("nan"),
        "max_rel_diff_se": float(d_se.max()) if len(b) else float("nan"),
        "max_abs_diff_neglog10p": float(d_p.max()) if len(d_p) else float("nan"),
        "n_compared_p": int(mask.sum()),
        "tolerances": {"beta_se_rtol": beta_se_rtol, "neglog10p_atol": neglog10p_atol},
        "criterion": "all matched variants within tolerances; zero unmatched variants",
        "executed_at": now(),
    }
    ok = (len(b) > 0 and unmatched == 0 and res["max_rel_diff_effect"] <= beta_se_rtol
          and res["max_rel_diff_se"] <= beta_se_rtol and (len(d_p) == 0 or res["max_abs_diff_neglog10p"] <= neglog10p_atol))
    res["status"] = PASSED if ok else FAILED
    return res


# ------------------------------------------------------------------------------------ V6
def v6_expression_score(eval_result: dict, alignment_report: dict, implementation_checks: dict,
                        alpha: float = 0.05) -> dict:
    """Report scoring-pipeline implementation and biological predictive evidence SEPARATELY.

    ``status`` reflects ONLY the implementation checks. A null biological result does not fail
    the build (amendment 3); it is reported in ``biological_evidence``.
    """
    impl_ok = bool(implementation_checks) and all(bool(v) for v in implementation_checks.values())
    p = eval_result["empirical_p_one_sided"]
    lo = eval_result["bootstrap_ci95_mean_r"][0]
    evidence = bool(p < alpha and lo > 0)
    return {
        "test": "V6",
        "implementation": {"checks": {k: bool(v) for k, v in implementation_checks.items()},
                           "alignment_report": alignment_report,
                           "status": PASSED if impl_ok else FAILED},
        "biological_evidence": {
            "observed_mean_r": eval_result["observed_mean_r"],
            "bootstrap_ci95_mean_r": eval_result["bootstrap_ci95_mean_r"],
            "observed_mean_r2": eval_result["observed_mean_r2"],
            "n_genes_scored": eval_result["n_genes_scored"], "n_samples": eval_result["n_samples"],
            "null_ci95": eval_result["null_ci95"], "n_permutations": eval_result["n_permutations"],
            "empirical_p_one_sided": p, "alpha": alpha,
            "evidence_above_permutation_null": evidence,
            "interpretation": ("Observed predictive performance exceeds the permutation null (empirical p < alpha) "
                               "with a bootstrap CI above zero." if evidence else
                               "No evidence of predictive performance above the permutation null at this alpha. "
                               "This is a valid null biological result and does not indicate an implementation failure."),
            "caveat": "Discovery tissue is blood/PBMC; target tissue is LCL: validates scoring machinery, not cross-tissue biology.",
        },
        "status": PASSED if impl_ok else FAILED,
        "executed_at": now(),
    }
