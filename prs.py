"""Polygenic score utilities for REAL complex traits only (leakage-controlled).

PRS terminology is reserved for complex-trait prediction. The eQTLGen -> GEUVADIS workflow is
a Genetic Expression Score (expression_score.py) and never uses this module's wording.

Audit Stage 7(d): a rigorous disease-trait PRS with fully open individual-level genotypes and
phenotypes is not currently possible; it needs a controlled-access target (e.g. UK Biobank).
No open-data PRS validation is run by VariantBridge; this module supports (i) the clearly
SIMULATED method demonstration in the app and (ii) future controlled-access use.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from .stats_utils import empirical_pvalue

DISEASE_PRS_LIMITATION = (
    "A rigorous disease-trait PRS requires independent discovery summary statistics and an "
    "independent, ancestry-matched target cohort with phenotypes. Fully open individual-level "
    "genotype+phenotype data of that kind do not currently exist, so VariantBridge makes no "
    "disease-PRS claim. Any score shown for simulated data is a method demonstration only."
)


def assert_disjoint_cohorts(discovery_ids, target_ids) -> None:
    """Fail loudly if any sample identifier occurs in both discovery and target sets."""
    overlap = set(map(str, discovery_ids)) & set(map(str, target_ids))
    if overlap:
        raise ValueError(f"discovery/target overlap ({len(overlap)} samples): leakage; refusing to score")


def polygenic_score(G, weights) -> np.ndarray:
    G = np.asarray(G, dtype=float)
    G = np.where(np.isfinite(G), G, np.nanmean(G, axis=0))
    return G @ np.asarray(weights, dtype=float)


def select_weights(discovery_results: pd.DataFrame, p_threshold: float, max_variants: int | None = None) -> pd.DataFrame:
    """Select discovery variants by a threshold FIXED A PRIORI (never tuned on the target)."""
    d = discovery_results[(discovery_results["status"] == "ok") & (discovery_results["p"] <= p_threshold)]
    d = d.sort_values("p")
    return d.head(max_variants) if max_variants else d


def evaluate_polygenic_score(score, y, covariates=None, n_permutations: int = 1000, n_bootstrap: int = 500,
                             seed: int = 0) -> dict:
    """Target-cohort performance with uncertainty and a permutation null.

    Quantitative y: incremental R^2 of the score over covariates. Binary y: AUC of the score.
    """
    rng = np.random.default_rng(seed)
    score, y = np.asarray(score, float), np.asarray(y, float)
    binary = set(np.unique(y)).issubset({0.0, 1.0})
    if covariates is None:
        cov = np.ones((len(y), 1))
    else:
        cov = np.column_stack([np.ones(len(y)), np.asarray(covariates, float)])

    def metric(s, yy, cc):
        if binary:
            return roc_auc_score(yy, s)
        r0 = yy - cc @ np.linalg.lstsq(cc, yy, rcond=None)[0]
        X = np.column_stack([cc, s])
        r1 = yy - X @ np.linalg.lstsq(X, yy, rcond=None)[0]
        return 1.0 - (r1 @ r1) / (r0 @ r0)

    obs = float(metric(score, y, cov))
    null = np.array([metric(score, y[rng.permutation(len(y))], cov) for _ in range(n_permutations)])
    boot = []
    for _ in range(n_bootstrap):
        i = rng.integers(0, len(y), len(y))
        if binary and len(np.unique(y[i])) < 2:
            continue
        boot.append(metric(score[i], y[i], cov[i]))
    return {"metric": "AUC" if binary else "incremental R2", "observed": obs,
            "bootstrap_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            "null_ci95": [float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))],
            "empirical_p_one_sided": empirical_pvalue(obs, null, "greater"),
            "n_permutations": n_permutations, "seed": seed}
