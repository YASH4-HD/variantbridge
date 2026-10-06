"""Sample/variant QC, exact HWE test, sex check, KING-robust kinship.

Genotype matrices are n_samples x n_variants arrays of additive dosage 0/1/2 with NaN for
missing. Thresholds are never hard-coded inside functions: they are explicit arguments, and
``QCThresholds`` records the project's chosen values together with their provenance label
([source ...] or [project design choice]) per audit section 5.

Capability tier: STATISTICALLY DEFENSIBLE (standard QC; Anderson 2010, Marees 2018).
Scalable production QC is done offline with PLINK2 (scripts/02_prepare_1000g_pca.py); these
functions are the reference implementations used for the app, small data, and unit tests.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class QCThresholds:
    sample_call_rate_min: float = 0.98      # [project design choice within 95-98% convention: Anderson 2010; Marees 2018]
    heterozygosity_sd: float = 3.0          # [project design choice: +/-3 SD convention, Anderson 2010; Marees 2018]
    variant_missing_max: float = 0.02       # [project design choice within 2-5% convention]
    maf_min: float = 0.01                   # [project design choice within 1-5% convention]
    hwe_p_min: float = 1e-6                 # [GWAS-only convention, controls only; NEVER pooled multi-ancestry, NEVER eQTL path]

    def provenance(self) -> dict:
        return asdict(self)


def sample_metrics(G) -> pd.DataFrame:
    G = np.asarray(G, dtype=float)
    called = np.isfinite(G)
    call_rate = called.mean(axis=1)
    het = np.where(called.sum(axis=1) > 0, (G == 1).sum(axis=1) / np.maximum(called.sum(axis=1), 1), np.nan)
    return pd.DataFrame({"call_rate": call_rate, "het_rate": het})


def variant_metrics(G) -> pd.DataFrame:
    G = np.asarray(G, dtype=float)
    called = np.isfinite(G)
    n_called = called.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        af = np.nansum(G, axis=0) / (2.0 * n_called)
    return pd.DataFrame({
        "missing_rate": 1.0 - n_called / G.shape[0],
        "alt_af": af,
        "maf": np.minimum(af, 1.0 - af),
        "n_called": n_called,
    })


def sample_qc_mask(metrics: pd.DataFrame, thr: QCThresholds) -> pd.Series:
    mu, sd = metrics["het_rate"].mean(), metrics["het_rate"].std(ddof=1)
    het_ok = (metrics["het_rate"] - mu).abs() <= thr.heterozygosity_sd * sd if sd > 0 else pd.Series(True, index=metrics.index)
    return (metrics["call_rate"] >= thr.sample_call_rate_min) & het_ok


def variant_qc_mask(metrics: pd.DataFrame, thr: QCThresholds, hwe_p: np.ndarray | None = None) -> np.ndarray:
    m = (metrics["missing_rate"] <= thr.variant_missing_max) & (metrics["maf"] >= thr.maf_min)
    if hwe_p is not None:
        m = m & (np.asarray(hwe_p) >= thr.hwe_p_min)
    return m.to_numpy()


# ------------------------------------------------------------------------------ HWE
def hwe_exact_pvalue(n_hom1: int, n_het: int, n_hom2: int) -> float:
    """Exact HWE test (Wigginton, Cutler & Abecasis 2005), two-sided."""
    n_hom1, n_het, n_hom2 = int(n_hom1), int(n_het), int(n_hom2)
    n = n_hom1 + n_het + n_hom2
    if n == 0:
        return float("nan")
    rare = 2 * min(n_hom1, n_hom2) + n_het
    n_rare_hom = min(n_hom1, n_hom2)
    n_common_hom = max(n_hom1, n_hom2)
    probs = np.zeros(rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if (mid % 2) != (rare % 2):
        mid += 1
    mid = min(mid, rare)
    probs[mid] = 1.0
    s = 1.0
    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets > 1:
        probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        s += probs[curr_hets - 2]
        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1
    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets <= rare - 2:
        probs[curr_hets + 2] = probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        s += probs[curr_hets + 2]
        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1
    probs /= s
    obs = n_het
    p = probs[probs <= probs[obs] * (1 + 1e-9)].sum()   # relative tolerance; extreme p may underflow to 0
    return float(min(1.0, p))


def hwe_pvalues(G) -> np.ndarray:
    G = np.asarray(G, dtype=float)
    out = np.empty(G.shape[1])
    for j in range(G.shape[1]):
        g = G[:, j]
        g = g[np.isfinite(g)]
        out[j] = hwe_exact_pvalue((g == 0).sum(), (g == 1).sum(), (g == 2).sum())
    return out


# ------------------------------------------------------------------------- sex check
def infer_sex_from_x_het(x_het_rate, threshold: float) -> np.ndarray:
    """Infer sex from X-chromosome heterozygosity: 2 = female (high het), 1 = male (low het).

    ``threshold`` is explicit by design: it is a [project design choice] that must be set and
    recorded by the caller from the observed bimodal distribution, not guessed here.
    """
    x = np.asarray(x_het_rate, dtype=float)
    return np.where(x >= threshold, 2, 1)


def sex_concordance(inferred, labelled) -> dict:
    inferred = np.asarray(inferred)
    labelled = np.asarray(labelled)
    ok = np.isin(labelled, [1, 2])
    n = int(ok.sum())
    k = int((inferred[ok] == labelled[ok]).sum())
    return {"n_compared": n, "n_concordant": k, "concordance": (k / n) if n else float("nan")}


# --------------------------------------------------------------------------- kinship
KING_BINS = [  # lower bounds of KING kinship coefficient (Manichaikul et al. 2010)
    (0.354, "duplicate/MZ twin"),
    (0.177, "1st-degree"),
    (0.0884, "2nd-degree"),
    (0.0442, "3rd-degree"),
]


def king_robust(G) -> np.ndarray:
    """KING-robust pairwise kinship (Manichaikul 2010): (N_AaAa - 2 N_AAaa) / (N_Aa(i) + N_Aa(j)).

    Het counts in the denominator are taken over SNPs genotyped in both samples. Diagonal is
    set to 0.5 by convention. Dense n x n; fine for n up to a few thousand.
    """
    G = np.asarray(G, dtype=float)
    called = np.isfinite(G).astype(np.float64)
    het = (G == 1).astype(np.float64)
    h0 = (G == 0).astype(np.float64)
    h2 = (G == 2).astype(np.float64)
    n_hh = het @ het.T
    n_opp = h0 @ h2.T
    n_opp = n_opp + n_opp.T
    n_het_i = het @ called.T          # hets of i over SNPs where j is also called (i,j)
    denom = n_het_i + n_het_i.T
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = (n_hh - 2.0 * n_opp) / denom
    np.fill_diagonal(phi, 0.5)
    return phi


def classify_kinship(phi: float) -> str:
    for lo, label in KING_BINS:
        if phi > lo:
            return label
    return "unrelated"


def related_pairs(phi: np.ndarray, sample_ids, min_kinship: float = 0.0442) -> pd.DataFrame:
    iu = np.triu_indices_from(phi, k=1)
    mask = phi[iu] > min_kinship
    return pd.DataFrame({
        "id1": np.asarray(sample_ids)[iu[0][mask]],
        "id2": np.asarray(sample_ids)[iu[1][mask]],
        "kinship": phi[iu][mask],
        "degree": [classify_kinship(v) for v in phi[iu][mask]],
    }).sort_values("kinship", ascending=False).reset_index(drop=True)
