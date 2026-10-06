"""Small, dependency-light statistical helpers shared across modules."""
from __future__ import annotations

import numpy as np
from scipy import stats

CHI2_1DF_MEDIAN = float(stats.chi2.ppf(0.5, 1))  # ~0.4549364


def p_to_chisq(p) -> np.ndarray:
    """Convert two-sided p-values to 1-df chi-square statistics."""
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1.0)
    return stats.chi2.isf(p, 1)


def genomic_inflation_lambda(p) -> float:
    """Genomic inflation factor (Devlin & Roeder 1999): median(chi2) / expected median.

    A diagnostic, not a correction. Returns NaN when no finite p-values are given.
    """
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return float("nan")
    return float(np.median(p_to_chisq(p)) / CHI2_1DF_MEDIAN)


def benjamini_hochberg(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values). NaNs are preserved."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * m / (np.arange(1, m + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def bonferroni(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = int(np.isfinite(p).sum())
    return np.minimum(p * m, 1.0)


def empirical_pvalue(observed: float, null, tail: str = "greater") -> float:
    """Permutation p-value with the +1 correction: (1 + #{null >= obs}) / (1 + B)."""
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or null.size == 0:
        return float("nan")
    if tail == "greater":
        k = np.sum(null >= observed)
    elif tail == "less":
        k = np.sum(null <= observed)
    else:
        raise ValueError("tail must be 'greater' or 'less'")
    return float((1 + k) / (1 + null.size))


def qq_points(p) -> tuple[np.ndarray, np.ndarray]:
    """Expected and observed -log10(p) for a QQ plot."""
    p = np.sort(np.asarray(p, dtype=float)[np.isfinite(p)])
    p = np.clip(p, 1e-300, 1.0)
    n = p.size
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    return expected, -np.log10(p)


def two_sided_t_p(t, df) -> np.ndarray:
    return 2.0 * stats.t.sf(np.abs(t), df)
