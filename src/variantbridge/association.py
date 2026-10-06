"""Association models with single-model statistical consistency.

RELEASE-BLOCKING INVARIANT (handoff section 4, audit section C): for every reported
variant, the effect (beta or log-odds), standard error, confidence interval and p-value
must come from the SAME fitted model, including covariate adjustment.

Two engines are provided:

* linear (OLS, t-test, df = n - k - 2): exactly the model fitted by ``plink2 --glm`` for
  quantitative traits with covariates, and the linear model used by Matrix eQTL.
* logistic (IRLS maximum likelihood, Wald z-test): the model fitted by ``plink2 --glm``
  for binary traits when Firth fallback is disabled (``no-firth``). No Firth correction is
  implemented here; separation / non-convergence is reported as status != "ok" instead of
  being silently patched.

``scan_linear`` is a vectorised implementation (Frisch-Waugh-Lovell residualisation) that is
algebraically identical to fitting one OLS model per variant. ``check_statistical_consistency``
re-derives CI and p from (effect, SE, df) and fails loudly on any mismatch.

Capability tier: STATISTICALLY DEFENSIBLE engine; end-to-end benchmark against PLINK2 is
validation test V5 (see validation.py) and is NOT RUN until executed.
"""
from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .stats_utils import benjamini_hochberg, bonferroni

RESULT_COLUMNS = [
    "variant", "n", "effect", "se", "ci_low", "ci_high", "stat", "df", "p",
    "model", "status", "odds_ratio", "or_ci_low", "or_ci_high",
]


def _as_matrix(covariates, n: int) -> np.ndarray:
    if covariates is None:
        return np.empty((n, 0))
    c = covariates.to_numpy(dtype=float) if hasattr(covariates, "to_numpy") else np.asarray(covariates, dtype=float)
    if c.ndim == 1:
        c = c[:, None]
    if c.shape[0] != n:
        raise ValueError(f"covariates have {c.shape[0]} rows, expected {n}")
    if not np.isfinite(c).all():
        raise ValueError("covariates must not contain NaN/inf; drop or impute samples first")
    return c


def _row(variant, n, effect, se, df, model, status, alpha, logistic=False) -> dict:
    """Build one result row. Every statistic is derived here from (effect, se, df)."""
    if not (np.isfinite(effect) and np.isfinite(se) and se > 0):
        return dict(variant=variant, n=n, effect=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan,
                    stat=np.nan, df=df, p=np.nan, model=model,
                    status=status if status != "ok" else "degenerate",
                    odds_ratio=np.nan, or_ci_low=np.nan, or_ci_high=np.nan)
    stat = effect / se
    if logistic:
        p = float(2.0 * stats.norm.sf(abs(stat)))
        crit = float(stats.norm.isf(alpha / 2.0))
    else:
        p = float(2.0 * stats.t.sf(abs(stat), df))
        crit = float(stats.t.isf(alpha / 2.0, df))
    lo, hi = effect - crit * se, effect + crit * se
    return dict(variant=variant, n=n, effect=effect, se=se, ci_low=lo, ci_high=hi, stat=stat, df=df,
                p=p, model=model, status=status,
                odds_ratio=float(np.exp(effect)) if logistic else np.nan,
                or_ci_low=float(np.exp(lo)) if logistic else np.nan,
                or_ci_high=float(np.exp(hi)) if logistic else np.nan)


# ----------------------------------------------------------------------------- linear
def fit_ols(y, x, covariates=None, variant: str = "variant", alpha: float = 0.05) -> dict:
    """Fit y ~ 1 + covariates + x by OLS on complete cases; return one result row."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    cov = _as_matrix(covariates, y.shape[0])
    ok = np.isfinite(y) & np.isfinite(x)
    y, x, cov = y[ok], x[ok], cov[ok]
    n = y.shape[0]
    X = np.column_stack([np.ones(n), cov, x])
    k = X.shape[1]
    df = n - k
    if df < 1 or np.ptp(x) == 0:
        return _row(variant, n, np.nan, np.nan, df, "OLS", "degenerate", alpha)
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    sigma2 = float(resid @ resid) / df
    se = float(np.sqrt(sigma2 * XtX_inv[-1, -1]))
    return _row(variant, n, float(beta[-1]), se, df, "OLS", "ok", alpha)


def scan_linear(G, y, covariates=None, variant_ids: Optional[Sequence[str]] = None,
                alpha: float = 0.05) -> pd.DataFrame:
    """OLS association of ``y`` on every column of genotype matrix ``G`` (n x m).

    Variants with no missing genotypes are handled in one vectorised pass; variants with
    missing genotypes are refit on complete cases (exactly what plink2 --glm does), so every
    row equals ``fit_ols`` on that variant.
    """
    G = np.asarray(G, dtype=float)
    y = np.asarray(y, dtype=float)
    n, m = G.shape
    if y.shape[0] != n:
        raise ValueError("y and G have different numbers of samples")
    if not np.isfinite(y).all():
        raise ValueError("y contains NaN; drop samples with missing phenotype before scanning")
    ids = list(variant_ids) if variant_ids is not None else [f"v{j}" for j in range(m)]
    cov = _as_matrix(covariates, n)
    A = np.column_stack([np.ones(n), cov])
    # residualise y once
    Q, _ = np.linalg.qr(A)
    yr = y - Q @ (Q.T @ y)
    rows: list[dict] = [None] * m  # type: ignore[list-item]
    complete = np.isfinite(G).all(axis=0)
    idx = np.flatnonzero(complete)
    if idx.size:
        Gc = G[:, idx]
        Gr = Gc - Q @ (Q.T @ Gc)
        sxx = np.einsum("ij,ij->j", Gr, Gr)
        sxy = Gr.T @ yr
        syy = float(yr @ yr)
        df = n - (A.shape[1] + 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = sxy / sxx
            rss = syy - beta * sxy
            se = np.sqrt((rss / df) / sxx)
        degenerate = (sxx <= 1e-10 * max(n, 1)) | (df < 1)
        for pos, j in enumerate(idx):
            if degenerate[pos]:
                rows[j] = _row(ids[j], n, np.nan, np.nan, df, "OLS", "degenerate", alpha)
            else:
                rows[j] = _row(ids[j], n, float(beta[pos]), float(se[pos]), df, "OLS", "ok", alpha)
    for j in np.flatnonzero(~complete):
        rows[j] = fit_ols(y, G[:, j], covariates=cov, variant=ids[j], alpha=alpha)
    return _finalise(pd.DataFrame(rows))


# --------------------------------------------------------------------------- logistic
def fit_logistic(y, x, covariates=None, variant: str = "variant", alpha: float = 0.05,
                 max_iter: int = 50, tol: float = 1e-10) -> dict:
    """Logistic regression y ~ 1 + covariates + x by IRLS; Wald z-test on the x coefficient.

    Reports status "ok", "no_convergence" or "separation" (never silently corrected).
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    cov = _as_matrix(covariates, y.shape[0])
    ok = np.isfinite(y) & np.isfinite(x)
    y, x, cov = y[ok], x[ok], cov[ok]
    n = y.shape[0]
    if not set(np.unique(y)).issubset({0.0, 1.0}):
        raise ValueError("logistic phenotype must be coded 0/1")
    if n < 5 or np.ptp(x) == 0 or y.min() == y.max():
        return _row(variant, n, np.nan, np.nan, np.nan, "logistic-Wald", "degenerate", alpha, True)
    X = np.column_stack([np.ones(n), cov, x])
    k = X.shape[1]
    b = np.zeros(k)
    b[0] = np.log(y.mean() / (1 - y.mean()))
    prev_dev = np.inf
    status = "no_convergence"
    cov_mat = None
    for _ in range(max_iter):
        eta = np.clip(X @ b, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = mu * (1 - mu)
        if np.any(w < 1e-10) and np.max(np.abs(b)) > 15:
            status = "separation"
            break
        H = X.T @ (X * w[:, None])
        try:
            step = np.linalg.solve(H, X.T @ (y - mu))
        except np.linalg.LinAlgError:
            status = "separation"
            break
        b = b + step
        eta = np.clip(X @ b, -30, 30)
        mu = np.clip(1.0 / (1.0 + np.exp(-eta)), 1e-15, 1 - 1e-15)
        dev = -2.0 * float(np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu)))
        if abs(prev_dev - dev) < tol * (abs(dev) + 0.1):
            status = "ok"
            w = mu * (1 - mu)
            cov_mat = np.linalg.pinv(X.T @ (X * w[:, None]))
            break
        prev_dev = dev
    if status != "ok" or cov_mat is None:
        return _row(variant, n, np.nan, np.nan, np.nan, "logistic-Wald", status, alpha, True)
    effect = float(b[-1])
    se = float(np.sqrt(cov_mat[-1, -1]))
    if abs(effect) > 12 or se > 1e3:
        return _row(variant, n, np.nan, np.nan, np.nan, "logistic-Wald", "separation", alpha, True)
    return _row(variant, n, effect, se, np.nan, "logistic-Wald", "ok", alpha, True)


def scan_logistic(G, y, covariates=None, variant_ids: Optional[Sequence[str]] = None,
                  alpha: float = 0.05) -> pd.DataFrame:
    G = np.asarray(G, dtype=float)
    n, m = G.shape
    ids = list(variant_ids) if variant_ids is not None else [f"v{j}" for j in range(m)]
    cov = _as_matrix(covariates, n)
    rows = [fit_logistic(y, G[:, j], covariates=cov, variant=ids[j], alpha=alpha) for j in range(m)]
    return _finalise(pd.DataFrame(rows))


# ------------------------------------------------------------------------- dispatcher
def is_binary(y) -> bool:
    vals = set(pd.Series(np.asarray(y)).dropna().unique())
    return len(vals) > 0 and vals.issubset({0, 1, 0.0, 1.0})


def scan(G, y, covariates=None, variant_ids=None, alpha: float = 0.05, model: str = "auto") -> pd.DataFrame:
    """Dispatch to the linear or logistic engine. ``model`` in {"auto","linear","logistic"}."""
    if model == "auto":
        model = "logistic" if is_binary(y) else "linear"
    if model == "logistic":
        return scan_logistic(G, y, covariates, variant_ids, alpha)
    if model == "linear":
        return scan_linear(G, y, covariates, variant_ids, alpha)
    raise ValueError(f"unknown model {model!r}")


def _finalise(df: pd.DataFrame) -> pd.DataFrame:
    df = df[RESULT_COLUMNS].copy()
    df["neglog10p"] = -np.log10(df["p"].clip(lower=1e-300))
    df["bonferroni"] = bonferroni(df["p"].to_numpy())
    df["fdr_bh"] = benjamini_hochberg(df["p"].to_numpy())
    return df


def annotate_with_metadata(results: pd.DataFrame, meta: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Left-join variant metadata (chrom, pos, annotations); fill chrom/pos if absent."""
    d = results.copy()
    if meta is not None and "variant" in meta:
        d = d.merge(meta, on="variant", how="left")
    if "chrom" not in d:
        d["chrom"] = "NA"
    if "pos" not in d:
        d["pos"] = np.arange(len(d))
    return d


# ---------------------------------------------------------------- consistency checking
def check_statistical_consistency(results: pd.DataFrame, alpha: float = 0.05, rtol: float = 1e-8) -> list[str]:
    """Re-derive stat, CI and p from (effect, se, df) and report every violation.

    Returns a list of human-readable violations (empty list == consistent). This is the
    internal half of V5; the external half (PLINK2 benchmark) is a separate, NOT RUN-until-
    executed test.
    """
    problems: list[str] = []
    d = results[results["status"] == "ok"]
    for r in d.itertuples(index=False):
        logistic = str(r.model).startswith("logistic")
        stat = r.effect / r.se
        if logistic:
            p = 2 * stats.norm.sf(abs(stat))
            crit = stats.norm.isf(alpha / 2)
        else:
            p = 2 * stats.t.sf(abs(stat), r.df)
            crit = stats.t.isf(alpha / 2, r.df)
        if not np.isclose(r.stat, stat, rtol=rtol, atol=1e-12):
            problems.append(f"{r.variant}: stat != effect/se")
        if not np.isclose(r.p, p, rtol=rtol, atol=1e-300):
            problems.append(f"{r.variant}: p-value not derived from effect/se/df")
        if not (np.isclose(r.ci_low, r.effect - crit * r.se, rtol=rtol, atol=1e-12)
                and np.isclose(r.ci_high, r.effect + crit * r.se, rtol=rtol, atol=1e-12)):
            problems.append(f"{r.variant}: CI not derived from effect/se/df")
        if logistic and not np.isclose(r.odds_ratio, np.exp(r.effect), rtol=rtol):
            problems.append(f"{r.variant}: OR != exp(effect)")
    return problems


def assert_statistical_consistency(results: pd.DataFrame, **kw) -> None:
    problems = check_statistical_consistency(results, **kw)
    if problems:
        raise AssertionError("statistical-consistency violations: " + "; ".join(problems[:5]))
