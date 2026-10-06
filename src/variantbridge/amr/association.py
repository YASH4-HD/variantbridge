"""Unbiased, genome-wide association with EXPLICIT MIC censoring (addendum Part 2, items 2-3).

Primary model (user-directed, replaces a censoring-naive linear model): interval-censored Gaussian
('Tobit-type') regression of log2 MIC,

    log2 MIC*_i = x_i' beta + eps_i,  eps ~ N(0, sigma^2),   observed in [lo_i, hi_i]
    exact: lo = hi;  left-censored (<= b): (-inf, b];  right-censored (>= b): [b, inf)

fitted per feature by Newton's method with analytic derivatives. Effect (beta_g, log2-MIC change per
allele), SE, CI, z and Wald p all come from the SAME fitted model (fixes the v0.1 mismatch).
A likelihood-ratio p is reported as a labelled DIAGNOSTIC column and never replaces the Wald p.

Rules enforced here:
  * No biological pre-filtering: the feature space is whatever is passed; this module takes no
    annotation / answer-key / pathway input and its signature is tested for that.
  * Only the discovery tier may be analysed (assert_tier); replication uses `replicate_features`.
  * Lineage control = ST fixed effects (structure.lineage_design); an unadjusted fit is also
    reported (`unadjusted_p`) for transparency.
  * Multiple testing: BH-FDR (primary) + Bonferroni over the features actually tested.
  * Breakpoint-based binary models run only via assign_binary (VERIFIED breakpoints), else
    BreakpointUnverified (NOT RUN).
"""
from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import log_ndtr

from ..association import scan as core_scan
from ..stats_utils import benjamini_hochberg, bonferroni
from .cohorts import assert_tier
from .config import threshold
from .phenotypes import assign_binary, censoring_summary

MODEL_NAME = "interval_censored_gaussian_Wald"
RESULT_COLUMNS = ["feature_id", "feature_type", "n", "n_present", "effect", "se", "ci_low", "ci_high", "z", "p_value",
                  "p_lrt_diagnostic", "unadjusted_effect", "unadjusted_p", "sigma", "model", "status",
                  "q_value", "bonferroni_p", "lineage_adjusted", "structure_control"]


# ----------------------------------------------------------------------------------------- likelihood
def _terms(mu, s, lo, hi, kind):
    with np.errstate(all="ignore"):       # extreme trial points in line searches may overflow; callers reject non-finite values
        return _terms_raw(mu, s, lo, hi, kind)


def _terms_raw(mu, s, lo, hi, kind):
    """Per-observation loglik, d/dmu, d/ds, d2/dmu2, d2/dmu ds, d2/ds2 (s = log sigma)."""
    sig = np.exp(s)
    ll = np.empty_like(mu); gm = np.empty_like(mu); gs = np.empty_like(mu)
    hmm = np.empty_like(mu); hms = np.empty_like(mu); hss = np.empty_like(mu)
    e = kind == 0
    if e.any():
        z = (lo[e] - mu[e]) / sig
        ll[e] = -s - 0.5 * z * z - 0.5 * math.log(2 * math.pi)
        gm[e] = z / sig; gs[e] = z * z - 1.0
        hmm[e] = -1.0 / sig ** 2; hms[e] = -2.0 * z / sig; hss[e] = -2.0 * z * z
    for code, sign in ((-1, -1.0), (1, 1.0)):       # left: v=(hi-mu)/sig ; right: v=(mu-lo)/sig
        m = kind == code
        if not m.any():
            continue
        b = hi[m] if code == -1 else lo[m]
        v = sign * (mu[m] - b) / sig if code == 1 else (b - mu[m]) / sig
        lg = log_ndtr(v)
        vv = np.clip(v, -37.0, 37.0)            # |v|>37: Mills ratio asymptote; avoids inf-inf
        lam = np.exp(-0.5 * vv * vv - 0.5 * math.log(2 * math.pi) - log_ndtr(vv))
        c = v + lam
        ll[m] = lg
        dv_dmu = (1.0 / sig) if code == 1 else (-1.0 / sig)
        gm[m] = lam * dv_dmu
        gs[m] = -lam * v
        hmm[m] = -lam * c / sig ** 2
        hms[m] = (-(lam / sig) * (1.0 - v * c)) if code == 1 else ((lam / sig) * (1.0 - v * c))
        hss[m] = lam * v * (1.0 - v * c)
    return ll, gm, gs, hmm, hms, hss


def _kinds(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    k = np.zeros(len(lo), dtype=int)
    k[np.isneginf(lo)] = -1
    k[np.isposinf(hi)] = 1
    if np.any(np.isneginf(lo) & np.isposinf(hi)):
        raise ValueError("observation censored on both sides (uninformative)")
    return k


def fit_censored(lo, hi, X, beta0=None, s0=None, max_iter: int = 100, tol: float = 1e-8) -> dict:
    """Maximum-likelihood fit. Returns beta, log_sigma, loglik, cov (for [beta, log_sigma]), status."""
    lo = np.asarray(lo, float); hi = np.asarray(hi, float); X = np.asarray(X, float)
    n, k = X.shape
    kind = _kinds(lo, hi)
    if beta0 is None:
        # start: OLS on exact rows, else on boundary midpoints of finite ends (START ONLY; never the reported estimate)
        y0 = np.where(kind == 0, lo, np.where(kind == -1, hi - 1.0, lo + 1.0))
        beta0 = np.linalg.lstsq(X, y0, rcond=None)[0]
        resid = y0 - X @ beta0
        s0 = math.log(max(float(np.std(resid)), 1e-3))
    theta = np.concatenate([beta0, [s0]]).astype(float)

    def evaluate(th):
        b, s = th[:k], th[k]
        mu = X @ b
        with np.errstate(all="ignore"):
            ll, gm, gs, hmm, hms, hss = _terms(mu, s, lo, hi, kind)
        L = ll.sum()
        with np.errstate(all="ignore"):
            g = np.concatenate([X.T @ gm, [gs.sum()]])
            H = np.empty((k + 1, k + 1))
            H[:k, :k] = (X * hmm[:, None]).T @ X
            H[:k, k] = H[k, :k] = X.T @ hms
            H[k, k] = hss.sum()
        return L, g, H

    L, g, H = evaluate(theta)
    status = "ok"
    lam_damp = 0.0
    for it in range(max_iter):
        if np.linalg.norm(g, np.inf) < tol * max(1.0, n):
            break
        A = -H + lam_damp * np.eye(k + 1)
        try:
            step = np.linalg.solve(A + 1e-10 * np.eye(k + 1), g)
        except np.linalg.LinAlgError:
            status = "degenerate"; break
        if not np.all(np.isfinite(step)):
            status = "degenerate"; break
        t = 1.0
        while t > 1e-6:
            th2 = theta + t * step
            L2, g2, H2 = evaluate(th2)
            if np.isfinite(L2) and np.all(np.isfinite(g2)) and np.all(np.isfinite(H2)) and L2 >= L - 1e-12:
                break
            t *= 0.5
        else:
            lam_damp = max(1e-3, lam_damp * 10)
            if lam_damp > 1e6:
                status = "no_convergence"; break
            continue
        theta, L, g, H = th2, L2, g2, H2
        lam_damp = max(0.0, lam_damp / 10)
    else:
        status = "no_convergence"
    if status == "ok" and np.linalg.norm(g, np.inf) > 1e-3 * max(1.0, n):
        status = "no_convergence"
    cov = np.full((k + 1, k + 1), np.nan)
    try:
        cov = np.linalg.inv(-H)
        if np.any(np.diag(cov) <= 0) or not np.all(np.isfinite(cov)):
            status = status if status != "ok" else "degenerate"
    except np.linalg.LinAlgError:
        status = "degenerate"
    if status == "ok" and (abs(theta[:k]).max() > 50 or abs(theta[k]) > 10):
        status = "separation_or_extreme"
    return {"beta": theta[:k], "log_sigma": float(theta[k]), "loglik": float(L), "cov": cov, "status": status, "n": n}


# ----------------------------------------------------------------------------------------- scanning
def _feature_is_binary(g: np.ndarray) -> bool:
    u = np.unique(g[np.isfinite(g)])
    return len(u) > 0 and set(u.tolist()).issubset({0.0, 1.0})


def run_association(features: np.ndarray, feature_ids: Sequence[str], feature_types: Sequence[str] | str,
                    pheno: pd.DataFrame, lineage: Optional[np.ndarray], cfg: dict,
                    answer_key_sha256: str, structure_control: str = "lineage_fixed_effects") -> pd.DataFrame:
    """Discovery-tier genome-wide association.

    `features`: (n_samples x n_features), rows aligned to `pheno` rows. `pheno` needs censoring
    columns (log2_lower, log2_upper, censoring) and cohort_tier == discovery for every row.
    `answer_key_sha256` must be supplied (computed from the committed pre-registered file BEFORE this
    call); it is stored on the result so the ordering is auditable. It is not used otherwise."""
    if not answer_key_sha256 or len(answer_key_sha256) != 64:
        raise ValueError("answer_key_sha256 (64 hex chars) is required: the answer key must be registered before association")
    assert_tier(pheno, "discovery", "association.run_association")
    G = np.asarray(features, float)
    if G.shape != (len(pheno), len(feature_ids)):
        raise ValueError(f"features shape {G.shape} != ({len(pheno)}, {len(feature_ids)})")
    ftypes = [feature_types] * len(feature_ids) if isinstance(feature_types, str) else list(feature_types)
    lo = pheno["log2_lower"].to_numpy(float); hi = pheno["log2_upper"].to_numpy(float)
    n = len(pheno)
    min_iso = threshold(cfg, "features.min_isolates_with_feature")
    maf_min = threshold(cfg, "features.maf_min")
    alpha = threshold(cfg, "association.alpha_fdr")
    cov = np.zeros((n, 0)) if lineage is None else np.asarray(lineage, float)
    base = np.column_stack([np.ones(n), cov])
    null_fit = fit_censored(lo, hi, base)
    if null_fit["status"] != "ok":
        raise RuntimeError(f"null (covariate-only) censored model failed: {null_fit['status']}")
    ones = np.ones((n, 1))
    null_unadj = fit_censored(lo, hi, ones)
    rows, filtered = [], []
    for j, fid in enumerate(feature_ids):
        g = G[:, j]
        if not np.all(np.isfinite(g)):
            filtered.append((fid, "missing_values")); continue
        if np.ptp(g) == 0:
            filtered.append((fid, "constant")); continue
        if _feature_is_binary(g):
            npres = int(g.sum())
            if npres < min_iso or (n - npres) < min_iso or min(npres, n - npres) / n < maf_min:
                filtered.append((fid, "frequency_filter")); continue
        else:
            npres = int((g != 0).sum())
        Xj = np.column_stack([base, g])
        f = fit_censored(lo, hi, Xj, beta0=np.concatenate([null_fit["beta"], [0.0]]), s0=null_fit["log_sigma"])
        fu = fit_censored(lo, hi, np.column_stack([ones, g]), beta0=np.array([null_unadj["beta"][0], 0.0]), s0=null_unadj["log_sigma"])
        bi, se_i = f["beta"][-1], math.sqrt(f["cov"][-2, -2]) if f["cov"][-2, -2] > 0 else float("nan")
        z = bi / se_i if se_i and np.isfinite(se_i) else float("nan")
        p = 2 * stats.norm.sf(abs(z)) if np.isfinite(z) else float("nan")
        zc = stats.norm.ppf(1 - 0.05 / 2)
        lrt = max(0.0, 2 * (f["loglik"] - null_fit["loglik"]))
        seu = math.sqrt(fu["cov"][-2, -2]) if fu["cov"][-2, -2] > 0 else float("nan")
        zu = fu["beta"][-1] / seu if np.isfinite(seu) and seu > 0 else float("nan")
        rows.append(dict(feature_id=fid, feature_type=ftypes[j], n=n, n_present=npres, effect=bi, se=se_i,
                         ci_low=bi - zc * se_i, ci_high=bi + zc * se_i, z=z, p_value=p,
                         p_lrt_diagnostic=stats.chi2.sf(lrt, 1), unadjusted_effect=fu["beta"][-1],
                         unadjusted_p=2 * stats.norm.sf(abs(zu)) if np.isfinite(zu) else float("nan"),
                         sigma=math.exp(f["log_sigma"]), model=MODEL_NAME, status=f["status"]))
    res = pd.DataFrame(rows, columns=[c for c in RESULT_COLUMNS if c in rows[0]] if rows else RESULT_COLUMNS)
    if res.empty:
        res = pd.DataFrame(columns=RESULT_COLUMNS)
    ok = res["status"] == "ok" if len(res) else pd.Series([], dtype=bool)
    res["q_value"] = np.nan; res["bonferroni_p"] = np.nan
    if ok.any():
        res.loc[ok, "q_value"] = benjamini_hochberg(res.loc[ok, "p_value"].to_numpy())
        res.loc[ok, "bonferroni_p"] = bonferroni(res.loc[ok, "p_value"].to_numpy())
    res["lineage_adjusted"] = cov.shape[1] > 0
    res["structure_control"] = structure_control if cov.shape[1] > 0 else "none"
    res.attrs.update({"filtered": filtered, "answer_key_sha256": answer_key_sha256, "alpha_fdr": alpha,
                      "censoring": censoring_summary(pheno), "n_features_input": len(feature_ids),
                      "n_features_tested": int(len(res)), "n_features_filtered": len(filtered),
                      "null_loglik": null_fit["loglik"], "correction_method": "bh_fdr"})
    return res.sort_values("p_value").reset_index(drop=True)


def check_amr_consistency(res: pd.DataFrame, rtol: float = 1e-6) -> list[str]:
    """Release-gate style check: z, p, CI recomputed from (effect, SE) must match the stored values."""
    bad = []
    ok = res[res["status"] == "ok"]
    zc = stats.norm.ppf(1 - 0.05 / 2)
    for _, r in ok.iterrows():
        z = r["effect"] / r["se"]
        if not np.isclose(z, r["z"], rtol=rtol):
            bad.append(f"{r['feature_id']}: z inconsistent")
        if not np.isclose(2 * stats.norm.sf(abs(z)), r["p_value"], rtol=1e-4, atol=1e-300):
            bad.append(f"{r['feature_id']}: p inconsistent with effect/SE")
        if not (np.isclose(r["effect"] - zc * r["se"], r["ci_low"], rtol=rtol) and np.isclose(r["effect"] + zc * r["se"], r["ci_high"], rtol=rtol)):
            bad.append(f"{r['feature_id']}: CI inconsistent")
        if not (r["ci_low"] <= r["effect"] <= r["ci_high"]):
            bad.append(f"{r['feature_id']}: effect outside CI")
    return bad


def run_binary_verified_breakpoint(features, feature_ids, pheno, lineage, cfg, cohort="discovery") -> pd.DataFrame:
    """Breakpoint-dependent logistic association. NOT RUN (BreakpointUnverified) unless config says VERIFIED."""
    assert_tier(pheno, "discovery", "association.run_binary_verified_breakpoint")
    lab = assign_binary(pheno, cfg, cohort)
    keep = lab.isin(["susceptible", "resistant"]).to_numpy()
    y = (lab[keep] == "resistant").astype(float).to_numpy()
    cov = None if lineage is None else np.asarray(lineage, float)[keep]
    res = core_scan(np.asarray(features, float)[keep], y, cov, variant_ids=list(feature_ids), model="logistic")
    res["censoring_handling"] = "binary_verified_breakpoint"
    return res


# ----------------------------------------------------------------------------------------- replication
def replicate_features(discovery: pd.DataFrame, features: np.ndarray, feature_ids: Sequence[str],
                       pheno: pd.DataFrame, lineage: Optional[np.ndarray], cfg: dict, q_max: Optional[float] = None) -> pd.DataFrame:
    """Test PRE-SPECIFIED discovery hits (q < q_max) in the replication tier with the same model.

    Reports direction concordance and effect presence; a hit is `replicated` if the direction matches
    and the replication p < 0.05 / (number of discovery hits tested) [Bonferroni over hits].
    Features absent from the replication feature set are `not_testable` (never dropped silently)."""
    assert_tier(pheno, "replication", "association.replicate_features")
    q_max = threshold(cfg, "association.alpha_fdr") if q_max is None else q_max
    hits = discovery[(discovery["status"] == "ok") & (discovery["q_value"] < q_max)]
    idx = {f: i for i, f in enumerate(feature_ids)}
    n = len(pheno)
    lo = pheno["log2_lower"].to_numpy(float); hi = pheno["log2_upper"].to_numpy(float)
    cov = np.zeros((n, 0)) if lineage is None else np.asarray(lineage, float)
    base = np.column_stack([np.ones(n), cov])
    nf = fit_censored(lo, hi, base)
    out = []
    zc = stats.norm.ppf(0.975)
    for _, h in hits.iterrows():
        fid = h["feature_id"]
        if fid not in idx:
            out.append(dict(feature_id=fid, discovery_effect=h["effect"], discovery_q=h["q_value"], status="not_testable")); continue
        g = np.asarray(features, float)[:, idx[fid]]
        if np.ptp(g) == 0 or not np.all(np.isfinite(g)):
            out.append(dict(feature_id=fid, discovery_effect=h["effect"], discovery_q=h["q_value"], status="not_testable")); continue
        f = fit_censored(lo, hi, np.column_stack([base, g]), beta0=np.concatenate([nf["beta"], [0.0]]), s0=nf["log_sigma"])
        b = f["beta"][-1]; se = math.sqrt(f["cov"][-2, -2]) if f["cov"][-2, -2] > 0 else float("nan")
        z = b / se if np.isfinite(se) and se > 0 else float("nan")
        out.append(dict(feature_id=fid, discovery_effect=h["effect"], discovery_q=h["q_value"], effect=b, se=se,
                        ci_low=b - zc * se, ci_high=b + zc * se, z=z, p_value=2 * stats.norm.sf(abs(z)) if np.isfinite(z) else float("nan"),
                        same_direction=bool(np.sign(b) == np.sign(h["effect"])), n_present=int((g != 0).sum()), status=f["status"]))
    r = pd.DataFrame(out)
    if len(r):
        thr = 0.05 / max(1, len(hits))
        r["replicated"] = (r.get("status") == "ok") & r.get("same_direction", False).astype(bool) & (r.get("p_value", np.nan) < thr)
        r.attrs["bonferroni_threshold_over_hits"] = thr
    return r


# ----------------------------------------------------------------------------------------- pyseer wrapper
class PyseerUnavailable(RuntimeError):
    """pyseer is not installed/on PATH => the pyseer stage is NOT RUN (never silently replaced)."""


PYSEER_MODES = ("binary_verified_breakpoint", "continuous_censoring_naive_SENSITIVITY_ONLY")


def pyseer_available() -> bool:
    import shutil
    return shutil.which("pyseer") is not None


def write_pyseer_phenotype(pheno: pd.DataFrame, cfg: dict, mode: str, path, acknowledge_censoring_naive: bool = False,
                           cohort: str = "discovery") -> dict:
    """Write pyseer's phenotype file (samples<TAB>phenotype). Returns {'path', 'censoring_handling', 'n'}.

    binary_verified_breakpoint: needs VERIFIED breakpoints (BreakpointUnverified otherwise); indeterminate/intermediate
        isolates are dropped and COUNTED.
    continuous_censoring_naive_SENSITIVITY_ONLY: boundary values stand in for censored MICs. This is the naive analysis the
        project forbids as a primary result; it is only available with acknowledge_censoring_naive=True, is labelled in the
        output, and an Evidence Card may not use it as association evidence."""
    assert_tier(pheno, "discovery", "association.run_pyseer")
    if mode == "binary_verified_breakpoint":
        lab = assign_binary(pheno, cfg, cohort)
        keep = lab.isin(["susceptible", "resistant"])
        out = pd.DataFrame({"samples": pheno.loc[keep, "genome_id"], "phenotype": (lab[keep] == "resistant").astype(int)})
        handling = "binary_verified_breakpoint"
        dropped = int((~keep).sum())
    elif mode == "continuous_censoring_naive_SENSITIVITY_ONLY":
        if not acknowledge_censoring_naive:
            raise ValueError("continuous pyseer mode treats censored MIC boundaries as measurements; pass acknowledge_censoring_naive=True (sensitivity only)")
        y = np.where(np.isfinite(pheno["log2_mic_exact"]), pheno["log2_mic_exact"], pheno["log2_boundary"])
        out = pd.DataFrame({"samples": pheno["genome_id"], "phenotype": y})
        handling = "naive_boundary_substitution_SENSITIVITY_ONLY"
        dropped = 0
    else:
        raise ValueError(f"mode must be one of {PYSEER_MODES}")
    out.to_csv(path, sep="\t", index=False)
    return {"path": str(path), "censoring_handling": handling, "n": int(len(out)), "dropped": dropped}


def build_pyseer_command(feature_kind: str, features_path: str, pheno_path: str, similarity_path: str | None = None,
                         covariates_path: str | None = None, min_af: float = 0.01, cpu: int = 1, continuous: bool = False) -> list[str]:
    """pyseer argument list. Flag names follow the pyseer documentation as recalled and are UNVERIFIED in this build
    (pyseer could not be installed); the first real run must be checked against `pyseer --help`."""
    if feature_kind not in ("unitig", "gene_pa"):
        raise ValueError("feature_kind must be 'unitig' or 'gene_pa'")
    cmd = ["pyseer", "--phenotypes", pheno_path, "--kmers" if feature_kind == "unitig" else "--pres", features_path,
           "--min-af", str(min_af), "--max-af", str(1 - min_af), "--cpu", str(cpu)]
    if similarity_path:
        cmd += ["--lmm", "--similarity", similarity_path]
    if covariates_path:
        cmd += ["--covariates", covariates_path, "--use-covariates"]
    if continuous:
        cmd += ["--continuous"]
    return cmd


def parse_pyseer_output(df: pd.DataFrame, feature_type: str, n: int, censoring_handling: str, structure_control: str) -> pd.DataFrame:
    """pyseer TSV -> association_results schema. Expected columns (UNVERIFIED): variant, af, filter-pvalue, lrt-pvalue,
    beta, beta-std-err, notes. Effect/SE/CI/p are made consistent: p_value is the WALD p from beta/SE (the LRT p is kept
    as a diagnostic), so all reported statistics describe the same estimate. Rows pyseer flagged in `notes` are kept with
    status 'pyseer_flag:<notes>' and excluded from FDR."""
    need = {"variant", "beta", "beta-std-err", "lrt-pvalue"}
    if not need <= set(df.columns):
        raise KeyError(f"pyseer output lacks columns {sorted(need - set(df.columns))} (format UNVERIFIED)")
    zc = stats.norm.ppf(0.975)
    out = pd.DataFrame({"feature_id": df["variant"].astype(str), "feature_type": feature_type, "n": n,
                        "effect": df["beta"].astype(float), "se": df["beta-std-err"].astype(float)})
    out["ci_low"] = out["effect"] - zc * out["se"]; out["ci_high"] = out["effect"] + zc * out["se"]
    out["z"] = out["effect"] / out["se"]
    out["p_value"] = 2 * stats.norm.sf(out["z"].abs())
    out["p_lrt_diagnostic"] = df["lrt-pvalue"].astype(float)
    notes = df["notes"].fillna("") if "notes" in df else pd.Series("", index=df.index)
    out["status"] = np.where(notes.astype(str).str.len() > 0, "pyseer_flag:" + notes.astype(str), "ok")
    out["model"] = "pyseer_" + ("LMM" if structure_control == "lmm" else "fixed") + "_Wald"
    out["censoring_handling"] = censoring_handling
    out["structure_control"] = structure_control
    ok = out["status"] == "ok"
    out["q_value"] = np.nan
    out.loc[ok, "q_value"] = benjamini_hochberg(out.loc[ok, "p_value"].to_numpy())
    return out.sort_values("p_value").reset_index(drop=True)


def run_pyseer(feature_kind: str, features_path: str, pheno: pd.DataFrame, cfg: dict, workdir, mode: str,
               similarity_path: str | None = None, acknowledge_censoring_naive: bool = False, runner=None) -> pd.DataFrame:
    """Subprocess wrapper (addendum interface: accepts ONLY the discovery tier). `runner` is injectable for tests.
    Raises PyseerUnavailable when pyseer is absent (=> NOT RUN)."""
    import subprocess
    from pathlib import Path
    assert_tier(pheno, "discovery", "association.run_pyseer")
    if runner is None:
        if not pyseer_available():
            raise PyseerUnavailable("pyseer not found on PATH: pyseer stage NOT RUN")
        runner = subprocess.run
    wd = Path(workdir); wd.mkdir(parents=True, exist_ok=True)
    meta = write_pyseer_phenotype(pheno, cfg, mode, wd / "pyseer_phenotype.tsv", acknowledge_censoring_naive)
    cmd = build_pyseer_command(feature_kind, features_path, meta["path"], similarity_path, min_af=threshold(cfg, "features.maf_min"),
                               continuous=mode.startswith("continuous"))
    out = wd / "pyseer_output.tsv"
    with open(out, "w") as fh:
        r = runner(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if getattr(r, "returncode", 0) != 0:
        raise RuntimeError(f"pyseer failed: {getattr(r, 'stderr', '')[-500:]}")
    return parse_pyseer_output(pd.read_csv(out, sep="\t"), feature_kind, meta["n"], meta["censoring_handling"], "lmm" if similarity_path else "none")
