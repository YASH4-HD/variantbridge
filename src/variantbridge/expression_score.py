"""Genetic Expression Score / Predicted Expression (eQTLGen weights -> GEUVADIS target).

TERMINOLOGY (auditor amendment 1): this workflow predicts GENE EXPRESSION from genotype. It is
NOT a polygenic risk score and must never be labelled "PRS" or "risk". PRS terminology is
reserved for real complex traits (see prs.py).

Design (audit Stage 7 route a):
  * discovery weights: eQTLGen cis summary statistics (z-scores) for blood/PBMC;
  * target: GEUVADIS EUR genotypes + measured expression (sample-disjoint from eQTLGen);
  * per gene: select weights with p <= p_threshold, align alleles, drop palindromic SNPs,
    clump greedily by p with r^2 <= r2_clump measured in the target genotypes (genotype-only,
    no phenotype), score = sum(z_aligned * dosage);
  * evaluation (V6): per-gene Pearson r against measured expression, aggregated statistic,
    bootstrap CI over samples, and a PERMUTATION NULL (sample labels of measured expression
    permuted jointly across genes) with an empirical p-value.

LEAKAGE CONTROLS: p_threshold and r2_clump are fixed a priori (arguments, recorded in the
manifest) and never tuned on target expression. Discovery/target disjointness is asserted
from cohort identifiers by ``prs.assert_disjoint_cohorts``-style check in the scripts.

Caveat for every report: discovery tissue is blood/PBMC, target tissue is LCL, so this
validates the scoring machinery, not cross-tissue biology.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np
import pandas as pd

from .io import GenotypeBlock
from .stats_utils import empirical_pvalue

PALINDROMIC = {frozenset("AT"), frozenset("CG")}
REQUIRED_WEIGHT_COLUMNS = ["gene_id", "variant_id", "pvalue", "zscore", "assessed_allele", "other_allele"]


def is_palindromic(a1: str, a2: str) -> bool:
    return frozenset((a1.upper(), a2.upper())) in PALINDROMIC


def align_weights(weights: pd.DataFrame, block: GenotypeBlock) -> tuple[pd.DataFrame, dict]:
    """Align discovery weights to the target block's counted allele.

    Returns (aligned weights with column ``w`` = z in units of the block's counted allele,
    report dict with the number dropped for each reason). Nothing is silently dropped.
    """
    missing = [c for c in REQUIRED_WEIGHT_COLUMNS if c not in weights.columns]
    if missing:
        raise ValueError(f"weights table missing columns: {missing}")
    idx = pd.Series(np.arange(len(block.variant_ids)), index=block.variant_ids)
    idx = idx[~idx.index.duplicated(keep=False)]  # ambiguous IDs are dropped, counted below
    w = weights.copy()
    rep = {"n_input": int(len(w))}
    w["j"] = w["variant_id"].map(idx)
    rep["dropped_not_in_target"] = int(w["j"].isna().sum())
    w = w.dropna(subset=["j"]).copy()
    w["j"] = w["j"].astype(int)
    counted = block.counted_allele[w["j"].to_numpy()]
    other = block.other_allele[w["j"].to_numpy()]
    pal = np.array([is_palindromic(a, b) for a, b in zip(counted, other)])
    rep["dropped_palindromic"] = int(pal.sum())
    w, counted, other = w[~pal].copy(), counted[~pal], other[~pal]
    a = w["assessed_allele"].str.upper().to_numpy()
    same = a == np.char.upper(counted.astype(str))
    flip = a == np.char.upper(other.astype(str))
    rep["dropped_allele_mismatch"] = int((~(same | flip)).sum())
    keep = same | flip
    w = w[keep].copy()
    w["w"] = np.where(same[keep], w["zscore"], -w["zscore"])
    rep["n_aligned"] = int(len(w))
    return w, rep


def clump_greedy(G: np.ndarray, order: np.ndarray, r2_threshold: float) -> list[int]:
    """Greedy LD clumping on columns of G, visiting columns in ``order`` (most significant first)."""
    Gc = np.where(np.isfinite(G), G, np.nanmean(G, axis=0))
    Gc = Gc - Gc.mean(axis=0)
    sd = np.sqrt((Gc ** 2).sum(axis=0))
    kept: list[int] = []
    for j in order:
        if sd[j] == 0:
            continue
        redundant = False
        for i in kept:
            r = float((Gc[:, i] @ Gc[:, j]) / (sd[i] * sd[j]))
            if r * r > r2_threshold:
                redundant = True
                break
        if not redundant:
            kept.append(int(j))
    return kept


@dataclass(frozen=True)
class ExpressionScoreParams:
    p_threshold: float = 1e-4   # [project design choice; fixed a priori, NOT tuned on target]
    r2_clump: float = 0.1       # [project design choice; Choi 2020 conventional range]
    min_snps: int = 1
    cis_window: int = 1_000_000

    def provenance(self) -> dict:
        return asdict(self)


def predict_expression(weights: pd.DataFrame, block: GenotypeBlock, gene_tss: pd.Series,
                       params: ExpressionScoreParams) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Predicted expression for every gene with usable weights.

    Returns (samples x genes predicted-expression frame, per-gene info table, alignment report).
    """
    aligned, report = align_weights(weights, block)
    aligned = aligned[aligned["pvalue"] <= params.p_threshold]
    report["n_after_p_threshold"] = int(len(aligned))
    preds, info = {}, []
    for gene, g in aligned.groupby("gene_id"):
        if gene not in gene_tss.index:
            continue
        j = g["j"].to_numpy()
        pos_ok = np.abs(block.positions[j] - int(gene_tss[gene])) <= params.cis_window
        g, j = g[pos_ok], j[pos_ok]
        if len(j) < params.min_snps:
            continue
        G = block.dosage[:, j].astype(float)
        order = np.argsort(g["pvalue"].to_numpy(), kind="stable")
        kept = clump_greedy(G, order, params.r2_clump)
        if not kept:
            continue
        Gk = G[:, kept]
        Gk = np.where(np.isfinite(Gk), Gk, np.nanmean(Gk, axis=0))
        score = Gk @ g["w"].to_numpy()[kept]
        if np.std(score) == 0:
            continue
        preds[gene] = score
        info.append(dict(gene_id=gene, n_weights=len(j), n_snps_in_score=len(kept)))
    pred = pd.DataFrame(preds, index=block.sample_ids)
    return pred, pd.DataFrame(info), report


def _standardise(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, dtype=float)
    return (M - M.mean(axis=0)) / M.std(axis=0, ddof=1)


def evaluate_expression_score(predicted: pd.DataFrame, measured: pd.DataFrame, n_permutations: int = 1000,
                              n_bootstrap: int = 1000, seed: int = 0) -> dict:
    """Observed performance + uncertainty + permutation null (V6 biological-evidence half).

    Primary statistic (pre-registered): mean per-gene Pearson r across scored genes.
    Also reported: median r, mean r^2, fraction of genes with r>0, per-gene Fisher-z 95% CIs
    (in the returned per-gene table), bootstrap percentile CI of the primary statistic over
    samples, the permutation-null distribution summary, and the one-sided empirical p-value.
    A null result is a valid scientific outcome and is reported as such.
    """
    genes = [g for g in predicted.columns if g in measured.columns]
    if not genes:
        raise ValueError("no genes in common between predicted and measured expression")
    if not predicted.index.equals(measured.index):
        measured = measured.loc[predicted.index]
    P = _standardise(predicted[genes].to_numpy())
    Y = _standardise(measured[genes].to_numpy())
    n = P.shape[0]
    r = (P * Y).sum(axis=0) / (n - 1)
    rng = np.random.default_rng(seed)
    obs = float(r.mean())
    null = np.empty(n_permutations)
    for b in range(n_permutations):
        perm = rng.permutation(n)
        null[b] = float(((P * Y[perm]).sum(axis=0) / (n - 1)).mean())
    boot = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        Pb, Yb = _standardise(P[idx]), _standardise(Y[idx])
        boot[b] = float(((Pb * Yb).sum(axis=0) / (n - 1)).mean())
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    per_gene = pd.DataFrame({"gene_id": genes, "r": r, "r2": r ** 2,
                             "r_ci_low": np.tanh(z - 1.96 * se), "r_ci_high": np.tanh(z + 1.96 * se)})
    emp = empirical_pvalue(obs, null, "greater")
    return {
        "n_samples": int(n), "n_genes_scored": int(len(genes)),
        "primary_statistic": "mean per-gene Pearson r",
        "observed_mean_r": obs, "observed_median_r": float(np.median(r)),
        "observed_mean_r2": float((r ** 2).mean()), "fraction_genes_r_positive": float((r > 0).mean()),
        "bootstrap_ci95_mean_r": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
        "n_permutations": int(n_permutations), "n_bootstrap": int(n_bootstrap), "seed": int(seed),
        "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
        "null_ci95": [float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))],
        "null_max": float(null.max()),
        "empirical_p_one_sided": emp,
        "per_gene": per_gene,
        "null_distribution": null,
    }
