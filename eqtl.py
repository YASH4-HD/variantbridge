"""GEUVADIS-native cis-eQTL engine (audit section 5, stages 5-6; Implementation Spec M6).

Model (source: Lappalainen 2013 SI sections 8 and 11, via audit): linear regression of
PEER-residualised, standard-normal-transformed expression on genotype dosage with covariates
(imputation status + PCs: 1-3 EUR, 1-2 YRI), variants within +/-1 Mb of the TSS, MAF>5% in
either population, autosomes only, populations analysed separately.

Multiple-testing calibration here is DIRECT PER-GENE PERMUTATION (best |t| over cis variants
vs. permuted-phenotype best), then Benjamini-Hochberg across genes. This follows the
FastQTL/tensorQTL "direct" scheme the audit lists as an acceptable equivalent [project design
choice]; it is NOT the exact permutation scheme of the 2013 paper (which used thresholds from
permuting a random gene subset). Differences are documented in docs/METHODS.md.

Not applied here by design: 5e-8 GWAS threshold, LD pruning of test variants, pooled EUR+YRI.
Capability tier: STATISTICALLY DEFENSIBLE (designed purpose of the data); replication of the
published answer key is test V4 and is NOT RUN until executed on the real data.
"""
from __future__ import annotations

from typing import Callable, Dict, Iterable, Optional

import numpy as np
import pandas as pd
from scipy import stats

from .io import GenotypeBlock
from .stats_utils import benjamini_hochberg, genomic_inflation_lambda
from . import association

DEFAULT_WINDOW = 1_000_000  # [source: GEUVADIS SI via audit]
DEFAULT_MAF = 0.05          # [source: GEUVADIS SI via audit]
PCS_BY_POPULATION = {"EUR": 3, "YRI": 2}  # [source: GEUVADIS SI via audit]


def cis_slice(positions: np.ndarray, tss: int, window: int = DEFAULT_WINDOW) -> slice:
    lo = int(np.searchsorted(positions, tss - window, side="left"))
    hi = int(np.searchsorted(positions, tss + window, side="right"))
    return slice(lo, hi)


def build_covariates(sample_ids, pcs: pd.DataFrame, n_pcs: int,
                     imputation_status: Optional[pd.Series] = None) -> pd.DataFrame:
    """Covariate table: PCs 1..n_pcs (+ imputation status when provided). Index = sample_id."""
    pc = pcs.set_index("sample_id").loc[list(sample_ids), [f"PC{i + 1}" for i in range(n_pcs)]]
    if imputation_status is not None:
        pc = pc.assign(imputed=imputation_status.loc[list(sample_ids)].astype(float).to_numpy())
    return pc


def maf_either_population(blocks: Iterable[GenotypeBlock], threshold: float = DEFAULT_MAF) -> np.ndarray:
    """Variant mask: MAF > threshold in ANY of the given (same-variant-set) blocks."""
    keep = None
    for b in blocks:
        af = np.nanmean(b.dosage.astype(float), axis=0) / 2.0
        m = np.minimum(af, 1 - af) > threshold
        keep = m if keep is None else (keep | m)
    return keep


def gene_permutation_test(y: np.ndarray, G: np.ndarray, covariates: Optional[np.ndarray],
                          n_perm: int, rng: np.random.Generator) -> dict:
    """Per-gene permutation test on the best cis variant (phenotype residuals permuted).

    Both y and G are residualised on [1, covariates]; the residualised phenotype is permuted.
    Because the nominal p-value is a monotone function of |r| at fixed df, the best p per
    permutation is the max |r|. Returns the observed best statistics and the empirical gene p.
    """
    n, m = G.shape
    cov = np.empty((n, 0)) if covariates is None else np.asarray(covariates, dtype=float)
    A = np.column_stack([np.ones(n), cov])
    Q, _ = np.linalg.qr(A)
    yr = y - Q @ (Q.T @ y)
    Gr = G - Q @ (Q.T @ G)
    sxx = np.einsum("ij,ij->j", Gr, Gr)
    ok = sxx > 1e-10 * n
    if not ok.any():
        return dict(empirical_p=np.nan, best_idx=-1, n_tested=0)
    Gr = Gr[:, ok]
    sxx = sxx[ok]
    syy = float(yr @ yr)
    r_obs = (Gr.T @ yr) / np.sqrt(sxx * syy)
    best_obs = float(np.max(np.abs(r_obs)))
    best_idx_local = int(np.argmax(np.abs(r_obs)))
    best_idx = int(np.flatnonzero(ok)[best_idx_local])
    if n_perm > 0:
        Yp = np.stack([yr[rng.permutation(n)] for _ in range(n_perm)], axis=1)
        R = (Gr.T @ Yp) / np.sqrt(sxx[:, None] * syy)
        best_perm = np.abs(R).max(axis=0)
        emp = float((1 + np.sum(best_perm >= best_obs)) / (1 + n_perm))
    else:
        emp = np.nan
    return dict(empirical_p=emp, best_idx=best_idx, n_tested=int(ok.sum()))


def run_cis_scan(phenotypes: pd.DataFrame, gene_meta: pd.DataFrame, blocks: Dict[str, GenotypeBlock],
                 covariates: pd.DataFrame, population: str, variant_mask: Optional[Dict[str, np.ndarray]] = None,
                 window: int = DEFAULT_WINDOW, n_perm: int = 0, seed: int = 0,
                 keep_pairs: Optional[set] = None, pair_sink: Optional[Callable[[pd.DataFrame], None]] = None,
                 allow_mean_impute: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the cis scan for every gene.

    phenotypes : samples x genes (index = sample_id, columns = gene_id), already normalised.
    gene_meta  : index gene_id; columns chrom, tss.
    blocks     : chrom -> GenotypeBlock (same sample order as ``phenotypes.index``).
    covariates : samples x covariates (same sample order).
    keep_pairs : optional set of (gene_id, variant_id) whose full statistics are returned.
    pair_sink  : optional callback receiving every gene's full cis table (for streaming to disk).

    Returns (per-gene best-hit table, retained pairs table).
    """
    rng = np.random.default_rng(seed)
    samples = phenotypes.index.to_numpy()
    for c, b in blocks.items():
        if not np.array_equal(b.sample_ids, samples):
            raise ValueError(f"sample order mismatch between phenotypes and genotype block chr{c}")
    if not np.array_equal(covariates.index.to_numpy(), samples):
        raise ValueError("sample order mismatch between phenotypes and covariates")
    cov = covariates.to_numpy(dtype=float)
    best_rows, pair_rows = [], []
    for gene in phenotypes.columns:
        if gene not in gene_meta.index:
            continue
        chrom = str(gene_meta.at[gene, "chrom"]).replace("chr", "")
        if chrom not in blocks:
            continue
        b = blocks[chrom]
        sl = cis_slice(b.positions, int(gene_meta.at[gene, "tss"]), window)
        idx = np.arange(b.positions.shape[0])[sl]
        if variant_mask is not None and chrom in variant_mask:
            idx = idx[variant_mask[chrom][idx]]
        if idx.size == 0:
            continue
        G = b.dosage[:, idx].astype(float)
        if np.isnan(G).any():
            if not allow_mean_impute:
                raise ValueError("missing genotypes in cis window; set allow_mean_impute=True explicitly "
                                 "(recorded deviation) or fix upstream")
            G = np.where(np.isnan(G), np.nanmean(G, axis=0), G)
        y = phenotypes[gene].to_numpy(dtype=float)
        res = association.scan_linear(G, y, cov, variant_ids=b.variant_ids[idx])
        res["gene_id"] = gene
        res["pos"] = b.positions[idx]
        res["tss_distance"] = res["pos"] - int(gene_meta.at[gene, "tss"])
        if pair_sink is not None:
            pair_sink(res)
        if keep_pairs:
            sel = [(gene, v) in keep_pairs for v in res["variant"]]
            if any(sel):
                pair_rows.append(res[sel])
        ok = res[res["status"] == "ok"]
        if ok.empty:
            continue
        top = ok.loc[ok["p"].idxmin()]
        perm = gene_permutation_test(y, G, cov, n_perm, rng) if n_perm > 0 else dict(empirical_p=np.nan)
        best_rows.append(dict(gene_id=gene, chrom=chrom, n_variants=int(len(ok)), best_variant=top["variant"],
                              best_pos=int(top["pos"]), tss_distance=int(top["tss_distance"]),
                              effect=top["effect"], se=top["se"], p_nominal=top["p"],
                              empirical_p=perm["empirical_p"]))
    best = pd.DataFrame(best_rows)
    if not best.empty and n_perm > 0:
        best["q_bh"] = benjamini_hochberg(best["empirical_p"].to_numpy())
    pairs = pd.concat(pair_rows, ignore_index=True) if pair_rows else pd.DataFrame()
    return best, pairs


def permuted_scan_lambda(phenotypes: pd.DataFrame, gene_meta: pd.DataFrame, blocks: Dict[str, GenotypeBlock],
                         covariates: pd.DataFrame, n_permutations: int = 1, seed: int = 0, window: int = DEFAULT_WINDOW,
                         variant_mask=None, max_genes: Optional[int] = None) -> list[float]:
    """V3 helper: genomic inflation lambda of ALL cis gene-variant p-values after permuting
    sample labels of the expression matrix (genotypes/covariates fixed). Expected ~1."""
    rng = np.random.default_rng(seed)
    lambdas = []
    genes = list(phenotypes.columns)[: max_genes or None]
    for _ in range(n_permutations):
        perm = rng.permutation(len(phenotypes))
        ph = pd.DataFrame(phenotypes.to_numpy()[perm], index=phenotypes.index, columns=phenotypes.columns)
        ps: list[np.ndarray] = []
        run_cis_scan(ph[genes], gene_meta, blocks, covariates, "perm", variant_mask, window, 0, 0,
                     pair_sink=lambda d: ps.append(d["p"].to_numpy()))
        lambdas.append(genomic_inflation_lambda(np.concatenate(ps)) if ps else float("nan"))
    return lambdas
