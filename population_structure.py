"""Genotype PCA, simple LD pruning, and region exclusion.

Production PCA on real cohorts is run offline with PLINK2 (``--pca``) or flashpca
(scripts/02_prepare_1000g_pca.py). ``genotype_pca`` here is the reference implementation of
the standard Patterson-style normalised PCA used by the app and unit tests.

Long-range LD regions (Price et al. 2008) are NOT hard-coded: the region table must be
supplied from a verified source (see docs/DATASETS.md). Omitting it is allowed only as an
explicit, manifest-recorded deviation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def genotype_pca(G, n_components: int = 10) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """PCA of additive genotypes. Returns (PCs n x k, variance_ratio k, used_variant_mask).

    Columns are centred by 2p and scaled by sqrt(2p(1-p)); missing values are set to the
    mean (i.e. 0 after centring); monomorphic variants are dropped. Signs are fixed so the
    largest-magnitude sample loading of each PC is positive, making output deterministic.
    """
    G = np.asarray(G, dtype=float)
    n, m = G.shape
    p = np.nanmean(G, axis=0) / 2.0
    keep = np.isfinite(p) & (p > 0) & (p < 1)
    Gk = G[:, keep]
    pk = p[keep]
    Z = (Gk - 2 * pk) / np.sqrt(2 * pk * (1 - pk))
    Z = np.nan_to_num(Z, nan=0.0)
    k = int(min(n_components, n - 1, Z.shape[1]))
    # n x n Gram matrix is cheaper than the SVD of n x m for m >> n
    gram = Z @ Z.T / Z.shape[1]
    evals, evecs = np.linalg.eigh(gram)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    pcs = evecs[:, :k] * np.sqrt(np.maximum(evals[:k], 0) * Z.shape[1])
    for j in range(k):
        if pcs[np.argmax(np.abs(pcs[:, j])), j] < 0:
            pcs[:, j] *= -1
    total = float(np.sum(np.maximum(evals, 0)))
    ratio = np.maximum(evals[:k], 0) / total if total > 0 else np.zeros(k)
    return pcs, ratio, keep


def pcs_frame(pcs: np.ndarray, sample_ids) -> pd.DataFrame:
    d = pd.DataFrame(pcs, columns=[f"PC{i + 1}" for i in range(pcs.shape[1])])
    d.insert(0, "sample_id", list(sample_ids))
    return d


def ld_prune_greedy(G, window_variants: int = 50, r2_threshold: float = 0.2) -> np.ndarray:
    """Greedy forward LD pruning on an ordered variant matrix (reference implementation).

    Keeps a variant unless r^2 with an already-kept variant within the preceding
    ``window_variants`` variants exceeds ``r2_threshold``. NOT identical to PLINK's
    ``--indep-pairwise`` (which also steps/removes by MAF); production pruning uses PLINK2.
    Parameter values are [project design choice] within the conventional range.
    """
    G = np.asarray(G, dtype=float)
    Gc = np.where(np.isfinite(G), G, np.nanmean(G, axis=0))
    Gc = Gc - Gc.mean(axis=0)
    sd = np.sqrt((Gc ** 2).sum(axis=0))
    keep: list[int] = []
    for j in range(G.shape[1]):
        if sd[j] == 0:
            continue
        redundant = False
        for i in reversed(keep):
            if j - i > window_variants:
                break
            r = float((Gc[:, i] @ Gc[:, j]) / (sd[i] * sd[j]))
            if r * r > r2_threshold:
                redundant = True
                break
        if not redundant:
            keep.append(j)
    mask = np.zeros(G.shape[1], dtype=bool)
    mask[keep] = True
    return mask


def regions_to_mask(chrom, pos, regions: pd.DataFrame) -> np.ndarray:
    """True for variants INSIDE any region; ``regions`` needs columns chrom,start,end (1-based inclusive)."""
    chrom = np.asarray([str(c).replace("chr", "") for c in chrom])
    pos = np.asarray(pos)
    inside = np.zeros(len(pos), dtype=bool)
    for r in regions.itertuples(index=False):
        c = str(r.chrom).replace("chr", "")
        inside |= (chrom == c) & (pos >= r.start) & (pos <= r.end)
    return inside


def feature_pca_v01(X, n_components: int = 5):
    """v0.1 PCA (StandardScaler + sklearn PCA), preserved UNCHANGED for the EvoResist-AI branch.

    Used for binary pathogen gene-presence/variant features, where the diploid 2p(1-p) scaling of
    ``genotype_pca`` does not apply. Not a validated microbial population-structure method.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    X = np.asarray(X, dtype=float)
    X = np.where(np.isfinite(X), X, np.nanmean(X, axis=0))
    X = np.nan_to_num(X, nan=0.0)
    p = PCA(n_components=min(n_components, X.shape[1], X.shape[0] - 1))
    a = p.fit_transform(StandardScaler().fit_transform(X))
    return a, p.explained_variance_ratio_
