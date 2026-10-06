"""Lineage / population-structure helpers.

NOT mash: distances here are computed from the supplied gene presence/absence matrix (Jaccard) because
mash/sketching tools are not part of this build; a mash-based distance can be supplied via
`distance_matrix=`. ST (MLST) labels are the lineage unit for fixed effects and for lineage-blocked CV."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import threshold


def lineage_groups(st: pd.Series, min_group: int, max_groups: int) -> pd.Series:
    """ST label -> analysis group. Missing ST => 'unknown'. Small STs pooled to 'other_small';
    only the `max_groups` largest retained groups keep their own label (rest -> 'other_small')."""
    s = st.astype(object).where(st.notna(), "unknown").astype(str)
    counts = s.value_counts()
    keep = [g for g in counts.index if counts[g] >= min_group][:max_groups]
    return s.where(s.isin(keep), "other_small")


def lineage_design(st: pd.Series, cfg: dict) -> tuple[np.ndarray, list[str]]:
    """Fixed-effect dummy matrix (reference = largest group). Returns (matrix, column names)."""
    grp = lineage_groups(st, threshold(cfg, "association.min_isolates_per_lineage_group"), threshold(cfg, "association.max_lineage_groups"))
    levels = grp.value_counts().index.tolist()
    ref, rest = levels[0], levels[1:]
    mat = np.column_stack([(grp == g).to_numpy(float) for g in rest]) if rest else np.zeros((len(grp), 0))
    return mat, [f"lineage[{g}]" for g in rest]


def jaccard_distance(pa: np.ndarray) -> np.ndarray:
    P = (np.asarray(pa) > 0).astype(float)
    inter = P @ P.T
    sizes = P.sum(1)
    union = sizes[:, None] + sizes[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        d = 1.0 - np.where(union > 0, inter / union, 1.0)
    np.fill_diagonal(d, 0.0)
    return d


def pcoa(dist: np.ndarray, n_components: int = 5) -> np.ndarray:
    """Classical MDS (deterministic sign: largest-|loading| positive)."""
    n = dist.shape[0]
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (dist ** 2) @ J
    w, v = np.linalg.eigh(B)
    order = np.argsort(w)[::-1][:n_components]
    w, v = np.clip(w[order], 0, None), v[:, order]
    for k in range(v.shape[1]):
        if v[np.argmax(np.abs(v[:, k])), k] < 0:
            v[:, k] = -v[:, k]
    return v * np.sqrt(w)


def lineage_structure_report(st: pd.Series, y_boundary_or_rank: pd.Series | None = None) -> dict:
    s = st.astype(object).where(st.notna(), "unknown").astype(str)
    c = s.value_counts()
    return {"n_isolates": int(len(s)), "n_lineages": int(len(c)), "largest_lineage": str(c.index[0]) if len(c) else None,
            "largest_lineage_fraction": float(c.iloc[0] / len(s)) if len(c) else float("nan"),
            "n_singletons": int((c == 1).sum())}
