"""Plotly figure builders (no Streamlit dependency)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .stats_utils import qq_points


def qq_figure(p, title="QQ plot"):
    e, o = qq_points(p)
    f = go.Figure([go.Scatter(x=e, y=o, mode="markers", name="Observed")])
    L = float(max(e.max(), o.max())) if len(e) else 1.0
    f.add_trace(go.Scatter(x=[0, L], y=[0, L], mode="lines", name="Expected"))
    f.update_layout(title=title, xaxis_title="Expected -log10(p)", yaxis_title="Observed -log10(p)")
    return f


def manhattan_figure(df: pd.DataFrame, title="Association overview"):
    return px.scatter(df, x="pos", y="neglog10p", color=df["chrom"].astype(str), hover_name="variant", title=title)


def pca_figure(pcs: pd.DataFrame, color: pd.Series, ratio=None, title="PCA"):
    d = pcs.copy()
    d["label"] = color.astype(str).to_numpy()
    sub = ""
    if ratio is not None and len(ratio) > 1:
        sub = f" - PC1 {ratio[0] * 100:.1f}% - PC2 {ratio[1] * 100:.1f}%"
    return px.scatter(d, x="PC1", y="PC2", color="label", hover_name="sample_id", title=title + sub)
