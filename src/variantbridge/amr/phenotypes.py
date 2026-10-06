"""MIC parsing with EXPLICIT censoring (user-critical requirement).

Raw MIC strings such as '<=0.25', '>4', '>=4', '1' are never turned into ordinary continuous
numbers. Each record keeps, in separate columns:

  mic_raw               original string, untouched
  mic_operator          '', '=', '<=', '<', '>=', '>'
  mic_boundary_mg_l     the number in the string (a BOUNDARY if censored)
  censoring             'exact' | 'left' | 'right' | 'invalid'
  bound_inclusive       whether the boundary value itself is allowed
  log2_mic_exact        log2(MIC) ONLY for exact records; NaN for censored records
  log2_boundary         log2(boundary) for censored records (NaN for exact)
  log2_lower/log2_upper interval the true log2 MIC is known to lie in ([-inf, b] left; [b, inf] right;
                        [v, v] exact)
  ordinal_code          dense rank of (boundary, left<exact<right) for ORDINAL models only

Left-censored: true MIC <= boundary.  Right-censored: true MIC >= boundary.
No column holds boundary values disguised as exact measurements, so a naive linear model on
`log2_mic_exact` would silently drop censored rows (visible), never silently treat them as exact.

Breakpoint-dependent labels are produced only through `assign_binary`, which refuses unless the
breakpoint is VERIFIED in config (see config.require_verified_breakpoints).
"""
from __future__ import annotations

import math
import re
from typing import Optional

import numpy as np
import pandas as pd

from .config import require_verified_breakpoints

_OPS = ["<=", ">=", "<", ">", "="]
_NORMALISE = {"≤": "<=", "≥": ">=", "=<": "<=", "=>": ">="}
_NUM = re.compile(r"^[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?$")

EXACT, LEFT, RIGHT, INVALID = "exact", "left", "right", "invalid"


def _split_operator(s: str) -> tuple[str, str]:
    for k, v in _NORMALISE.items():
        s = s.replace(k, v)
    s = s.strip().replace(" ", "")
    for op in _OPS:
        if s.startswith(op):
            return op, s[len(op):]
    return "", s


def parse_mic(raw, sign: Optional[str] = None, value=None) -> dict:
    """Parse one MIC. `raw` is the string ('<=0.25'); `sign`/`value` are optional separate fields
    (BV-BRC has measurement_sign / measurement_value). If both forms are present and disagree the
    record is INVALID (never silently choosing one)."""
    out = {"mic_raw": None if raw is None or (isinstance(raw, float) and math.isnan(raw)) else str(raw),
           "mic_operator": None, "mic_boundary_mg_l": np.nan, "censoring": INVALID, "bound_inclusive": None,
           "log2_mic_exact": np.nan, "log2_boundary": np.nan, "log2_lower": np.nan, "log2_upper": np.nan,
           "mic_parse_ok": False, "mic_parse_note": ""}
    cand = []
    if out["mic_raw"] not in (None, "", "nan", "None"):
        op, num = _split_operator(out["mic_raw"])
        cand.append((op, num, "raw"))
    if value is not None and not (isinstance(value, float) and math.isnan(value)) and str(value).strip() not in ("", "nan", "None"):
        op2, num2 = _split_operator(str(sign or "") + str(value))
        cand.append((op2, num2, "sign+value"))
    if not cand:
        out["mic_parse_note"] = "no MIC string"
        return out
    parsed = []
    for op, num, src in cand:
        if not _NUM.match(num):
            out["mic_parse_note"] = f"non-numeric MIC {num!r} ({src})"
            return out
        v = float(num)
        if not math.isfinite(v) or v <= 0:
            out["mic_parse_note"] = f"MIC must be positive finite, got {num!r}"
            return out
        parsed.append((("=" if op == "" else op), v))
    if len(set(parsed)) > 1:
        out["mic_parse_note"] = f"conflicting representations {parsed}"
        return out
    op, v = parsed[0]
    out["mic_operator"], out["mic_boundary_mg_l"] = op, v
    lg = math.log2(v)
    if op == "=":
        out.update(censoring=EXACT, bound_inclusive=True, log2_mic_exact=lg, log2_lower=lg, log2_upper=lg)
    elif op in ("<=", "<"):
        out.update(censoring=LEFT, bound_inclusive=(op == "<="), log2_boundary=lg, log2_lower=-np.inf, log2_upper=lg)
    else:
        out.update(censoring=RIGHT, bound_inclusive=(op == ">="), log2_boundary=lg, log2_lower=lg, log2_upper=np.inf)
    out["mic_parse_ok"] = True
    return out


def parse_mic_table(df: pd.DataFrame, raw_col: str = "measurement", sign_col: str = "measurement_sign",
                    value_col: str = "measurement_value", grid_tol_log2: float = 0.10) -> pd.DataFrame:
    """Vector wrapper: appends the explicit censoring columns. Row order preserved."""
    rows = []
    for i in range(len(df)):
        raw = df[raw_col].iloc[i] if raw_col in df else None
        sg = df[sign_col].iloc[i] if sign_col in df else None
        vl = df[value_col].iloc[i] if value_col in df else None
        rows.append(parse_mic(raw, sg, vl))
    parsed = pd.DataFrame(rows, index=df.index)
    out = pd.concat([df.drop(columns=[c for c in parsed.columns if c in df.columns]), parsed], axis=1)
    out["off_grid"] = flag_off_grid(out, grid_tol_log2)
    out["ordinal_code"] = ordinal_codes(out)
    return out


def flag_off_grid(df: pd.DataFrame, tol: float = 0.10) -> pd.Series:
    """True if log2(boundary) is not (within tol) on the two-fold dilution grid. Flag only; never altered."""
    lg = np.log2(df["mic_boundary_mg_l"].astype(float))
    frac = np.abs(lg - np.round(lg))
    return (frac > tol) & df["mic_parse_ok"].astype(bool)


def ordinal_codes(df: pd.DataFrame) -> pd.Series:
    """Dense integer rank of (log2 boundary, left < exact < right). For ordinal models only. Invalid -> -1."""
    rank = {LEFT: -1, EXACT: 0, RIGHT: 1}
    ok = df["mic_parse_ok"].astype(bool)
    keys = [(round(float(b), 6), rank[c]) for b, c in zip(np.log2(df.loc[ok, "mic_boundary_mg_l"].astype(float)), df.loc[ok, "censoring"])]
    levels = {k: i for i, k in enumerate(sorted(set(keys)))}
    codes = pd.Series(-1, index=df.index, dtype=int)
    codes.loc[ok] = [levels[k] for k in keys]
    return codes


def censoring_summary(df: pd.DataFrame) -> dict:
    """Counts needed to judge how censored the data are; reported in every manifest and card."""
    n = len(df)
    c = df["censoring"].value_counts().to_dict() if n else {}
    return {"n": int(n), "exact": int(c.get(EXACT, 0)), "left": int(c.get(LEFT, 0)), "right": int(c.get(RIGHT, 0)),
            "invalid": int(c.get(INVALID, 0)),
            "fraction_censored": float((c.get(LEFT, 0) + c.get(RIGHT, 0)) / n) if n else float("nan"),
            "off_grid": int(df["off_grid"].sum()) if "off_grid" in df else 0}


def mic_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """Value counts by (operator, boundary) as reported; no pooling of censored with exact."""
    g = df[df["mic_parse_ok"]].groupby(["mic_operator", "mic_boundary_mg_l"]).size().rename("n").reset_index()
    return g.sort_values(["mic_boundary_mg_l", "mic_operator"]).reset_index(drop=True)


def _ok(df: pd.DataFrame) -> np.ndarray:
    """Parse-ok mask; derived from `censoring` when the compact artifact omits mic_parse_ok."""
    if "mic_parse_ok" in df:
        return df["mic_parse_ok"].astype(bool).to_numpy()
    return (df["censoring"] != INVALID).to_numpy()


def assign_binary(df: pd.DataFrame, cfg: dict, cohort: str) -> pd.Series:
    """Breakpoint-dependent S/I/R from the censored interval. REFUSES unless the breakpoint is VERIFIED.

    Rows whose interval straddles a breakpoint are 'indeterminate' (never forced to a class)."""
    s_max, r_min = require_verified_breakpoints(cfg, cohort)
    lo, hi = df["log2_lower"].to_numpy(float), df["log2_upper"].to_numpy(float)
    ls, lr = math.log2(s_max), math.log2(r_min)
    lab = np.full(len(df), "indeterminate", dtype=object)
    ok = _ok(df)
    lab[ok & (hi <= ls)] = "susceptible"
    lab[ok & (lo >= lr)] = "resistant"
    lab[ok & (lo > ls) & (hi < lr)] = "intermediate"
    lab[~ok] = "invalid"
    return pd.Series(lab, index=df.index, name="binary_from_verified_breakpoint")


def extreme_contrast(df: pd.DataFrame) -> pd.Series:
    """BREAKPOINT-FREE extreme-class label (0 = 'low', 1 = 'high', NaN = not in either class).

    low  = left-censored records (true MIC <= b) whose bound b lies strictly below every right-censored bound;
    high = right-censored records (true MIC >= b) whose bound b lies strictly above every retained low bound.
    Every retained low isolate therefore has a true MIC strictly below every retained high isolate, with no
    assumption about clinical breakpoints. This is NOT an S/R call and must never be described as one."""
    ok = pd.Series(_ok(df), index=df.index)
    out = pd.Series(np.nan, index=df.index, name="extreme_contrast", dtype=float)
    left = df[ok & (df["censoring"] == LEFT)]
    right = df[ok & (df["censoring"] == RIGHT)]
    if left.empty or right.empty:
        return out
    low = left[left["log2_boundary"] < right["log2_boundary"].min()]
    if low.empty:
        return out
    high = right[right["log2_boundary"] > low["log2_boundary"].max()]
    out.loc[low.index] = 0.0
    out.loc[high.index] = 1.0
    return out
