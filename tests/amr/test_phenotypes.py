import math

import numpy as np
import pandas as pd
import pytest

from variantbridge.amr import phenotypes as P
from variantbridge.amr.config import BreakpointUnverified


@pytest.mark.parametrize("raw,op,cens,incl,lo,hi", [
    ("<=0.25", "<=", "left", True, -np.inf, -2.0),
    ("<0.25", "<", "left", False, -np.inf, -2.0),
    (">4", ">", "right", False, 2.0, np.inf),
    (">=4", ">=", "right", True, 2.0, np.inf),
    ("1", "=", "exact", True, 0.0, 0.0),
    ("=2", "=", "exact", True, 1.0, 1.0),
    ("≤0.25", "<=", "left", True, -np.inf, -2.0),
    ("≥4", ">=", "right", True, 2.0, np.inf),
    (" <= 0.25 ", "<=", "left", True, -np.inf, -2.0),
])
def test_parse_cases(raw, op, cens, incl, lo, hi):
    d = P.parse_mic(raw)
    assert d["mic_parse_ok"] and d["mic_operator"] == op and d["censoring"] == cens and d["bound_inclusive"] == incl
    assert d["log2_lower"] == lo and d["log2_upper"] == hi
    assert d["mic_raw"] == raw


@pytest.mark.parametrize("raw", ["abc", "", None, "<=-1", "0", "<=", "1/2", "nan", "inf", "0,25"])
def test_invalid_is_flagged_never_guessed(raw):
    d = P.parse_mic(raw)
    assert not d["mic_parse_ok"] and d["censoring"] == "invalid" and np.isnan(d["log2_mic_exact"])


def test_censored_never_becomes_exact_continuous():
    """CRITICAL: boundary values must not appear in the exact-measurement column."""
    for raw in ("<=0.25", ">4", ">=4", "<0.03"):
        d = P.parse_mic(raw)
        assert np.isnan(d["log2_mic_exact"]), raw
        assert not np.isnan(d["log2_boundary"])
    assert P.parse_mic("1")["log2_mic_exact"] == 0.0 and np.isnan(P.parse_mic("1")["log2_boundary"])


def test_sign_and_value_fields_and_conflicts():
    assert P.parse_mic("<=0.25", "<=", "0.25")["censoring"] == "left"
    assert P.parse_mic(None, "<=", 0.25)["censoring"] == "left"
    bad = P.parse_mic("<=0.25", ">", "0.25")
    assert bad["censoring"] == "invalid" and "conflicting" in bad["mic_parse_note"]


def test_table_preserves_raw_and_columns():
    df = pd.DataFrame({"measurement": ["<=0.25", "1", ">4", "bogus"], "genome_id": list("abcd")})
    out = P.parse_mic_table(df)
    assert out["mic_raw"].tolist() == ["<=0.25", "1", ">4", "bogus"]
    assert out["censoring"].tolist() == ["left", "exact", "right", "invalid"]
    for c in ("mic_operator", "mic_boundary_mg_l", "log2_mic_exact", "log2_boundary", "log2_lower", "log2_upper", "ordinal_code", "off_grid"):
        assert c in out
    assert out["ordinal_code"].tolist()[3] == -1


def test_ordinal_order_left_exact_right():
    df = P.parse_mic_table(pd.DataFrame({"measurement": [">4", "<=0.25", "4", "0.25", "1", "2"]}))
    order = df.sort_values("ordinal_code")["mic_raw"].tolist()
    assert order == ["<=0.25", "0.25", "1", "2", "4", ">4"]


def test_off_grid_flag_only():
    df = P.parse_mic_table(pd.DataFrame({"measurement": ["0.03", "0.12", "0.25", "3", "0.3"]}))
    assert df["off_grid"].tolist() == [False, False, False, True, True]
    assert df.loc[3, "mic_boundary_mg_l"] == 3.0     # never altered


def test_censoring_summary_and_distribution_do_not_pool():
    df = P.parse_mic_table(pd.DataFrame({"measurement": ["<=0.25"] * 3 + ["0.25"] * 2 + [">4"] + ["x"]}))
    s = P.censoring_summary(df)
    assert (s["left"], s["exact"], s["right"], s["invalid"], s["n"]) == (3, 2, 1, 1, 7)
    assert s["fraction_censored"] == pytest.approx(4 / 7)
    dist = P.mic_distribution(df)
    assert len(dist) == 3 and set(dist["mic_operator"]) == {"<=", "=", ">"}


def test_binary_refuses_without_verified_breakpoint(cfg):
    df = P.parse_mic_table(pd.DataFrame({"measurement": ["<=0.25", ">4"]}))
    with pytest.raises(BreakpointUnverified):
        P.assign_binary(df, cfg, "discovery")
    with pytest.raises(BreakpointUnverified):
        P.assign_binary(df, cfg, "replication")


def test_status_verified_but_empty_fields_still_refuses(cfg):
    import copy
    c = copy.deepcopy(cfg)
    c["breakpoints"]["ciprofloxacin_ecoli_discovery"]["status"] = "VERIFIED"      # numbers/provenance still null
    with pytest.raises(BreakpointUnverified):
        P.assign_binary(P.parse_mic_table(pd.DataFrame({"measurement": ["1"]})), c, "discovery")


def test_binary_interval_logic_with_test_breakpoints(verified_cfg):
    """TEST-ONLY breakpoints (S<=0.5, R>=2): intervals that straddle a breakpoint are indeterminate."""
    raws = ["<=0.25", "<=0.5", "<=1", "0.5", "1", "2", ">=2", ">4", ">=1", "bogus"]
    df = P.parse_mic_table(pd.DataFrame({"measurement": raws}))
    lab = P.assign_binary(df, verified_cfg, "discovery").tolist()
    assert lab == ["susceptible", "susceptible", "indeterminate", "susceptible", "intermediate", "resistant", "resistant", "resistant", "indeterminate", "invalid"]


def test_extreme_contrast_strict_separation_and_no_breakpoint_needed():
    df = P.parse_mic_table(pd.DataFrame({"measurement": ["<=0.25"] * 3 + ["1", "2"] + [">4"] * 2}))
    ec = P.extreme_contrast(df)
    assert ec.tolist()[:3] == [0.0] * 3 and np.isnan(ec[3]) and np.isnan(ec[4]) and ec.tolist()[5:] == [1.0, 1.0]
    # overlapping censoring bounds can never be placed in opposite classes
    df2 = P.parse_mic_table(pd.DataFrame({"measurement": ["<=4", ">=2", "<=0.25", ">4"]}))
    e2 = P.extreme_contrast(df2)
    lows = df2[e2 == 0.0]; highs = df2[e2 == 1.0]
    assert (lows["log2_upper"].max() < highs["log2_lower"].min()) if len(lows) and len(highs) else True
    assert e2.iloc[0] != 0.0 or e2.iloc[1] != 1.0   # '<=4' and '>=2' overlap: not both retained
