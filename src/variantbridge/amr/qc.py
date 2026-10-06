"""Genome QC (audit step [3]); thresholds come from config with provenance, never constants.

Missing QC metadata is NOT treated as passing. Species identity (k-mer/ANI) is not computed here:
if no `species_ok` column is supplied the check is reported as NOT RUN."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .config import threshold


def _col(df: pd.DataFrame, *names: str):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def genome_qc(meta: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Return meta + boolean gate columns + `qc_pass` + `qc_fail_reasons` (semicolon list)."""
    comp_min = threshold(cfg, "genome_qc.checkm_completeness_min")
    cont_max = threshold(cfg, "genome_qc.checkm_contamination_max")
    lo_mb, hi_mb = threshold(cfg, "genome_qc.genome_length_mb")
    max_contigs = threshold(cfg, "genome_qc.max_contigs")
    lab_re = re.compile(threshold(cfg, "genome_qc.lab_strain_regex"))
    out = meta.copy()
    c_comp = _col(out, "checkm_completeness", "completeness")
    c_cont = _col(out, "checkm_contamination", "contamination")
    c_len = _col(out, "genome_length", "length")
    c_ctg = _col(out, "contigs", "n_contigs")
    c_name = _col(out, "genome_name", "name")

    def num(c):
        return pd.to_numeric(out[c], errors="coerce") if c else pd.Series(np.nan, index=out.index)

    comp, cont, glen, ctg = num(c_comp), num(c_cont), num(c_len), num(c_ctg)
    gates = {
        "missing_checkm": comp.isna() | cont.isna(),
        "low_completeness": comp < comp_min,
        "high_contamination": cont > cont_max,
        "missing_length": glen.isna(),
        "length_out_of_range": (glen / 1e6 < lo_mb) | (glen / 1e6 > hi_mb),
        "too_many_contigs": ctg > max_contigs,
        "lab_or_reference_strain_name": (out[c_name].astype(str).apply(lambda s: bool(lab_re.search(s))) if c_name else pd.Series(False, index=out.index)),
    }
    if "species_ok" in out:
        gates["species_mismatch"] = ~out["species_ok"].astype(bool)
    for k, v in gates.items():
        out[f"qc_{k}"] = v.fillna(False).astype(bool)
    out["qc_fail_reasons"] = [";".join(k for k in gates if out.loc[i, f"qc_{k}"]) for i in out.index]
    out["qc_pass"] = out["qc_fail_reasons"] == ""
    return out


def qc_report(qc: pd.DataFrame) -> dict:
    gate_cols = [c for c in qc.columns if c.startswith("qc_") and c not in ("qc_pass", "qc_fail_reasons")]
    return {"n": int(len(qc)), "n_pass": int(qc["qc_pass"].sum()),
            "per_gate_failures": {c[3:]: int(qc[c].sum()) for c in gate_cols},
            "species_check": "RUN (species_ok supplied)" if "qc_species_mismatch" in qc else "NOT RUN (no species_ok column / ANI not computed)"}
