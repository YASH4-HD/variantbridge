"""Capability labels and live validation status (audit section 1).

Every analysis carries one tier:
  TECHNICALLY_POSSIBLE   data/format/licence permit running it
  STATISTICALLY_DEFENSIBLE  design, power and precedent support interpreting it as science
  PORTFOLIO_DEMO_ONLY    runnable and instructive, too underpowered/mismatched to interpret
  SIMULATED              synthetic data; method demonstration only (never validation)

Validation status is read from artifacts/VALIDATION_STATUS.json, which is written ONLY by
the validation scripts after real execution. Absent => "NOT RUN". Nothing here fabricates a
result.
"""
from __future__ import annotations

import json
from pathlib import Path

TECHNICALLY_POSSIBLE = "TECHNICALLY POSSIBLE"
STATISTICALLY_DEFENSIBLE = "STATISTICALLY DEFENSIBLE"
PORTFOLIO_DEMO_ONLY = "PORTFOLIO DEMONSTRATION ONLY"
SIMULATED = "SIMULATED DATA - METHOD DEMONSTRATION ONLY"

NOT_RUN = "NOT RUN"
UNVERIFIED = "UNVERIFIED"
PASSED = "PASSED"
FAILED = "FAILED"

VALIDATION_IDS = {
    "V1": "Population structure: leading PCs recover 1000 Genomes superpopulations",
    "V2": "Sex check concordance",
    "V3": "Null calibration: lambda ~ 1 on permuted phenotypes",
    "V4": "GEUVADIS answer-key replication (EUR373, YRI89) [RELEASE-BLOCKING]",
    "V5": "Statistical consistency + PLINK2 --glm benchmark [RELEASE-BLOCKING]",
    "V6": "Genetic Expression Score: pipeline implementation + biological predictive evidence (reported separately)",
}

CAPABILITIES = {
    "Sample/variant QC": (STATISTICALLY_DEFENSIBLE, "V2"),
    "Relatedness (KING)": (STATISTICALLY_DEFENSIBLE, None),
    "PCA / population structure": (STATISTICALLY_DEFENSIBLE, "V1"),
    "GEUVADIS cis-eQTL association": (STATISTICALLY_DEFENSIBLE, "V4"),
    "Association engine (beta/SE/CI/p single model)": (STATISTICALLY_DEFENSIBLE, "V5"),
    "Null calibration (lambda)": (STATISTICALLY_DEFENSIBLE, "V3"),
    "Genetic Expression Score (eQTLGen -> GEUVADIS)": (STATISTICALLY_DEFENSIBLE, "V6"),
    "Height score on PGP-UK (optional route c)": (PORTFOLIO_DEMO_ONLY, None),
    "Disease-trait PRS": ("NOT POSSIBLE with open data (needs controlled-access target)", None),
    "Built-in demo cohorts / AMR demo": (SIMULATED, None),
    "EvoResist-AI exploratory AMR prediction": (SIMULATED, None),
}


def status_path(root: str | Path = "artifacts") -> Path:
    return Path(root) / "VALIDATION_STATUS.json"


def load_status(root: str | Path = "artifacts") -> dict:
    """Return {id: record}. Missing file or missing id => status NOT RUN."""
    p = status_path(root)
    data = {}
    if p.exists():
        with open(p) as fh:
            data = json.load(fh)
    out = {}
    for vid in VALIDATION_IDS:
        out[vid] = data.get(vid, {"status": NOT_RUN})
    return out


def record_status(vid: str, record: dict, root: str | Path = "artifacts") -> Path:
    if vid not in VALIDATION_IDS:
        raise KeyError(vid)
    if "status" not in record or "executed_at" not in record:
        raise ValueError("record must include 'status' and 'executed_at'")
    p = status_path(root)
    p.parent.mkdir(parents=True, exist_ok=True)
    data = json.loads(p.read_text()) if p.exists() else {}
    data[vid] = record
    p.write_text(json.dumps(data, indent=2, sort_keys=True, default=str))
    return p
