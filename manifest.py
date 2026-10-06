"""Reproducibility manifests (handoff section 12; audit Implementation Spec part E)."""
from __future__ import annotations

import datetime as _dt
import hashlib
import importlib
import json
import platform
import sys
from pathlib import Path

NOT_RECORDED = "NOT RECORDED"
NOT_APPLICABLE = "NOT APPLICABLE"

REQUIRED_FIELDS = [
    "dataset_name", "source", "citation", "license_access", "access_date", "phenotype_definition",
    "inclusion_exclusion", "qc_thresholds", "covariates", "model_specification", "random_seed",
    "software_versions", "analysis_date", "limitations",
]

_PACKAGES = ["numpy", "scipy", "pandas", "scikit-learn", "statsmodels", "streamlit", "plotly"]


def software_versions(extra: dict | None = None) -> dict:
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for p in _PACKAGES:
        try:
            from importlib import metadata
            out[p] = metadata.version(p)
        except Exception:
            out[p] = "not installed"
    if extra:
        out.update(extra)
    return out


def file_sha256(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def build_manifest(**fields) -> dict:
    """Assemble a manifest. Every required field must be supplied explicitly; if a value truly
    is unknown pass NOT_RECORDED / NOT_APPLICABLE rather than leaving it out (fail loudly)."""
    fields.setdefault("software_versions", software_versions())
    fields.setdefault("analysis_date", _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"))
    missing = [k for k in REQUIRED_FIELDS if k not in fields or fields[k] is None]
    if missing:
        raise ValueError(f"manifest missing required fields: {missing}")
    return dict(fields)


def write_manifest(path, manifest: dict) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True, default=str)


def read_manifest(path) -> dict:
    with open(path) as fh:
        return json.load(fh)


def validate_manifest(m: dict) -> list[str]:
    return [k for k in REQUIRED_FIELDS if k not in m]
