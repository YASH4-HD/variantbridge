"""Shared helpers for the EvoResist-AI offline scripts (never imported by Streamlit)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from variantbridge.amr import config as C  # noqa: E402
from variantbridge.manifest import file_sha256, read_manifest, write_manifest  # noqa: E402

DEFAULT_ART = ROOT / "artifacts" / "amr"


def write_table(df: pd.DataFrame, art: Path, name: str, manifest: dict) -> Path:
    art.mkdir(parents=True, exist_ok=True)
    p = art / name
    df.to_csv(p, index=False)
    manifest = dict(manifest)
    manifest["output_file"] = name
    manifest["output_sha256"] = file_sha256(p)
    write_manifest(art / (name.split(".")[0] + ".manifest.json"), manifest)
    return p


def read_verified(art: Path, name: str):
    p = art / name
    mp = art / (name.split(".")[0] + ".manifest.json")
    if not p.exists() or not mp.exists():
        raise SystemExit(f"NOT RUN upstream: {p} or its manifest is missing")
    m = read_manifest(mp)
    if m.get("output_sha256") != file_sha256(p):
        raise SystemExit(f"checksum mismatch for {name}: artifact modified after creation")
    return pd.read_csv(p), m
