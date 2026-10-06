"""Loading compact validated artifacts produced by the offline scripts.

The Streamlit app never downloads or processes raw genomic data. If an artifact or its
manifest is missing, loaders return (None, reason) and the UI shows NOT RUN."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .manifest import file_sha256, read_manifest, validate_manifest

ARTIFACT_DIR = Path("artifacts")


def load_table_artifact(name: str, root: Path | str = ARTIFACT_DIR) -> tuple[Optional[pd.DataFrame], Optional[dict], str]:
    """Return (dataframe, manifest, message). message is "" on success."""
    root = Path(root)
    path = root / name
    mpath = root / (name.split(".")[0] + ".manifest.json")
    if not path.exists():
        return None, None, f"artifact {name} not found (offline preparation NOT RUN)"
    if not mpath.exists():
        return None, None, f"artifact {name} has no manifest; refusing to display unverifiable data"
    m = read_manifest(mpath)
    missing = validate_manifest(m)
    if missing:
        return None, m, f"manifest incomplete (missing {missing})"
    expected = m.get("output_sha256")
    if expected and expected != file_sha256(path):
        return None, m, f"checksum mismatch for {name}; artifact modified after preparation"
    return pd.read_csv(path), m, ""
