"""Compact feature-matrix artifact (presence/absence or dosage) produced OFFLINE by external callers
(unitig-caller, panaroo/Roary, variant/AMRFinderPlus point-mutation tables). Format (npz + manifest):
  X (uint8 or float32, samples x features), sample_ids (genome_id strings), feature_ids, feature_types,
  and a feature-metadata CSV (feature_id, feature_type, gene, aa_ref, aa_pos, aa_alt)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..manifest import file_sha256, read_manifest, write_manifest


def save_feature_matrix(path: str | Path, X: np.ndarray, sample_ids, feature_ids, feature_meta: pd.DataFrame,
                        manifest_fields: dict, data_scope: str) -> Path:
    if data_scope not in ("real", "synthetic"):
        raise ValueError("data_scope must be 'real' or 'synthetic'")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    X = np.asarray(X)
    if X.shape != (len(sample_ids), len(feature_ids)):
        raise ValueError("matrix shape does not match ids")
    Xs = X.astype(np.uint8) if set(np.unique(X)).issubset({0, 1}) else X.astype(np.float32)
    np.savez_compressed(path, X=Xs, sample_ids=np.asarray(sample_ids, dtype=str), feature_ids=np.asarray(feature_ids, dtype=str))
    meta_path = path.with_suffix("").with_suffix(".feature_meta.csv.gz")
    feature_meta.to_csv(meta_path, index=False)
    m = dict(manifest_fields)
    m.update(data_scope=data_scope, output_sha256=file_sha256(path), feature_meta_file=meta_path.name,
             feature_meta_sha256=file_sha256(meta_path), n_samples=int(X.shape[0]), n_features=int(X.shape[1]))
    write_manifest(path.with_suffix("").with_suffix(".manifest.json"), m)
    return path


def load_feature_matrix(path: str | Path):
    """Return (X float, sample_ids, feature_ids, feature_meta, manifest). Refuses checksum mismatch / missing manifest."""
    path = Path(path)
    mp = path.with_suffix("").with_suffix(".manifest.json")
    if not mp.exists():
        raise FileNotFoundError(f"{mp.name} missing: refusing to load an unverifiable feature matrix")
    m = read_manifest(mp)
    if m.get("output_sha256") != file_sha256(path):
        raise ValueError("feature matrix checksum mismatch")
    z = np.load(path, allow_pickle=False)
    meta = pd.read_csv(path.parent / m["feature_meta_file"])
    return z["X"].astype(float), list(z["sample_ids"]), list(z["feature_ids"]), meta, m
