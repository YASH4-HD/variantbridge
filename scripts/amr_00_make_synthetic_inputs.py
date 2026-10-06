#!/usr/bin/env python
"""Write SYNTHETIC EvoResist inputs (records.csv, meta.csv, features.npz + manifest) for plumbing
demonstrations and tests. Everything is simulated; manifests say data_scope=synthetic and every downstream
script propagates that. Never put these files under artifacts/amr/."""
import argparse
from pathlib import Path

import _amr_common as K  # noqa: F401
from variantbridge.amr import synthetic
from variantbridge.amr.features_io import save_feature_matrix
from variantbridge.manifest import build_manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    s = synthetic.make_synthetic(seed=a.seed)
    s["records"].to_csv(out / "records.csv", index=False)
    s["meta"].to_csv(out / "meta.csv", index=False)
    mf = build_manifest(dataset_name="SYNTHETIC", source="synthetic generator", citation="n/a", license_access="n/a", access_date="n/a",
                        phenotype_definition="synthetic log2 MIC with planted effects", inclusion_exclusion="n/a", qc_thresholds="n/a",
                        covariates="n/a", model_specification="n/a", random_seed=a.seed, limitations="SIMULATED - no biological meaning")
    save_feature_matrix(out / "features.npz", s["features"].to_numpy(), list(s["features"].index), list(s["features"].columns),
                        s["feature_meta"], mf, data_scope="synthetic")
    print(f"synthetic inputs written to {out} (data_scope=synthetic)")


if __name__ == "__main__":
    main()
