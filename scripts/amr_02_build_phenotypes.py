#!/usr/bin/env python
"""Build the compact phenotype artifact: explicit MIC censoring, cohort tiers, duplicate exclusion, genome QC.

Real data:      python scripts/amr_02_build_phenotypes.py --raw-dir data/raw/bvbrc --meta-csv data/raw/bvbrc/genome_metadata.csv
Plumbing demo:  python scripts/amr_02_build_phenotypes.py --records-csv R.csv --meta-csv M.csv --data-scope synthetic --artifact-dir OUT

Outputs (artifact dir): amr_phenotypes.csv.gz (+manifest), amr_exclusions.csv.gz. Raw MIC, parsed boundary, censoring
direction and transformed values are separate columns. Breakpoints are never applied here."""
import argparse
import sys
from pathlib import Path

import pandas as pd

import _amr_common as K
from variantbridge.amr import cohorts, data_io as D, phenotypes as P, qc as Q
from variantbridge.amr import reporting, validation as V

ADDENDUM_DISCOVERY_COUNTS = {("<=", 0.25): 764, ("=", 1.0): 253, ("=", 2.0): 23, ("=", 4.0): 14, (">", 4.0): 1375}   # addendum 1.1 (queried 2026-10-05)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--raw-dir", default=None)
    ap.add_argument("--records-csv", default=None)
    ap.add_argument("--meta-csv", required=True)
    ap.add_argument("--data-scope", choices=["real", "synthetic"], default=None, help="required with --records-csv (no default: never mislabel)")
    ap.add_argument("--artifact-dir", default=str(K.DEFAULT_ART))
    a = ap.parse_args()
    cfg = K.C.load_config(a.config)
    art = Path(a.artifact_dir)
    if a.records_csv:
        if a.data_scope is None:
            sys.exit("--data-scope is required with --records-csv")
        scope, rec_in = a.data_scope, pd.read_csv(a.records_csv, dtype={"genome_id": str, "study_pmid": str})
    else:
        scope = "real"
        try:
            rec_in = D.load_cached_records(D.fetch_genome_amr(cfg, a.raw_dir or "data/raw/bvbrc", None, False))
        except D.DataNotAvailable as e:
            sys.exit(f"NOT RUN: {e}")
    rec = D.fetch_phenotypes(cfg, records=rec_in)
    meta = pd.read_csv(a.meta_csv, dtype={"genome_id": str})
    kept, excl = cohorts.resolve_cohorts(rec, cfg)
    qcdf = Q.genome_qc(meta, cfg)
    bad = qcdf.loc[~qcdf["qc_pass"], ["genome_id", "qc_fail_reasons"]]
    missing_meta = kept.loc[~kept["genome_id"].isin(qcdf["genome_id"]), "genome_id"]
    qc_excl = kept[kept["genome_id"].isin(bad["genome_id"])][["genome_id", "study_pmid", "cohort_tier"]].merge(bad, on="genome_id").rename(columns={"cohort_tier": "record_tier", "qc_fail_reasons": "reason"})
    qc_excl["reason"] = "qc:" + qc_excl["reason"]
    nm = kept[kept["genome_id"].isin(missing_meta)][["genome_id", "study_pmid", "cohort_tier"]].rename(columns={"cohort_tier": "record_tier"})
    nm["reason"] = "qc:no_genome_metadata"
    drop_ids = set(bad["genome_id"]) | set(missing_meta)
    if "mlst" not in qcdf:
        qcdf["mlst"] = None
    kept = kept[~kept["genome_id"].isin(drop_ids)].merge(qcdf[["genome_id", "mlst"]], on="genome_id", how="left")
    excl = pd.concat([excl, qc_excl, nm], ignore_index=True)
    cols = ["genome_id", "study_pmid", "cohort_tier", "mic_raw", "mic_operator", "mic_boundary_mg_l", "mic_parse_ok", "censoring", "bound_inclusive",
            "log2_mic_exact", "log2_boundary", "log2_lower", "log2_upper", "ordinal_code", "off_grid", "laboratory_typing_method",
            "resistant_phenotype", "testing_standard", "mlst"]
    out = kept[[c for c in cols if c in kept.columns]]
    per_tier = {t: P.censoring_summary(out[out["cohort_tier"] == t]) for t in K.C.TIERS}
    m = reporting.amr_manifest(
        cfg, phenotype_definition="log2 MIC with explicit censoring (left/right/exact); R/S calls NOT used (breakpoints unverified)",
        inclusion_exclusion={"exclusion_reason_counts": excl["reason"].value_counts().to_dict(), "tier_counts": cohorts.tier_counts(out),
                             "overlap_matrix_pmids": cohorts.overlap_matrix(cohorts.normalise_ids(rec_in)).to_dict()},
        covariates="ST (MLST) lineage groups", model_specification={"stage": "phenotype build", "censoring_summary_per_tier": per_tier,
                                                                    "qc": Q.qc_report(qcdf)},
        random_seed=None, limitations=["BV-BRC field names UNVERIFIED vs live API", "breakpoint for discovery cohort UNVERIFIED", "species check NOT RUN unless species_ok supplied"],
        data_scope=scope)
    K.write_table(out, art, "amr_phenotypes.csv.gz", m)
    K.write_table(excl, art, "amr_exclusions.csv.gz", dict(m, note="exclusion table"))
    print(f"phenotypes: {len(out)} genomes; tiers {cohorts.tier_counts(out)}; exclusions {len(excl)}; scope={scope}")
    d = out[out["cohort_tier"] == "discovery"]
    if scope == "real" and len(d):
        dist = {(r.mic_operator, float(r.mic_boundary_mg_l)): int(r.n) for r in P.mic_distribution(d).itertuples()}
        diffs = {str(k): {"expected": v, "observed": dist.get(k, 0)} for k, v in ADDENDUM_DISCOVERY_COUNTS.items()}
        ok = all(v["expected"] == v["observed"] for v in diffs.values())
        V.record_status("A5", {"status": V.PASSED if ok else V.FAILED, "executed_at": V.now(), "data_scope": "real", "detail": diffs,
                               "note": "informational reconciliation with addendum 1.1 counts; mismatch is reported, not auto-corrected"}, art)


if __name__ == "__main__":
    main()
