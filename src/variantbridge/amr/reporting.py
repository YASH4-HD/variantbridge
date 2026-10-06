"""Provenance manifest (handoff section 12 fields) and methods report for AMR runs."""
from __future__ import annotations

import json
from pathlib import Path

from ..manifest import NOT_APPLICABLE, NOT_RECORDED, build_manifest, write_manifest
from .config import breakpoint_block, threshold_provenance_table


def amr_manifest(cfg: dict, *, phenotype_definition: str, inclusion_exclusion: dict, covariates, model_specification: dict,
                 random_seed, limitations, data_scope: str, extra: dict | None = None, **kw) -> dict:
    """Wraps core build_manifest; adds data_scope, cohort/threshold provenance and breakpoint status."""
    if data_scope not in ("real", "synthetic"):
        raise ValueError("data_scope must be 'real' or 'synthetic'")
    random_seed = NOT_APPLICABLE if random_seed is None else random_seed
    m = build_manifest(dataset_name=cfg["dataset"]["name"], source=cfg["dataset"]["source_api"], citation=cfg["dataset"]["citation"],
                       license_access=cfg["dataset"]["license"], access_date=kw.pop("access_date", NOT_RECORDED),
                       phenotype_definition=phenotype_definition, inclusion_exclusion=inclusion_exclusion,
                       qc_thresholds=threshold_provenance_table(cfg), covariates=covariates, model_specification=model_specification,
                       random_seed=random_seed, limitations=limitations, **kw)
    m["data_scope"] = data_scope
    m["breakpoint_status"] = {c: breakpoint_block(cfg, c).get("status") for c in ("discovery", "replication")}
    m["cohort_tier_mapping"] = cfg["cohorts"]["pmid_to_tier"]
    if extra:
        m.update(extra)
    return m


def methods_report(cfg: dict) -> str:
    """Plain-text methods summary that prints the provenance class of EVERY threshold (addendum item 1)."""
    L = ["EvoResist-AI methods report (configuration view)", "=" * 48, "",
         "Thresholds (config.yaml). None is immutable; provenance class shown:"]
    for r in threshold_provenance_table(cfg):
        L.append(f"  - {r['threshold']} = {r['value']}  [{r['provenance']}; {r['date_or_citation']}]")
    L += ["", "Breakpoints:"]
    for c in ("discovery", "replication"):
        b = breakpoint_block(cfg, c)
        L.append(f"  - {c}: {b['status']} (standard stated in source: {b.get('standard_stated_in_source')})")
    L += ["", "Primary phenotype model: " + cfg["phenotype"]["primary_model"] + " (censoring explicit; no naive linear model on boundary values)."]
    return "\n".join(L)
