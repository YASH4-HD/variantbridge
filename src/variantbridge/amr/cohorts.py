"""Cohort tiers, duplicate exclusion and cohort-separation guards (addendum Part 1).

Tier = f(PMID). A genome that appears in more than one tier is assigned to the HIGHEST tier only
(discovery > replication > exploratory); its records from lower-tier studies are dropped. The
38219757 / 28720578 overlap is resolved by keeping the superset study's record. Remaining
duplicate records of a genome are collapsed only if identical; conflicting MIC records exclude the
genome (never averaged, never picked arbitrarily).
"""
from __future__ import annotations

import pandas as pd

from .config import TIERS, pmid_to_tier


class CohortSeparationError(RuntimeError):
    pass


def normalise_ids(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["genome_id"] = out["genome_id"].astype(str).str.strip()
    pm = out["study_pmid"] if "study_pmid" in out else out.get("pmid")
    if pm is None:
        raise KeyError("phenotype records need a study_pmid (or pmid) column")
    out["study_pmid"] = pm.astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    return out


def assign_record_tier(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    out = normalise_ids(df)
    out["record_tier"] = [pmid_to_tier(cfg, p) for p in out["study_pmid"]]
    return out


def resolve_cohorts(records: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (kept, exclusions). `kept` has one MIC row per genome and a `cohort_tier` column.

    exclusions columns: genome_id, study_pmid, record_tier, reason. Reasons are exhaustive and
    deterministic: superset_duplicate, lower_tier_duplicate_genome, non_mic_method, invalid_mic,
    exact_duplicate_record, conflicting_duplicate_mic."""
    rec = assign_record_tier(records, cfg)
    order = {t: i for i, t in enumerate(TIERS)}
    excl: list[pd.DataFrame] = []

    def drop(mask, reason):
        nonlocal rec
        if mask.any():
            e = rec.loc[mask, ["genome_id", "study_pmid", "record_tier"]].copy()
            e["reason"] = reason
            excl.append(e)
            rec = rec.loc[~mask]

    # 1. superset study preference (e.g. keep 38219757, drop 28720578 for shared genomes)
    for rule in cfg["cohorts"].get("superset_preference", []):
        keep_g = set(rec.loc[rec["study_pmid"] == rule["keep_pmid"], "genome_id"])
        drop((rec["study_pmid"] == rule["drop_pmid"]) & rec["genome_id"].isin(keep_g), "superset_duplicate")
    # 2. genome in >1 tier -> highest tier only (computed over ALL remaining records, any method)
    best = rec.groupby("genome_id")["record_tier"].agg(lambda s: min(s, key=order.get))
    rec = rec.assign(cohort_tier=rec["genome_id"].map(best))
    drop(rec["record_tier"] != rec["cohort_tier"], "lower_tier_duplicate_genome")
    # 3. primary analyses use quantitative MIC records only
    method = rec.get("laboratory_typing_method", pd.Series("MIC", index=rec.index)).astype(str)
    drop(method.str.upper() != cfg["phenotype"]["primary_method"].upper(), "non_mic_method")
    if "mic_parse_ok" in rec:
        drop(~rec["mic_parse_ok"].astype(bool), "invalid_mic")
    # 4. duplicate records within a genome
    if len(rec):
        key = rec["mic_operator"].astype(str) + "|" + rec["mic_boundary_mg_l"].astype(str)
        n_distinct = key.groupby(rec["genome_id"]).transform("nunique")
        drop(n_distinct > 1, "conflicting_duplicate_mic")
        dup = rec.duplicated(subset=["genome_id"], keep="first")
        drop(dup, "exact_duplicate_record")
    exclusions = pd.concat(excl, ignore_index=True) if excl else pd.DataFrame(columns=["genome_id", "study_pmid", "record_tier", "reason"])
    kept = rec.sort_values("genome_id").reset_index(drop=True)
    assert_cohort_separation(kept, cfg)
    return kept, exclusions


def assert_cohort_separation(df: pd.DataFrame, cfg: dict | None = None) -> None:
    """Hard guard: a genome in exactly one tier, one row; discovery/replication contain only their PMIDs."""
    if df["genome_id"].duplicated().any():
        raise CohortSeparationError(f"genome_id appears more than once: {df.loc[df['genome_id'].duplicated(), 'genome_id'].head().tolist()}")
    t = df.groupby("genome_id")["cohort_tier"].nunique()
    if (t > 1).any():
        raise CohortSeparationError("genome assigned to more than one tier")
    if cfg is not None:
        for pm, tier in cfg["cohorts"]["pmid_to_tier"].items():
            bad = df[(df["study_pmid"] == pm) & (df["cohort_tier"] != tier)]
            if len(bad):
                raise CohortSeparationError(f"PMID {pm} records outside tier {tier}")
        non_mapped = df[df["cohort_tier"].isin(["discovery", "replication"]) & ~df["study_pmid"].isin(cfg["cohorts"]["pmid_to_tier"].keys())]
        if len(non_mapped):
            raise CohortSeparationError("non-registered PMID in discovery/replication tier")


def assert_tier(df: pd.DataFrame, allowed: str, who: str) -> None:
    """Used by association (discovery only) and prediction training (discovery) etc."""
    tiers = set(df["cohort_tier"].unique())
    if tiers != {allowed}:
        raise CohortSeparationError(f"{who} accepts only tier {allowed!r}; got {sorted(tiers)}")


def overlap_matrix(records: pd.DataFrame) -> pd.DataFrame:
    """Shared genome_id counts between PMIDs (addendum 1.2) - a QC report, computed from data."""
    sets = {p: set(g["genome_id"]) for p, g in records.groupby("study_pmid")}
    pm = sorted(sets)
    return pd.DataFrame([[len(sets[a] & sets[b]) if a != b else len(sets[a]) for b in pm] for a in pm], index=pm, columns=pm)


def tier_counts(kept: pd.DataFrame) -> dict:
    return kept["cohort_tier"].value_counts().reindex(list(TIERS), fill_value=0).to_dict()
