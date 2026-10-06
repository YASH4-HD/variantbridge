import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from variantbridge.amr import cohorts, config, qc, reporting


# ------------------------------------------------------------------ config / provenance
def test_every_threshold_has_provenance(cfg):
    rows = config.threshold_provenance_table(cfg)
    assert len(rows) >= 10
    assert all(r["provenance"] in ("source_derived", "pre_registered") for r in rows)
    assert all(r["date_or_citation"] and r["rationale"] for r in rows)


@pytest.mark.parametrize("mutate,msg", [
    (lambda c: c["thresholds"]["genome_qc"]["checkm_completeness_min"].pop("provenance"), "provenance"),
    (lambda c: c["thresholds"]["genome_qc"]["checkm_completeness_min"].update(provenance="vibes"), "provenance"),
    (lambda c: c["thresholds"]["genome_qc"]["checkm_completeness_min"].update(provenance="source_derived"), "citation"),
    (lambda c: c["thresholds"]["genome_qc"]["checkm_completeness_min"].pop("date"), "date"),
    (lambda c: c["thresholds"]["genome_qc"]["checkm_completeness_min"].pop("value"), "value"),
])
def test_invalid_threshold_blocks_load(cfg, tmp_path, mutate, msg):
    c = copy.deepcopy(cfg); mutate(c)
    p = tmp_path / "c.yaml"; p.write_text(yaml.safe_dump(c))
    with pytest.raises(config.ConfigError, match=msg):
        config.load_config(p)


def test_methods_report_prints_every_provenance_class(cfg):
    txt = reporting.methods_report(cfg)
    for r in config.threshold_provenance_table(cfg):
        assert r["threshold"] in txt
    assert txt.count("[pre_registered") >= 10 and "UNVERIFIED" in txt and "interval_censored_gaussian" in txt


def test_shipped_config_has_unverified_breakpoints_and_no_inferred_numbers(cfg):
    for tier in ("discovery", "replication"):
        b = config.breakpoint_block(cfg, tier)
        assert b["status"] == "UNVERIFIED" and b["susceptible_max_mg_l"] is None and b["resistant_min_mg_l"] is None
        with pytest.raises(config.BreakpointUnverified):
            config.require_verified_breakpoints(cfg, tier)


def test_cohort_mapping(cfg):
    assert config.pmid_to_tier(cfg, "38052776") == "discovery"
    assert config.pmid_to_tier(cfg, 34485958) == "replication"
    assert config.pmid_to_tier(cfg, "38219757") == "exploratory"
    assert config.pmid_to_tier(cfg, "28720578") == "exploratory"
    assert config.pmid_to_tier(cfg, "unknown") == "exploratory"


def test_answer_key_pre_registered_hash_and_guards(tmp_path):
    key, sha = config.load_answer_key()
    assert key["status"] == "pre_registered" and key["registered_on"]
    assert sha == hashlib.sha256((config.CONFIG_DIR / "answer_key.yaml").read_bytes()).hexdigest()
    genes = {(d["gene"], d["position"]) for d in key["primary_determinants"]}
    assert {("gyrA", 83), ("gyrA", 87), ("parC", 80), ("parC", 84)} <= genes
    bad = tmp_path / "k.yaml"; bad.write_text(yaml.safe_dump(dict(key, status="draft")))
    with pytest.raises(config.ConfigError):
        config.load_answer_key(bad)


# ------------------------------------------------------------------ cohorts
def test_tiers_and_exclusion_reasons(syn):
    kept, ex = syn["kept"], syn["exclusions"]
    assert kept["cohort_tier"].value_counts().to_dict() == {"discovery": 499, "replication": 250, "exploratory": 150}
    r = ex["reason"].value_counts().to_dict()
    assert r == {"superset_duplicate": 60, "lower_tier_duplicate_genome": 1, "non_mic_method": 1, "exact_duplicate_record": 1, "conflicting_duplicate_mic": 2}


def test_discovery_and_replication_contain_only_their_pmids(syn):
    k = syn["kept"]
    assert set(k[k["cohort_tier"] == "discovery"]["study_pmid"]) == {"38052776"}
    assert set(k[k["cohort_tier"] == "replication"]["study_pmid"]) == {"34485958"}


def test_genome_in_two_tiers_goes_to_highest(syn):
    disc_first = syn["kept"]["genome_id"].iloc[0]
    g0 = [g for g in syn["features"].index if g.startswith("562.1")][0]
    row = syn["kept"][syn["kept"]["genome_id"] == g0].iloc[0]
    assert row["cohort_tier"] == "discovery" and row["study_pmid"] == "38052776"   # lower-tier 38219757 duplicate dropped


def test_superset_rule_keeps_38219757_drops_28720578(syn):
    k = syn["kept"]
    assert "28720578" not in set(k["study_pmid"])
    ex = syn["exclusions"]; assert (ex[ex["reason"] == "superset_duplicate"]["study_pmid"] == "28720578").all()


def test_conflicting_duplicates_excluded_entirely_not_averaged(syn):
    conflicted = syn["exclusions"][syn["exclusions"]["reason"] == "conflicting_duplicate_mic"]["genome_id"].unique()
    assert len(conflicted) == 1 and conflicted[0] not in set(syn["kept"]["genome_id"])


def test_one_row_per_genome_and_separation_guard(syn, cfg):
    k = syn["kept"]
    assert not k["genome_id"].duplicated().any()
    cohorts.assert_cohort_separation(k, cfg)
    bad = pd.concat([k, k.iloc[[0]].assign(cohort_tier="exploratory")], ignore_index=True)
    with pytest.raises(cohorts.CohortSeparationError):
        cohorts.assert_cohort_separation(bad, cfg)
    wrong = k.copy(); wrong.loc[wrong["study_pmid"] == "38052776", "cohort_tier"] = "replication"
    with pytest.raises(cohorts.CohortSeparationError):
        cohorts.assert_cohort_separation(wrong, cfg)


def test_assert_tier(syn):
    k = syn["kept"]
    cohorts.assert_tier(k[k["cohort_tier"] == "discovery"], "discovery", "t")
    with pytest.raises(cohorts.CohortSeparationError):
        cohorts.assert_tier(k, "discovery", "t")
    with pytest.raises(cohorts.CohortSeparationError):
        cohorts.assert_tier(k[k["cohort_tier"] == "exploratory"], "discovery", "t")


def test_overlap_matrix_counts(syn):
    m = cohorts.overlap_matrix(cohorts.normalise_ids(syn["records"]))
    assert m.loc["38219757", "28720578"] == 60 and m.loc["38052776", "34485958"] == 0
    assert m.loc["38052776", "38219757"] == 1


def test_invalid_mic_records_excluded(cfg):
    rec = pd.DataFrame({"genome_id": ["a", "b"], "study_pmid": ["38052776"] * 2, "laboratory_typing_method": ["MIC"] * 2, "measurement": ["<=0.25", "garbage"]})
    from variantbridge.amr import data_io
    k, e = cohorts.resolve_cohorts(data_io.fetch_phenotypes(cfg, records=rec), cfg)
    assert k["genome_id"].tolist() == ["a"] and e["reason"].tolist() == ["invalid_mic"]


# ------------------------------------------------------------------ QC
def _meta(**over):
    base = dict(genome_id="g", genome_name="Escherichia coli x", checkm_completeness=99.0, checkm_contamination=0.5, genome_length=5_000_000, contigs=100)
    base.update(over)
    return pd.DataFrame([base])


@pytest.mark.parametrize("over,gate", [
    ({"checkm_completeness": 90.0}, "low_completeness"), ({"checkm_contamination": 3.0}, "high_contamination"),
    ({"genome_length": 4_000_000}, "length_out_of_range"), ({"genome_length": 6_000_000}, "length_out_of_range"),
    ({"contigs": 600}, "too_many_contigs"), ({"genome_name": "Escherichia coli K-12 MG1655"}, "lab_or_reference_strain_name"),
    ({"checkm_completeness": np.nan}, "missing_checkm"), ({"genome_length": np.nan}, "missing_length"),
])
def test_each_qc_gate_excludes(cfg, over, gate):
    out = qc.genome_qc(_meta(**over), cfg)
    assert not out["qc_pass"].iloc[0] and gate in out["qc_fail_reasons"].iloc[0]


def test_good_genome_passes_and_boundaries_inclusive(cfg):
    assert qc.genome_qc(_meta(), cfg)["qc_pass"].iloc[0]
    assert qc.genome_qc(_meta(checkm_completeness=95.0, checkm_contamination=2.0), cfg)["qc_pass"].iloc[0]


def test_qc_thresholds_come_from_config_not_constants(cfg):
    c = copy.deepcopy(cfg); c["thresholds"]["genome_qc"]["checkm_completeness_min"]["value"] = 99.5
    assert not qc.genome_qc(_meta(), c)["qc_pass"].iloc[0]


def test_species_check_reported_not_run_unless_supplied(cfg, syn):
    r = qc.qc_report(qc.genome_qc(syn["meta"], cfg))
    assert "NOT RUN" in r["species_check"] and r["n_pass"] == r["n"] - 6
    m = syn["meta"].copy(); m["species_ok"] = True; m.loc[0, "species_ok"] = False
    r2 = qc.qc_report(qc.genome_qc(m, cfg)); assert r2["per_gate_failures"]["species_mismatch"] == 1
