import copy
import json

import numpy as np
import pandas as pd
import pytest

from variantbridge.amr import annotation as AN, evidence as E

CTX = {"analysis_date": "2026-10-05", "pipeline_version": "EvoResist-AI 0.2.0", "dataset_name": "synthetic", "dataset_citation": "n/a", "random_seed": 0,
       "cohort_tier": "discovery", "censoring_handling": "interval_censored_gaussian", "breakpoint_status": "UNVERIFIED", "answer_key_sha256": "0" * 64,
       "data_scope": "real", "censoring_summary": {"n": 100, "exact": 50, "left": 30, "right": 20}, "alpha_fdr": 0.05}


def _hit(**o):
    h = {"feature_id": "gyrA_S83L", "feature_type": "amino_acid_substitution", "n": 500, "n_present": 200, "effect": 3.0, "se": 0.2, "ci_low": 2.6, "ci_high": 3.4,
         "p_value": 1e-50, "q_value": 1e-48, "unadjusted_p": 1e-30, "unadjusted_effect": 2.5, "lineage_adjusted": True, "structure_control": "lineage_fixed_effects", "model": "interval_censored_gaussian_Wald"}
    h.update(o); return h


ANN = {"role": "primary_QRDR", "in_database": False, "database": None, "is_answer_key_determinant": True}
PRED = {"cv_strategy": "nested_lineage_blocked", "validation_scope": "lineage_blocked_cv_within_discovery", "outer_fold_auc": 0.9, "coefficient": 1.2, "rank": 1}
REP = {"cohort_pmid": "34485958", "effect": 2.0, "standard_error": 0.3, "p_value": 1e-9, "same_direction": True, "replicated": True, "status": "ok"}


@pytest.fixture
def card():
    return E.build_card(_hit(), "EC-0001", CTX, ANN, PRED, REP, None, {"gene": "gyrA"})


def test_valid_real_card_builds_and_validates(card):
    E.validate_card(card)
    assert card["causal_status"]["status"] == "not_established"
    assert any("after lineage adjustment" in c for c in card["claim_currently_justified"])
    assert card["literature_support"] == []                           # never fabricated
    assert card["phenotype"]["n_resistant"] is None                   # no S/R counts while breakpoint unverified


def test_four_dimensions_are_separate_top_level_fields_and_no_composite(card):
    for d in E.DIMENSIONS:
        assert d in card
    assert not [k for k in card if "score" in k.lower()]
    sections = {"association_evidence", "predictive_contribution", "biological_plausibility", "causal_status"}
    assert sections <= set(card)
    # plausibility text does not mention association statistics / prediction and vice versa
    assert "p_value" not in json.dumps(card["biological_plausibility"]) and "auc" not in json.dumps(card["biological_plausibility"]).lower()
    assert "plausible_role" not in card["association_evidence"] and "beta" not in card["biological_plausibility"]


def test_schema_v1_1_is_additive_over_supplied_v1():
    v1 = json.load(open(E.CONFIG_DIR / "evidence_card_schema_v1.json")); v11 = E.load_schema()
    assert set(v1["properties"]) <= set(v11["properties"]) and set(v1["required"]) <= set(v11["required"])
    assert {"causal_status", "predictive_contribution"} <= set(v11["required"])
    for k, sub in v1["properties"].items():
        if sub.get("type") == "object" and "required" in sub:
            assert set(sub["required"]) <= set(v11["properties"][k].get("required", [])), k


@pytest.mark.parametrize("path", ["uncertainty", "evidence_still_missing", "claim_currently_justified", "claim_not_justified", "causal_status",
                                  "association_evidence", "population_adjusted_evidence", "known_amr_database_evidence", "biological_plausibility",
                                  "suggested_experimental_validation", "literature_support", "phenotype", "candidate", "predictive_contribution"])
def test_schema_rejects_missing_required_fields(card, path):
    bad = copy.deepcopy(card); del bad[path]
    with pytest.raises(E.EvidenceCardError, match="schema"):
        E.validate_card(bad)


def _viol(card, mutate, needle):
    bad = copy.deepcopy(card); mutate(bad)
    with pytest.raises(E.EvidenceCardError) as ei:
        E.validate_card(bad)
    assert any(needle in p for p in ei.value.problems), ei.value.problems


def test_causal_and_mechanistic_claims_rejected(card):
    for txt in ("This mutation causes resistance.", "This confers high-level resistance.", "The causal variant is gyrA S83L.",
                "This is the mechanism of resistance.", "Validated resistance marker.", "Useful for clinical treatment decisions."):
        _viol(card, lambda c, t=txt: c["claim_currently_justified"].append(t), "claim not justified")


def test_novel_determinant_claim_rejected(card):
    _viol(card, lambda c: c["claim_currently_justified"].append("This is a novel resistance determinant."), "novel")
    h = E.build_card(_hit(feature_id="gene_042"), "EC-0002", CTX, None, None, None, None, None)      # not in answer key
    _viol(h, lambda c: c["claim_currently_justified"].append("This gene is a determinant of resistance."), "determinant")


def test_combined_dimension_statement_is_violation(card):
    _viol(card, lambda c: c["claim_currently_justified"].append("It is predictive and biologically plausible, therefore causal."), "dimensions")
    _viol(card, lambda c: c["claim_currently_justified"].append("The associated feature is predictive of MIC."), "dimensions")
    _viol(card, lambda c: c["claim_currently_justified"].append("Associated with MIC and annotated as an efflux regulator."), "dimensions")


def test_composite_score_fields_rejected(card):
    for key in ("composite_score", "priority_score", "overall_rank", "combined_evidence_score"):
        _viol(card, lambda c, k=key: c.__setitem__(k, 0.9), "composite")
    _viol(card, lambda c: c["association_evidence"].__setitem__("total_score", 1.0), "composite")


def test_adjustment_honesty(card):
    unadj = E.build_card(_hit(lineage_adjusted=False, structure_control="none", q_value=1e-5), "EC-0003", CTX, None, None, None, None, None)
    assert not any("after lineage adjustment" in c for c in unadj["claim_currently_justified"])
    assert any("before lineage adjustment only" in c for c in unadj["claim_currently_justified"])
    _viol(unadj, lambda c: c["claim_currently_justified"].append("Associated with MIC after lineage adjustment."), "lineage")
    ns = E.build_card(_hit(q_value=0.4), "EC-0004", CTX, None, None, None, None, None)
    assert ns["claim_currently_justified"][0].startswith("No association")


def test_predictive_claims_need_honest_cv(card):
    bad = E.build_card(_hit(), "EC-0005", CTX, ANN, dict(PRED, cv_strategy="random_demo_leakage_only", validation_scope="random_split_leakage_demo"), None, None, None)
    assert not any("predictive" in c for c in bad["claim_currently_justified"])
    _viol(bad, lambda c: c["claim_currently_justified"].append("Retained in the predictive baseline."), "predictive claim")
    nopred = E.build_card(_hit(), "EC-0006", CTX, ANN, None, None, None, None)
    assert nopred["predictive_contribution"] is None and not any("predictive" in c for c in nopred["claim_currently_justified"])


def test_causal_status_consistency(card):
    _viol(card, lambda c: c["causal_status"].update(status="experimentally_supported"), "requires non-empty")
    _viol(card, lambda c: c["causal_status"].update(experimental_evidence="made up"), "must not carry")
    ok = copy.deepcopy(card); ok["causal_status"] = {"status": "experimentally_supported", "experimental_evidence": "allele replacement: MIC x16 (test fixture)"}
    ok["claim_not_justified"] = [c for c in ok["claim_not_justified"] if "caus" not in c.lower()] + ["Any clinical resistance category (S/I/R) or treatment statement."]
    ok["claim_currently_justified"].append("Allele replacement shifts MIC (experimental, test fixture).")
    E.validate_card(ok)       # causal wording allowed only with experimental evidence


def test_mandatory_negative_claims(card):
    _viol(card, lambda c: c.__setitem__("claim_not_justified", ["Something else."]), "forbid causal")
    _viol(card, lambda c: c.__setitem__("claim_not_justified", [x for x in c["claim_not_justified"] if "clinical" not in x.lower()]), "clinical")


def test_synthetic_cards_carry_no_biological_claim():
    c = E.build_card(_hit(), "EC-0007", dict(CTX, data_scope="synthetic"), ANN, PRED, REP, None, None)
    assert c["claim_currently_justified"] == [E.SYNTHETIC_CLAIM]
    _viol(c, lambda x: x["claim_currently_justified"].append("Associated with MIC."), "synthetic")


def test_non_discovery_tier_cannot_make_association_claims(card):
    _viol(card, lambda c: c["generated"].__setitem__("cohort_tier", "exploratory"), "discovery tier")


def test_plausibility_never_changes_statistical_ranking():
    a = pd.DataFrame({"feature_id": ["a", "b", "c", "d"], "status": "ok", "q_value": [0.01, 0.2, 0.001, 0.2], "z": [3, 1, 5, -2.0], "effect": [1, .1, 2, -.3]})
    ann1 = pd.DataFrame({"feature_id": list("abcd"), "in_answer_key": [False, True, False, True], "answer_key_role": [None, "primary_QRDR", None, "secondary_efflux_regulator"]})
    ann2 = pd.DataFrame({"feature_id": list("abcd"), "in_answer_key": [True, False, True, False], "answer_key_role": ["secondary_efflux_regulator", None, "secondary_PMQR", None]})
    r0, r1, r2 = E.rank_candidates(a), E.rank_candidates(a, ann1), E.rank_candidates(a, ann2)
    assert r0["feature_id"].tolist() == r1["feature_id"].tolist() == r2["feature_id"].tolist() == ["c", "a", "d", "b"]
    assert not [c for c in r1.columns if "score" in c.lower()]                       # decomposable: no merged score
    assert {"association_rank", "q_value", "z"} <= set(r1.columns) and {"in_answer_key", "answer_key_role"} <= set(r1.columns)


def test_answer_key_labels_are_after_the_fact_and_do_not_alter_stats(assoc, syn, answer_key):
    lab = AN.label_answer_key(assoc, syn["feature_meta"], answer_key[0])
    for c in ("effect", "se", "p_value", "q_value"):
        assert np.array_equal(lab[c].to_numpy(), assoc[c].to_numpy())
    r = lab.set_index("feature_id")
    assert r.loc["gyrA_S83L", "answer_key_role"] == "primary_QRDR" and r.loc["parC_S80I", "answer_key_role"] == "primary_QRDR"
    assert r.loc["acrR_Q15*", "answer_key_role"] == "secondary_efflux_regulator" and r.loc["qnrS1", "answer_key_role"] == "secondary_PMQR"
    assert not r.loc["gene_000", "in_answer_key"] and not r.loc["lineage_marker_X", "in_answer_key"]
    assert AN.parse_aa_feature("gyrA_S83L") == {"gene": "gyrA", "ref": "S", "pos": 83, "alt": "L"} and AN.parse_aa_feature("gene_001") is None
    assert AN.qrdr_residue_check("MSDLAREITPVNIEEE", 6, "R") and not AN.qrdr_residue_check("MSD", 10, "R")


def test_amrfinder_parser_fails_loudly_on_unexpected_format():
    with pytest.raises(KeyError):
        AN.parse_amrfinder_point_mutations(pd.DataFrame({"x": [1]}))
    ok = AN.parse_amrfinder_point_mutations(pd.DataFrame({"Name": ["s1", "s1"], "Gene symbol": ["gyrA_S83L", "qnrS1"]}))
    assert ok.to_dict("records") == [{"sample": "s1", "gene": "gyrA", "ref": "S", "pos": 83, "alt": "L"}]
