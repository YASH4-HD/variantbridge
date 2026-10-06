"""Evidence Cards (schema v1.1) with claim restrictions (addendum Part 2, item 4).

Four INDEPENDENT dimensions are carried and rendered separately and are never summed, averaged or
chained:  association_evidence | predictive_contribution | biological_plausibility | causal_status.
`validate_card` enforces the schema AND semantic rules; a card that violates them cannot be built."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

import jsonschema
import numpy as np
import pandas as pd

from .config import CONFIG_DIR

SCHEMA_PATH = CONFIG_DIR / "evidence_card_schema_v1_1.json"
DIMENSIONS = ("association_evidence", "predictive_contribution", "biological_plausibility", "causal_status")
SYNTHETIC_CLAIM = "SYNTHETIC DATA: software demonstration only; no biological claim is made."

_FORBIDDEN_UNLESS_CAUSAL = re.compile(
    r"\b(causes?|caused|causal(ly)?|causative|mechanis\w+|drives?|driven by|confers?|responsible for|proves?|proven|"
    r"validated|treatment|therapy|clinical(ly)?|diagnos\w+)\b", re.I)
_NOVEL = re.compile(r"\b(novel|new)\b.*\b(determinant|resistance gene|marker|resistance mechanism)\b", re.I)
_COMPOSITE_KEY = re.compile(r"(composite|combined|overall|total|priority|final|integrated).*(score|rank)|(score|rank).*(composite|combined|overall|total)", re.I)
_DIM_TERMS = {
    "association": re.compile(r"\bassociat\w+", re.I),
    "predictive": re.compile(r"\b(predict\w*|importance|AUC)\b", re.I),
    "plausibility": re.compile(r"\b(plausib\w+|efflux|pathway|target-site|annotat\w+)\b", re.I),
    "causal": re.compile(r"\b(causes?|causal\w*|mechanis\w+|therefore)\b", re.I),
}


class EvidenceCardError(ValueError):
    def __init__(self, problems: list[str]):
        super().__init__("; ".join(problems))
        self.problems = problems


def load_schema(path: Path | str = SCHEMA_PATH) -> dict:
    with open(path) as fh:
        return json.load(fh)


def _walk_keys(o, path=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield f"{path}.{k}" if path else k, k
            yield from _walk_keys(v, f"{path}.{k}" if path else k)
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from _walk_keys(v, f"{path}[{i}]")


def validate_card(card: dict, schema: Optional[dict] = None) -> None:
    """Raise EvidenceCardError listing EVERY problem (schema + claim restrictions)."""
    schema = schema or load_schema()
    problems: list[str] = []
    for e in sorted(jsonschema.Draft202012Validator(schema).iter_errors(card), key=lambda e: list(e.path)):
        problems.append(f"schema: {'/'.join(map(str, e.path)) or '<root>'}: {e.message}")
    if problems:                                   # semantic checks need a structurally valid card
        raise EvidenceCardError(problems)
    gen, causal = card["generated"], card["causal_status"]
    allowed = card["claim_currently_justified"]
    not_ok = card["claim_not_justified"]
    cs, ev = causal["status"], causal.get("experimental_evidence")
    # 1. causal_status consistency
    if cs == "experimentally_supported" and not (ev and str(ev).strip()):
        problems.append("causal_status experimentally_supported requires non-empty experimental_evidence")
    if cs == "not_established" and ev:
        problems.append("causal_status not_established must not carry experimental_evidence")
    # 2. no composite score anywhere
    for path, key in _walk_keys(card):
        if _COMPOSITE_KEY.search(key):
            problems.append(f"composite score/rank field {path!r} merges evidence dimensions (forbidden)")
    # 3. synthetic cards carry no biological claim
    if gen["data_scope"] == "synthetic" and allowed != [SYNTHETIC_CLAIM]:
        problems.append("synthetic card: claim_currently_justified must be exactly the synthetic disclaimer")
    # 4. forbidden vocabulary in allowed claims
    for c in allowed:
        if c == SYNTHETIC_CLAIM:
            continue
        if cs != "experimentally_supported" and _FORBIDDEN_UNLESS_CAUSAL.search(c):
            problems.append(f"claim not justified without experimental evidence: {c!r}")
        if _NOVEL.search(c):
            problems.append(f"'novel determinant/resistance gene' claims are never justified from association: {c!r}")
        if re.search(r"\bdeterminant\b", c, re.I) and not card["known_amr_database_evidence"]["is_answer_key_determinant"]:
            problems.append(f"'determinant' only allowed for pre-registered answer-key features: {c!r}")
        hit = [d for d, rx in _DIM_TERMS.items() if rx.search(c)]
        if len(hit) > 1:
            problems.append(f"one claim spans several evidence dimensions {hit}; dimensions must not be combined: {c!r}")
    # 5. adjustment honesty
    pa = card["population_adjusted_evidence"]
    if (pa["structure_control"] == "none" or not pa["survives_adjustment"]) and any(re.search(r"after lineage adjustment|lineage-adjusted|independent of lineage", c, re.I) for c in allowed):
        problems.append("claims 'after lineage adjustment' but the association does not survive / was not lineage-adjusted")
    # 6. predictive claims only from honest CV
    pc = card.get("predictive_contribution")
    pred_claim = any(_DIM_TERMS["predictive"].search(c) for c in allowed if c != SYNTHETIC_CLAIM)
    if pred_claim and (pc is None or pc.get("cv_strategy") == "random_demo_leakage_only"):
        problems.append("predictive claim without lineage-blocked / replication evidence (random-split is leakage-demo only)")
    # 7. association tier
    if gen["cohort_tier"] != "discovery" and any(_DIM_TERMS["association"].search(c) for c in allowed if c != SYNTHETIC_CLAIM):
        problems.append("primary association claims may only come from the discovery tier")
    # 8. mandatory negatives
    if cs != "experimentally_supported" and not any(re.search(r"caus", c, re.I) for c in not_ok):
        problems.append("claim_not_justified must explicitly forbid causal claims")
    if gen.get("breakpoint_status", "UNVERIFIED") != "VERIFIED" and not any(re.search(r"clinical|S/R|susceptib", c, re.I) for c in not_ok):
        problems.append("breakpoint not VERIFIED: claim_not_justified must forbid clinical/S-R category claims")
    if problems:
        raise EvidenceCardError(problems)


# ------------------------------------------------------------------------------------------------ builder
_ROLE_TO_PLAUSIBLE = {"primary_QRDR": "target_site_mutation", "secondary_efflux_regulator": "efflux_regulation", "secondary_PMQR": "horizontal_gene"}
_ROLE_TO_VALIDATION = {"target_site_mutation": "allele_replacement", "efflux_regulation": "efflux_phenotyping_assay", "horizontal_gene": "targeted_sequencing_cohort"}


def build_card(hit: pd.Series | dict, card_id: str, context: dict, annotation: Optional[dict] = None,
               prediction: Optional[dict] = None, replication: Optional[dict] = None,
               literature: Optional[list] = None, feature_meta: Optional[dict] = None,
               schema: Optional[dict] = None) -> dict:
    """Assemble and VALIDATE a card. `context`: analysis_date, pipeline_version, dataset_name, dataset_citation,
    random_seed, cohort_tier, censoring_handling, breakpoint_status, answer_key_sha256, data_scope,
    censoring_summary, alpha_fdr. `annotation`: {role, in_database, database, database_entry, is_answer_key_determinant, match}.
    `prediction`: None or {model, cv_strategy, validation_scope, outer_fold_auc, coefficient, rank, ...}."""
    h = dict(hit)
    ann = annotation or {}
    fm = feature_meta or {}
    alpha = context.get("alpha_fdr", 0.05)
    significant = bool(np.isfinite(h["q_value"]) and h["q_value"] < alpha)
    adjusted = bool(h.get("lineage_adjusted", False))
    survives = bool(significant and adjusted)
    role = ann.get("role")
    plaus = _ROLE_TO_PLAUSIBLE.get(role, "unknown")
    synthetic = context["data_scope"] == "synthetic"
    expected_pos = ann.get("is_answer_key_determinant")
    card = {
        "card_id": card_id,
        "generated": {k: context.get(k) for k in ("analysis_date", "pipeline_version", "dataset_name", "dataset_citation", "random_seed",
                                                    "cohort_tier", "censoring_handling", "breakpoint_status", "answer_key_sha256", "data_scope") if k in context},
        "candidate": {"feature_id": str(h["feature_id"]), "feature_type": h.get("feature_type", "gene_pa"), "gene": fm.get("gene"),
                      "product": fm.get("product"), "locus_tag": fm.get("locus_tag"), "reference_coordinates": fm.get("reference_coordinates"),
                      "allele_frequency_in_dataset": (float(h["n_present"]) / float(h["n"])) if h.get("n") else None},
        "phenotype": {"antibiotic": "ciprofloxacin", "phenotype_type": "MIC", "breakpoint_standard": None,
                      "n_resistant": None, "n_susceptible": None,
                      "class_balance_note": ("R/S counts withheld: breakpoint is " + str(context.get("breakpoint_status")) + "; censoring summary: " + json.dumps(context.get("censoring_summary", {})))},
        "association_evidence": {
            "model": h.get("model", "interval_censored_gaussian_Wald"), "beta": float(h["effect"]), "standard_error": float(h["se"]),
            "ci_low": float(h["ci_low"]), "ci_high": float(h["ci_high"]), "p_value": float(h["p_value"]),
            "correction_method": context.get("correction_method", "bh_fdr"), "q_value": float(h["q_value"]),
            "significant_after_correction": significant,
            "effect_direction_consistent_with_literature": (bool(h["effect"] > 0) if expected_pos else None),
            "effect_scale": "log2 MIC difference per feature present (interval-censored Gaussian model)",
            "censoring_summary": context.get("censoring_summary")},
        "population_adjusted_evidence": {
            "structure_control": h.get("structure_control", "lineage_fixed_effects" if adjusted else "none"),
            "survives_adjustment": survives,
            "unadjusted_p_value": (float(h["unadjusted_p"]) if pd.notna(h.get("unadjusted_p", np.nan)) else None),
            "lineage_association_note": ("unadjusted effect %.3f vs lineage-adjusted %.3f" % (h["unadjusted_effect"], h["effect"])) if pd.notna(h.get("unadjusted_effect", np.nan)) else None,
            "homoplasy_evidence": None},
        "known_amr_database_evidence": {"in_database": bool(ann.get("in_database", False)), "database": ann.get("database"),
                                         "database_entry": ann.get("database_entry"),
                                         "database_note": "AMR database annotation NOT RUN in this build" if not ann.get("in_database") and ann.get("database") is None else ann.get("database_note"),
                                         "is_answer_key_determinant": bool(ann.get("is_answer_key_determinant", False))},
        "predictive_contribution": None if prediction is None else {
            "model": prediction.get("model", "elastic_net_logistic_extreme_contrast"),
            "cv_strategy": prediction["cv_strategy"], "validation_scope": prediction.get("validation_scope", "lineage_blocked_cv_within_discovery"),
            "outer_fold_auc": prediction.get("outer_fold_auc"), "feature_coefficient_or_importance": prediction.get("coefficient"),
            "rank_in_model": prediction.get("rank"),
            "importance_interpretation_warning": "Feature importance reflects predictive utility under this CV design, not causal contribution."},
        "biological_plausibility": {"plausible_role": plaus,
                                    "rationale": ("Interpretation layer only (annotation.py, pre-registered answer-key label); it did not enter feature selection, the association model or any ranking."
                                                  if role else "No annotation-derived plausibility identified; absence of annotation is not evidence against association."),
                                    "pathway_membership": None, "plausibility_is_not_mechanism": True},
        "causal_status": {"status": "not_established", "experimental_evidence": None,
                          "note": "No allele-replacement or equivalent experimental evidence exists for this feature."},
        "replication_evidence": replication,
        "literature_support": list(literature or []),
        "uncertainty": _uncertainty(context, h, survives, replication),
        "evidence_still_missing": _missing(context, replication, plaus),
        "claim_currently_justified": [SYNTHETIC_CLAIM] if synthetic else _allowed_claims(h, significant, survives, replication, prediction, plaus, role),
        "claim_not_justified": _forbidden(context, ann),
        "suggested_experimental_validation": _suggest(plaus),
    }
    validate_card(card, schema)
    return card


def _allowed_claims(h, significant, survives, replication, prediction, plaus, role) -> list[str]:
    out = []
    if survives:
        s = "This feature is associated with ciprofloxacin MIC in the discovery cohort after lineage adjustment and is prioritised for follow-up."
        out.append(s)
        if replication and replication.get("replicated"):
            out.append("The direction of this association was observed again in the independent replication cohort (association evidence only).")
    elif significant:
        out.append("This feature is associated with ciprofloxacin MIC in the discovery cohort before lineage adjustment only; it may mark lineage rather than resistance.")
    else:
        out.append("No association with ciprofloxacin MIC is established for this feature in the discovery cohort.")
    if prediction and prediction.get("cv_strategy") in ("nested_lineage_blocked", "discovery_to_replication"):
        out.append("The feature was retained in the predictive baseline under lineage-blocked validation; this describes predictive utility only.")
    if role:
        out.append(f"Annotation: the feature matches the pre-registered key entry ({role}); this is an interpretation label.")
    return out


def _forbidden(context, ann) -> list[str]:
    f = ["This feature causes ciprofloxacin resistance.", "This feature acts through an established resistance mechanism.",
         "Predictive importance or biological plausibility implies causality."]
    if not ann.get("is_answer_key_determinant"):
        f.append("This feature is a novel resistance determinant.")
    if context.get("breakpoint_status") != "VERIFIED":
        f.append("Any clinical resistance category (S/I/R) or treatment statement: the breakpoint is not verified.")
    else:
        f.append("Any treatment or clinical-guidance statement.")
    return f


def _uncertainty(context, h, survives, replication) -> list[str]:
    u = ["Single discovery cohort; sampling frame enriched for 3GC-resistant / reduced-carbapenem-susceptibility isolates (addendum 1.1).",
         "Phenotype is censored at the assay limits; effect sizes are model-based under a Gaussian assumption on log2 MIC.",
         "Linked features cannot be separated; reported effects belong to a linkage cluster.",
         "Rare-variant filtering removes features below the frequency threshold (config).",
         f"Censoring summary: {json.dumps(context.get('censoring_summary', {}))}"]
    if context.get("breakpoint_status") != "VERIFIED":
        u.append("Breakpoint/AST standard for the discovery cohort is not verified from the source; no S/R analysis was performed.")
    if not survives:
        u.append("Association does not survive lineage adjustment at the pre-registered FDR.")
    if replication is None:
        u.append("Not tested in the replication cohort.")
    return u


def _missing(context, replication, plaus) -> list[str]:
    m = ["Literature review for this candidate (not performed by the pipeline).", "Allele replacement or equivalent experiment (causal evidence).",
         "MIC dose-response / isogenic comparison."]
    if replication is None or not replication.get("replicated"):
        m.insert(0, "Independent replication of direction and effect presence.")
    if context.get("breakpoint_status") != "VERIFIED":
        m.append("Verified breakpoint for any S/R-based analysis.")
    return m


def _suggest(plaus) -> dict:
    cls = _ROLE_TO_VALIDATION.get(plaus, "in_silico_replication_only")
    desc = {"allele_replacement": "Introduce the candidate allele into an isogenic E. coli background and measure the ciprofloxacin MIC shift (broth microdilution).",
            "efflux_phenotyping_assay": "Compare ethidium-bromide accumulation/efflux with and without an efflux-pump inhibitor in isogenic strains differing at the candidate.",
            "targeted_sequencing_cohort": "Screen an independent local isolate panel for the candidate and relate it to MIC.",
            "in_silico_replication_only": "No wet-lab experiment is proposed until the association replicates."}[cls]
    return {"validation_class": cls, "experiment_description": desc,
            "required_capability": "E. coli BSL-2 culture; CLSI broth microdilution; recombineering/cloning as applicable",
            "expected_result_if_hypothesis_true": "A reproducible MIC shift attributable to the candidate allele in the isogenic background.",
            "feasibility_note": "UNVERIFIED: laboratory capability, isolate source and sequencing budget must be confirmed with the supervisor."}


# --------------------------------------------------------------------------------- prioritisation (decomposable)
def rank_candidates(assoc: pd.DataFrame, annotation_cols: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Pure statistical ranking (BH q ascending, then |z| descending, then feature_id). Annotation columns are
    appended UNCHANGED beside it and never alter the order, and there is no combined score column."""
    ok = assoc[assoc["status"] == "ok"].copy()
    ok["_absz"] = ok["z"].abs()
    ok = ok.sort_values(["q_value", "_absz", "feature_id"], ascending=[True, False, True]).drop(columns="_absz").reset_index(drop=True)
    ok.insert(0, "association_rank", np.arange(1, len(ok) + 1))
    if annotation_cols is not None:
        add = annotation_cols.drop(columns=[c for c in annotation_cols.columns if c in ok.columns and c != "feature_id"], errors="ignore")
        ok = ok.merge(add, on="feature_id", how="left", sort=False)
    return ok
