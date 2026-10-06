"""Evidence Cards (handoff section 8): teach the difference between association, prediction,
mechanism and causality. A card never claims causation or clinical utility."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
import pandas as pd

NOT_ASSESSED = "not assessed"

ALWAYS_NOT_JUSTIFIED = [
    "causality of this variant/feature for the phenotype",
    "a confirmed molecular mechanism",
    "clinical, diagnostic or therapeutic utility",
]


@dataclass
class EvidenceCard:
    candidate: str
    phenotype: str
    effect: float
    effect_scale: str
    ci_low: float
    ci_high: float
    p: float
    multiple_testing: str
    adjustment: str
    data_tier: str
    replication_status: str = "none (single dataset)"
    known_database_evidence: str = NOT_ASSESSED
    predictive_importance: str = NOT_ASSESSED
    biological_plausibility: str = NOT_ASSESSED
    uncertainty: str = ""
    evidence_missing: list = field(default_factory=list)
    suggested_validation_class: str = "independent replication, then functional/experimental follow-up"
    claim_justified: str = ""
    claims_not_justified: list = field(default_factory=lambda: list(ALWAYS_NOT_JUSTIFIED))

    def as_dict(self) -> dict:
        return asdict(self)


def _fmt(x):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.4g}"


def make_card(row: pd.Series, phenotype: str, adjustment: str, data_tier: str, p_adjusted_col: str = "bonferroni",
              extra_missing: Optional[list] = None, predictive_importance: Optional[float] = None,
              known_annotation: Optional[str] = None) -> EvidenceCard:
    logistic = str(row.get("model", "")).startswith("logistic")
    scale = "log-odds (OR = exp(effect))" if logistic else "beta (phenotype units per allele/feature unit)"
    padj = float(row.get(p_adjusted_col, np.nan))
    survives = np.isfinite(padj) and padj < 0.05
    mt = f"{p_adjusted_col} = {_fmt(padj)} ({'survives' if survives else 'does not survive'} 0.05 correction across scanned features)"
    missing = ["independent replication", "domain-appropriate association model (relatedness/lineage)",
               "functional interpretation", "experimental validation"] + (extra_missing or [])
    just = (f"{row['variant']} is associated with {phenotype} in this dataset"
            + (" and survives multiple-testing correction" if survives else " (nominal; not corrected-significant)")
            + "; it is prioritised for follow-up.")
    return EvidenceCard(
        candidate=str(row["variant"]), phenotype=phenotype, effect=float(row["effect"]), effect_scale=scale,
        ci_low=float(row["ci_low"]), ci_high=float(row["ci_high"]), p=float(row["p"]), multiple_testing=mt,
        adjustment=adjustment, data_tier=data_tier,
        known_database_evidence=known_annotation or NOT_ASSESSED,
        predictive_importance=_fmt(predictive_importance) if predictive_importance is not None else NOT_ASSESSED,
        uncertainty=f"95% CI [{_fmt(row['ci_low'])}, {_fmt(row['ci_high'])}] from the same fitted model as the effect and p-value",
        evidence_missing=missing, claim_justified=just)


def render_text(c: EvidenceCard) -> str:
    lines = [
        f"EVIDENCE CARD - {c.candidate}", f"Data tier: {c.data_tier}", f"Phenotype: {c.phenotype}",
        f"Effect: {_fmt(c.effect)} [{c.effect_scale}]  95% CI [{_fmt(c.ci_low)}, {_fmt(c.ci_high)}]",
        f"p-value: {_fmt(c.p)}", f"Multiple testing: {c.multiple_testing}", f"Adjustment: {c.adjustment}",
        f"Replication: {c.replication_status}", f"Known database evidence: {c.known_database_evidence}",
        f"Predictive importance: {c.predictive_importance}", f"Biological plausibility: {c.biological_plausibility}",
        f"Uncertainty: {c.uncertainty}", "Evidence missing: " + "; ".join(c.evidence_missing),
        f"Suggested validation class: {c.suggested_validation_class}",
        f"CLAIM CURRENTLY JUSTIFIED: {c.claim_justified}",
        "CLAIMS NOT JUSTIFIED: " + "; ".join(c.claims_not_justified),
    ]
    return "\n".join(lines) + "\n"
