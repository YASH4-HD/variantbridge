"""Loads compact validated EvoResist artifacts for display. No computation of results here; the Streamlit
app never runs association/prediction. Missing/invalid artifacts are reported (NOT RUN), never substituted."""
from __future__ import annotations

import json
from pathlib import Path

from ..artifacts import load_table_artifact
from . import evidence as E
from .validation import load_status

NAMES = {"phenotypes": "amr_phenotypes.csv.gz", "exclusions": "amr_exclusions.csv.gz", "association": "amr_association_discovery.csv.gz",
         "replication": "amr_replication.csv.gz", "importance": "amr_prediction_importance.csv.gz", "folds": "amr_prediction_folds.csv.gz"}


def load_all(root: str | Path) -> dict:
    root = Path(root)
    out = {"root": root, "status": load_status(root), "messages": {}}
    for k, n in NAMES.items():
        df, m, msg = load_table_artifact(n, root)
        out[k], out[k + "_manifest"] = df, m
        if msg:
            out["messages"][k] = msg
    out["cards"], out["cards_rejected"], out["cards_scope"] = [], [], None
    cp = root / "amr_evidence_cards.json"
    if cp.exists():
        doc = json.loads(cp.read_text())
        out["cards_scope"] = doc.get("data_scope")
        schema = E.load_schema()
        for c in doc.get("cards", []):
            try:                                   # re-validate at load: a card that violates the rules is never displayed as valid
                E.validate_card(c, schema)
                out["cards"].append(c)
            except E.EvidenceCardError as e:
                out["cards_rejected"].append({"card_id": c.get("card_id"), "problems": e.problems})
        out["cards_rejected"] += doc.get("rejected", [])
    else:
        out["messages"]["cards"] = "amr_evidence_cards.json not found (script amr_06 NOT RUN)"
    ak = root / "amr_answer_key_validation.json"
    out["answer_key_validation"] = json.loads(ak.read_text()) if ak.exists() else None
    return out


def scopes(data: dict) -> set:
    s = set()
    for k in ("phenotypes", "association", "replication", "importance"):
        m = data.get(k + "_manifest")
        if m and m.get("data_scope"):
            s.add(m["data_scope"])
    if data.get("cards_scope"):
        s.add(data["cards_scope"])
    return s
