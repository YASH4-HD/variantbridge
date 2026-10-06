"""Downloadable methods / status report with capability labels. Never reports a result that
was not executed: status comes only from artifacts/VALIDATION_STATUS.json."""
from __future__ import annotations

import datetime as _dt
import json

from . import __version__
from .capability import CAPABILITIES, VALIDATION_IDS, load_status


def methods_report(root="artifacts", manifests: dict | None = None) -> str:
    st = load_status(root)
    out = [f"# VariantBridge v{__version__} - methods and validation status",
           f"Generated {_dt.datetime.now(_dt.timezone.utc).isoformat(timespec='seconds')}", "",
           "## Capability tiers", ""]
    for k, (tier, vid) in CAPABILITIES.items():
        out.append(f"- **{k}**: {tier}" + (f" (validation {vid}: {st[vid]['status']})" if vid else ""))
    out += ["", "## Validation status (NOT RUN until executed on real data)", ""]
    for vid, desc in VALIDATION_IDS.items():
        rec = st[vid]
        out.append(f"- **{vid}** {desc}: **{rec['status']}**" + (f" (executed {rec.get('executed_at')})" if "executed_at" in rec else ""))
    if manifests:
        out += ["", "## Manifests", "", "```json", json.dumps(manifests, indent=2, default=str), "```"]
    out += ["", "No clinical, diagnostic or causal claims are made. Association is not causation; prediction is not mechanism."]
    return "\n".join(out) + "\n"
