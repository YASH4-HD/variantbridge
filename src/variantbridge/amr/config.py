"""Configuration loading with enforced threshold provenance (addendum Part 2, item 1)."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import yaml

CONFIG_DIR = Path(__file__).resolve().parents[3] / "config" / "amr"
PROVENANCE_CLASSES = ("source_derived", "pre_registered")
TIERS = ("discovery", "replication", "exploratory")


class ConfigError(ValueError):
    pass


class BreakpointUnverified(RuntimeError):
    """Raised by any breakpoint-dependent analysis while the breakpoint is not VERIFIED."""


def load_config(path: str | Path | None = None) -> dict:
    p = Path(path) if path else CONFIG_DIR / "config.yaml"
    with open(p) as fh:
        cfg = yaml.safe_load(fh)
    validate_thresholds(cfg)
    return cfg


def _walk_thresholds(block: dict, prefix: str = ""):
    for k, v in block.items():
        name = f"{prefix}{k}"
        if isinstance(v, dict) and ("provenance" in v or "value" in v):      # a leaf WITHOUT provenance must be caught, not skipped
            yield name, v
        elif isinstance(v, dict):
            yield from _walk_thresholds(v, name + ".")


def validate_thresholds(cfg: dict) -> None:
    """Every threshold must carry provenance; no threshold may be immutable-by-omission."""
    if "thresholds" not in cfg:
        raise ConfigError("config has no 'thresholds' block")
    n = 0
    for name, t in _walk_thresholds(cfg["thresholds"]):
        n += 1
        if "value" not in t:
            raise ConfigError(f"threshold {name}: missing 'value'")
        prov = t.get("provenance")
        if prov not in PROVENANCE_CLASSES:
            raise ConfigError(f"threshold {name}: provenance must be one of {PROVENANCE_CLASSES}, got {prov!r}")
        if prov == "source_derived" and not t.get("citation"):
            raise ConfigError(f"threshold {name}: source_derived requires a 'citation'")
        if prov == "pre_registered" and not (t.get("date") and t.get("rationale")):
            raise ConfigError(f"threshold {name}: pre_registered requires 'date' and 'rationale'")
    if n == 0:
        raise ConfigError("no thresholds found")


def threshold(cfg: dict, dotted: str) -> Any:
    node: Any = cfg["thresholds"]
    for part in dotted.split("."):
        node = node[part]
    return node["value"]


def threshold_provenance_table(cfg: dict) -> list[dict]:
    """Rows printed in the methods report: name, value, provenance class, date/citation, rationale."""
    rows = []
    for name, t in _walk_thresholds(cfg["thresholds"]):
        rows.append({"threshold": name, "value": t["value"], "provenance": t["provenance"],
                     "date_or_citation": t.get("date") or t.get("citation"), "rationale": t.get("rationale", "")})
    return rows


def pmid_to_tier(cfg: dict, pmid) -> str:
    m = cfg["cohorts"]["pmid_to_tier"]
    return m.get(str(pmid).strip(), cfg["cohorts"]["default_tier"])


def breakpoint_block(cfg: dict, cohort: str = "discovery") -> dict:
    key = f"ciprofloxacin_ecoli_{cohort}"
    try:
        return cfg["breakpoints"][key]
    except KeyError as exc:
        raise ConfigError(f"no breakpoint block {key!r}") from exc


def require_verified_breakpoints(cfg: dict, cohort: str = "discovery") -> tuple[float, float]:
    """Return (susceptible_max, resistant_min) ONLY if status == VERIFIED with provenance fields.

    Never infers a breakpoint. Any other state raises BreakpointUnverified (=> analysis NOT RUN)."""
    b = breakpoint_block(cfg, cohort)
    if b.get("status") != "VERIFIED":
        raise BreakpointUnverified(
            f"breakpoint for {cohort} cohort is {b.get('status')!r}; breakpoint-dependent analysis NOT RUN / UNVERIFIED")
    missing = [k for k in ("verified_by", "source_excerpt", "susceptible_max_mg_l", "resistant_min_mg_l") if b.get(k) in (None, "")]
    if missing:
        raise BreakpointUnverified(f"status says VERIFIED but fields are empty: {missing}; refusing")
    return float(b["susceptible_max_mg_l"]), float(b["resistant_min_mg_l"])


def file_sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def load_answer_key(path: str | Path | None = None) -> tuple[dict, str]:
    """Return (key, sha256). Refuses a key that is not marked pre_registered with a date."""
    p = Path(path) if path else CONFIG_DIR / "answer_key.yaml"
    with open(p) as fh:
        key = yaml.safe_load(fh)
    if key.get("status") != "pre_registered" or not key.get("registered_on"):
        raise ConfigError("answer_key.yaml must be status: pre_registered with registered_on date before association")
    if "primary_determinants" not in key:
        raise ConfigError("answer_key.yaml has no primary_determinants")
    return key, file_sha256(p)
