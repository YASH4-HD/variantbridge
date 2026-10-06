"""BV-BRC acquisition and local caching. OFFLINE use only (never called by Streamlit).

Network access is injected (`http_get`) so the logic is testable without a network and so a
blocked network yields an explicit DataNotAvailable (=> NOT RUN) instead of silent fallbacks.
BV-BRC field names and paging syntax below follow the audit text where stated and are otherwise
UNVERIFIED against the live API in this build (the build environment could not reach BV-BRC);
`FIELD_ALIASES` makes the first real run fail loudly on a mismatch rather than guess.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Callable, Iterable, Optional

import pandas as pd

from .cohorts import assign_record_tier
from .config import load_config
from .phenotypes import parse_mic_table

API_AMR = "https://www.bv-brc.org/api/genome_amr/"
API_GENOME = "https://www.bv-brc.org/api/genome/"
API_SEQ = "https://www.bv-brc.org/api/genome_sequence/"
PAGE = 25000

# canonical name -> acceptable raw names (UNVERIFIED for the live API; resolved case-insensitively)
FIELD_ALIASES = {
    "genome_id": ["genome_id"],
    "genome_name": ["genome_name"],
    "antibiotic": ["antibiotic"],
    "measurement": ["measurement"],
    "measurement_sign": ["measurement_sign"],
    "measurement_value": ["measurement_value"],
    "measurement_unit": ["measurement_unit"],
    "resistant_phenotype": ["resistant_phenotype"],
    "laboratory_typing_method": ["laboratory_typing_method"],
    "testing_standard": ["testing_standard"],
    "testing_standard_year": ["testing_standard_year"],
    "evidence": ["evidence"],
    "pmid": ["pmid", "source_pmid"],
}


class DataNotAvailable(RuntimeError):
    """Required raw data is absent and network access was not allowed/possible => NOT RUN."""


def build_rql(cfg: dict, offset: int = 0, limit: Optional[int] = None) -> str:
    limit = PAGE if limit is None else limit
    return f'{cfg["dataset"]["query"]}&limit({limit},{offset})'


def fetch_genome_amr(cfg: dict, raw_dir: str | Path, http_get: Optional[Callable[[str], bytes]] = None,
                     allow_network: bool = False, max_pages: int = 50) -> list[Path]:
    """Resumable paged download of genome_amr records into raw_dir/genome_amr_page_NNNN.json.

    Existing pages are never re-downloaded. Stops at the first short page. Logs SHA-256 per page."""
    raw = Path(raw_dir)
    raw.mkdir(parents=True, exist_ok=True)
    pages: list[Path] = []
    for n in range(max_pages):
        p = raw / f"genome_amr_page_{n:04d}.json"
        if not p.exists():
            if not allow_network or http_get is None:
                if n == 0:
                    raise DataNotAvailable("no cached BV-BRC genome_amr pages and network not allowed (NOT RUN)")
                break
            body = http_get(API_AMR + "?" + build_rql(cfg, n * PAGE))
            tmp = p.with_suffix(".tmp")
            tmp.write_bytes(body)
            tmp.replace(p)
            _log(raw, p, body)
        pages.append(p)
        if len(json.loads(p.read_text())) < PAGE:
            break
    return pages


def _log(raw: Path, p: Path, body: bytes) -> None:
    lp = raw / "DOWNLOAD_LOG.json"
    log = json.loads(lp.read_text()) if lp.exists() else {}
    log[p.name] = {"sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body), "downloaded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    lp.write_text(json.dumps(log, indent=2))


def _resolve(df: pd.DataFrame) -> pd.DataFrame:
    low = {c.lower(): c for c in df.columns}
    ren = {}
    for canon, alts in FIELD_ALIASES.items():
        for a in alts:
            if a.lower() in low:
                ren[low[a.lower()]] = canon
                break
    out = df.rename(columns=ren)
    missing = [c for c in ("genome_id", "measurement", "pmid") if c not in out.columns]
    if missing:
        raise KeyError(f"BV-BRC records lack expected fields {missing}; field names are UNVERIFIED - check FIELD_ALIASES")
    return out


def load_cached_records(pages: Iterable[Path]) -> pd.DataFrame:
    rows = []
    for p in pages:
        rows.extend(json.loads(Path(p).read_text()))
    if not rows:
        raise DataNotAvailable("cached BV-BRC pages contain no records")
    return _resolve(pd.DataFrame(rows))


def fetch_phenotypes(config: dict | None = None, raw_dir: str | Path = "data/raw/bvbrc", http_get=None,
                     allow_network: bool = False, records: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Addendum interface: returns one row per RECORD with explicit censoring columns and `cohort_tier`.

    (`cohorts.resolve_cohorts` then applies duplicate exclusion.) Pass `records` to bypass I/O."""
    cfg = config or load_config()
    if records is None:
        pages = fetch_genome_amr(cfg, raw_dir, http_get, allow_network)
        records = load_cached_records(pages)
    records = records.copy()
    if "pmid" in records and "study_pmid" not in records:
        records = records.rename(columns={"pmid": "study_pmid"})
    parsed = parse_mic_table(records, grid_tol_log2=cfg["phenotype"]["dilution_grid_tolerance_log2"])
    out = assign_record_tier(parsed, cfg)
    out["batch"] = out["study_pmid"]
    return out.rename(columns={"record_tier": "cohort_tier"})


def contigs_json_to_fasta(contigs: list[dict], genome_id: str, width: int = 80) -> str:
    """Reconstruct FASTA from genome_sequence JSON (fields `sequence_id`/`accession`, `sequence`)."""
    lines = []
    for i, c in enumerate(contigs):
        name = c.get("sequence_id") or c.get("accession") or f"contig{i + 1}"
        seq = str(c["sequence"]).upper()
        bad = set(seq) - set("ACGTNRYKMSWBDHV")
        if bad:
            raise ValueError(f"{genome_id}/{name}: invalid sequence characters {sorted(bad)}")
        lines.append(f">{genome_id}|{name}")
        lines.extend(seq[j:j + width] for j in range(0, len(seq), width))
    return "\n".join(lines) + "\n"


META_FIELDS = ["genome_id", "genome_name", "mlst", "checkm_completeness", "checkm_contamination", "genome_length",
               "contigs", "collection_date", "isolation_country", "sra_accession", "genome_status"]   # names UNVERIFIED vs live API


def fetch_genome_metadata(genome_ids: list[str], raw_dir: str | Path, http_get: Optional[Callable[[str], bytes]] = None,
                          allow_network: bool = False, batch: int = 200) -> pd.DataFrame:
    """Resumable batched download of genome metadata (MLST, CheckM, length, contigs, ...). RQL syntax UNVERIFIED."""
    raw = Path(raw_dir)
    raw.mkdir(parents=True, exist_ok=True)
    ids = sorted(set(map(str, genome_ids)))
    rows = []
    for i in range(0, len(ids), batch):
        p = raw / f"genome_meta_{i // batch:05d}.json"
        if not p.exists():
            if not allow_network or http_get is None:
                raise DataNotAvailable(f"genome metadata batch {p.name} not cached and network not allowed (NOT RUN)")
            q = f'in(genome_id,({",".join(ids[i:i + batch])}))&select({",".join(META_FIELDS)})&limit({batch})'
            body = http_get(API_GENOME + "?" + q)
            tmp = p.with_suffix(".tmp"); tmp.write_bytes(body); tmp.replace(p)
            _log(raw, p, body)
        rows.extend(json.loads(p.read_text()))
    return pd.DataFrame(rows)
