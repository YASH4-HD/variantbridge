"""INTERPRETATION layer only. Nothing in this module feeds feature selection, the association
model, or a statistical ranking (addendum Part 2, items 2-3).

* label_answer_key: marks features that match the pre-registered key, AFTER the association run.
* parse_amrfinder_point_mutations: reads AMRFinderPlus output (tool NOT run in this build; the
  column names are those documented for AMRFinderPlus and are UNVERIFIED here).
* qrdr_residue_check: sanity check that a supplied reference protein carries the expected residue at
  a QRDR position (guards against off-by-one coordinate errors). No reference sequence is embedded
  here - sequences must be supplied from the reference accession GCF_000005845.2."""
from __future__ import annotations

import re
from typing import Optional

import pandas as pd

_AA = re.compile(r"^(?P<gene>[A-Za-z0-9()'\-]+?)[_:\s]+(?P<ref>[A-Z])(?P<pos>\d+)(?P<alt>[A-Z*])$")


def parse_aa_feature(feature_id: str) -> Optional[dict]:
    """'gyrA_S83L' / 'gyrA:S83L' -> {gene, ref, pos, alt}; None if the id is not an AA-substitution id."""
    m = _AA.match(str(feature_id).strip())
    return None if not m else {"gene": m["gene"], "ref": m["ref"], "pos": int(m["pos"]), "alt": m["alt"]}


def key_entries(key: dict) -> list[dict]:
    rows = []
    for d in key["primary_determinants"]:
        rows.append({"role": "primary_QRDR", "gene": d["gene"], "ref": d["reference_residue"], "pos": int(d["position"]), "alts": list(d["alt_residues"])})
    sec = key.get("secondary_determinants", {})
    for g in sec.get("pmqr_genes", []):
        rows.append({"role": "secondary_PMQR", "gene": g, "ref": None, "pos": None, "alts": []})
    for g in sec.get("efflux_regulator_genes", []):
        rows.append({"role": "secondary_efflux_regulator", "gene": g, "ref": None, "pos": None, "alts": []})
    return rows


def label_answer_key(assoc: pd.DataFrame, feature_meta: Optional[pd.DataFrame], key: dict) -> pd.DataFrame:
    """Add `in_answer_key`, `answer_key_role`, `answer_key_match` columns. feature_meta (optional) has
    feature_id, gene, aa_ref, aa_pos, aa_alt. Presence/absence of PMQR genes matches by gene name; efflux
    regulators match any feature whose `gene` is the regulator. Pure labelling; row order, p and q unchanged."""
    out = assoc.copy()
    meta = feature_meta.set_index("feature_id") if feature_meta is not None else pd.DataFrame()
    entries = key_entries(key)
    role, match = [], []
    for fid in out["feature_id"]:
        r = m = None
        info = parse_aa_feature(fid)
        gene = ref = pos = alt = None
        if info:
            gene, ref, pos, alt = info["gene"], info["ref"], info["pos"], info["alt"]
        if fid in meta.index:
            row = meta.loc[fid]
            gene = row.get("gene", gene) if pd.notna(row.get("gene", None)) else gene
            if pd.notna(row.get("aa_pos", None)):
                ref, pos, alt = row.get("aa_ref"), int(row["aa_pos"]), row.get("aa_alt")
        for e in entries:
            if e["role"] == "primary_QRDR":
                if gene == e["gene"] and pos == e["pos"] and ref == e["ref"] and alt in e["alts"]:
                    r, m = e["role"], f"{e['gene']} {e['ref']}{e['pos']}{alt}"; break
            else:
                by_gene = gene is not None and str(gene).lower() == e["gene"].lower()
                by_id = re.fullmatch(re.escape(e["gene"]) + r"\d*([_:].*)?", str(fid), flags=re.I) is not None
                if by_gene or by_id:
                    r, m = e["role"], e["gene"]; break
        role.append(r); match.append(m)
    out["answer_key_role"] = role
    out["answer_key_match"] = match
    out["in_answer_key"] = out["answer_key_role"].notna()
    return out


def qrdr_residue_check(protein_seq: str, position: int, expected_residue: str) -> bool:
    """1-based position check on a supplied reference protein sequence."""
    return 0 < position <= len(protein_seq) and protein_seq[position - 1].upper() == expected_residue.upper()


def parse_amrfinder_point_mutations(df: pd.DataFrame) -> pd.DataFrame:
    """Expect AMRFinderPlus TSV columns including 'Gene symbol', 'Element type', 'Element subtype',
    and a name/sample column. Returns tidy rows (sample, gene, ref, pos, alt) for POINT mutations.
    Raises on missing columns (UNVERIFIED format - fail loudly)."""
    need = {"Gene symbol"}
    missing = need - set(df.columns)
    if missing:
        raise KeyError(f"AMRFinderPlus table missing columns {sorted(missing)}")
    sample_col = next((c for c in ("Name", "Contig id", "sample") if c in df.columns), None)
    if sample_col is None:
        raise KeyError("AMRFinderPlus table needs a sample/Name column")
    rows = []
    for _, r in df.iterrows():
        sym = str(r["Gene symbol"])
        m = re.match(r"^(?P<gene>[A-Za-z0-9]+)_(?P<ref>[A-Z])(?P<pos>\d+)(?P<alt>[A-Z*])$", sym)
        if m:
            rows.append({"sample": r[sample_col], "gene": m["gene"], "ref": m["ref"], "pos": int(m["pos"]), "alt": m["alt"]})
    return pd.DataFrame(rows, columns=["sample", "gene", "ref", "pos", "alt"])
