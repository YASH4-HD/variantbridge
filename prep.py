"""Helpers for the OFFLINE preparation scripts (never imported by the Streamlit app).

File-format parsing is defensive: column names are resolved through explicit alias lists, and
any unresolved/ambiguous header makes the script fail loudly with the observed header instead
of guessing. Format assumptions that the audit did not verify are marked UNVERIFIED in
docs/DATASETS.md and must be checked against the real files on first use.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


class PrepError(RuntimeError):
    pass


def find_plink2(explicit: Optional[str] = None) -> str:
    cand = explicit or os.environ.get("PLINK2") or shutil.which("plink2")
    if not cand or not (Path(cand).exists() or shutil.which(cand)):
        raise PrepError("plink2 not found. Install PLINK 2 (https://www.cog-genomics.org/plink/2.0/), "
                        "put it on PATH or set the PLINK2 environment variable / --plink2 argument.")
    return cand


def plink2_version(plink2: str) -> str:
    out = subprocess.run([plink2, "--version"], capture_output=True, text=True)
    return out.stdout.strip() or out.stderr.strip()


def run_plink2(plink2: str, args: list[str], log=None) -> subprocess.CompletedProcess:
    cmd = [plink2, *args]
    cp = subprocess.run(cmd, capture_output=True, text=True)
    if log is not None:
        log.append({"cmd": " ".join(cmd), "returncode": cp.returncode})
    if cp.returncode != 0:
        raise PrepError(f"plink2 failed ({cp.returncode}): {' '.join(cmd)}\n{cp.stdout[-1500:]}\n{cp.stderr[-1500:]}")
    return cp


def resolve_columns(df_columns: Iterable[str], aliases: dict[str, list[str]], context: str) -> dict[str, str]:
    """Map canonical names to observed columns via case-insensitive alias lists; fail loudly."""
    cols = {c.lower().strip(): c for c in df_columns}
    out, problems = {}, []
    for canon, names in aliases.items():
        hit = list(dict.fromkeys(cols[n.lower()] for n in names if n.lower() in cols))  # case-insensitive, de-duplicated
        if len(hit) == 1:
            out[canon] = hit[0]
        elif len(hit) == 0:
            problems.append(f"{canon}: none of {names} found")
        else:
            problems.append(f"{canon}: ambiguous {hit}")
    if problems:
        raise PrepError(f"cannot resolve columns for {context}: {problems}. Observed header: {list(df_columns)}")
    return out


def strip_version(gene_id: str) -> str:
    return re.sub(r"\.\d+$", "", str(gene_id))


KEY_ALIASES = {
    "snp_id": ["SNP_ID", "snp_id", "SNP", "snp", "rsid"],
    "gene_id": ["GENE_ID", "gene_id", "gene", "TargetID"],
    "pvalue": ["pvalue", "p-value", "p_value", "pval", "P"],
    "rho": ["rho", "spearman_rho", "spearman", "rs", "corr", "correlation"],
}


def parse_answer_key(path) -> pd.DataFrame:
    """Parse a published GEUVADIS cis-eQTL result file (EUR373/YRI89 *.cis.FDR5.{all,best}.rs137.txt.gz)."""
    df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)
    m = resolve_columns(df.columns, KEY_ALIASES, f"answer key {path}")
    out = pd.DataFrame({
        "snp_id": df[m["snp_id"]], "gene_id": df[m["gene_id"]].map(strip_version),
        "pvalue": pd.to_numeric(df[m["pvalue"]], errors="coerce"),
        "rho": pd.to_numeric(df[m["rho"]], errors="coerce"),
    })
    if out[["pvalue", "rho"]].isna().any().any():
        raise PrepError(f"non-numeric pvalue/rho entries in {path}")
    return out


EXPR_META_ALIASES = {
    "gene_id": ["Gene_Symbol", "TargetID", "Gene", "gene_id"],
    "chrom": ["Chr", "chrom", "chromosome"],
    "tss": ["Coord", "TSS", "tss", "Start", "start"],
}


def parse_expression_matrix(path, sample_ids: Optional[Iterable[str]] = None,
                            gene_col: Optional[str] = None, chrom_col: Optional[str] = None,
                            tss_col: Optional[str] = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse the PEER-residual expression matrix.

    Returns (samples x genes expression, gene_meta indexed by gene_id with chrom, tss).
    Metadata column names may be given explicitly; otherwise resolved by alias and the script
    FAILS if ambiguous. Whether ``Coord`` is the TSS (strand-aware) is UNVERIFIED and must be
    checked against the GEUVADIS README on first real use.
    """
    df = pd.read_csv(path, sep="\t", dtype={0: str})
    if gene_col or chrom_col or tss_col:
        m = {"gene_id": gene_col, "chrom": chrom_col, "tss": tss_col}
        bad = [k for k, v in m.items() if v is None or v not in df.columns]
        if bad:
            raise PrepError(f"specified metadata columns missing/unspecified: {bad}; header: {list(df.columns)[:12]}")
    else:
        m = resolve_columns(df.columns[:6], EXPR_META_ALIASES, "expression matrix metadata")
    meta_cols = set(m.values())
    sample_cols = [c for c in df.columns if c not in meta_cols and c != "Gene_Symbol"]
    if sample_ids is not None:
        want = list(sample_ids)
        absent = [s for s in want if s not in sample_cols]
        if absent:
            raise PrepError(f"{len(absent)} requested samples not in expression matrix, e.g. {absent[:5]}")
        sample_cols = want
    gid = df[m["gene_id"]].map(strip_version)
    if gid.duplicated().any():
        raise PrepError("duplicate gene IDs in expression matrix")
    expr = pd.DataFrame(df[sample_cols].to_numpy(dtype=float).T, index=sample_cols, columns=gid.to_numpy())
    meta = pd.DataFrame({"chrom": df[m["chrom"]].astype(str).str.replace("chr", "", regex=False).to_numpy(),
                         "tss": pd.to_numeric(df[m["tss"]]).astype(np.int64).to_numpy()}, index=gid.to_numpy())
    meta.index.name = "gene_id"
    return expr, meta


def read_id_converter(path, from_col: int = 0, to_col: int = 1) -> pd.DataFrame:
    """Read the GEUVADIS dbSNP137 ID converter. Column order is UNVERIFIED: defaults to the first two
    columns (old, new); inspect the file head and override --converter-cols if different."""
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python", dtype=str, comment="#")
    if df.shape[1] < 2:
        raise PrepError("ID converter has fewer than 2 columns")
    return df[[from_col, to_col]].rename(columns={from_col: "old", to_col: "new"})


def read_panel(path) -> pd.DataFrame:
    """1000 Genomes sample panel (sample, pop, super_pop, gender)."""
    df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)
    m = resolve_columns(df.columns, {"sample": ["sample", "Sample"], "pop": ["pop", "Population"],
                                     "super_pop": ["super_pop", "super_population", "Superpopulation"],
                                     "gender": ["gender", "sex", "Sex"]}, "1000G panel")
    out = df[[m["sample"], m["pop"], m["super_pop"], m["gender"]]].copy()
    out.columns = ["sample", "pop", "super_pop", "gender"]
    return out


EUR_POPS = {"CEU", "FIN", "GBR", "TSI"}
YRI_POPS = {"YRI"}


def geuvadis_population(pop_code: str) -> Optional[str]:
    if pop_code in EUR_POPS:
        return "EUR"
    if pop_code in YRI_POPS:
        return "YRI"
    return None


def read_pvar(path) -> pd.DataFrame:
    """Read a PLINK2 .pvar (skips '##' meta lines; header starts with '#CHROM')."""
    with open(path) as fh:
        n_meta = 0
        for line in fh:
            if line.startswith("##"):
                n_meta += 1
            else:
                break
    df = pd.read_csv(path, sep="\t", skiprows=n_meta, dtype=str)
    df.columns = [c.lstrip("#") for c in df.columns]
    need = {"CHROM", "POS", "ID", "REF", "ALT"}
    if not need <= set(df.columns):
        raise PrepError(f"unexpected .pvar header in {path}: {list(df.columns)}")
    return df
