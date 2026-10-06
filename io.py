"""Input validation and compact genotype-block storage.

Streamlit never reads raw genomic files. Offline scripts write compact ``GenotypeBlock``
files (npz) and tables (csv.gz); the app and the scientific modules consume those.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class GenotypeBlock:
    """Additive dosage for one chromosome. dosage is n_samples x n_variants, NaN = missing."""
    chrom: str
    positions: np.ndarray          # int64, sorted ascending
    variant_ids: np.ndarray        # str
    counted_allele: np.ndarray     # str, allele whose dosage is stored
    other_allele: np.ndarray       # str
    dosage: np.ndarray             # float32
    sample_ids: np.ndarray         # str

    def __post_init__(self):
        # store strings as fixed-width unicode (not object) so npz files load with allow_pickle=False
        for name in ("variant_ids", "counted_allele", "other_allele", "sample_ids"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=str))
        self.positions = np.asarray(self.positions, dtype=np.int64)
        self.dosage = np.asarray(self.dosage, dtype=np.float32)
        n, m = self.dosage.shape
        assert len(self.sample_ids) == n, "sample_ids length mismatch"
        for name in ("positions", "variant_ids", "counted_allele", "other_allele"):
            assert len(getattr(self, name)) == m, f"{name} length mismatch"
        if m > 1 and np.any(np.diff(self.positions) < 0):
            raise ValueError("positions must be sorted ascending")

    def subset_samples(self, sample_ids) -> "GenotypeBlock":
        idx = pd.Index(self.sample_ids).get_indexer(sample_ids)
        if (idx < 0).any():
            raise KeyError("some requested samples are absent from the genotype block")
        return GenotypeBlock(self.chrom, self.positions, self.variant_ids, self.counted_allele,
                             self.other_allele, self.dosage[idx], np.asarray(sample_ids))

    def subset_variants(self, mask) -> "GenotypeBlock":
        mask = np.asarray(mask)
        return GenotypeBlock(self.chrom, self.positions[mask], self.variant_ids[mask],
                             self.counted_allele[mask], self.other_allele[mask],
                             self.dosage[:, mask], self.sample_ids)

    def save(self, path) -> None:
        np.savez_compressed(path, chrom=self.chrom, positions=self.positions, variant_ids=self.variant_ids,
                            counted_allele=self.counted_allele, other_allele=self.other_allele,
                            dosage=self.dosage, sample_ids=self.sample_ids)

    @staticmethod
    def load(path) -> "GenotypeBlock":
        z = np.load(path, allow_pickle=False)
        return GenotypeBlock(str(z["chrom"]), z["positions"], z["variant_ids"], z["counted_allele"],
                             z["other_allele"], z["dosage"], z["sample_ids"])


def read_plink_raw(raw_path, pvar_path=None) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    """Parse a PLINK ``--export A`` .raw file. Returns (sample_ids, dosage, variant_ids, counted_alleles).

    The counted allele is read from the column header (``<id>_<allele>``); it is never assumed.
    """
    df = pd.read_csv(raw_path, sep=r"\s+", dtype=str)
    iid = df["IID"].to_numpy()
    geno_cols = [c for c in df.columns if c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    vids, alleles = [], []
    for c in geno_cols:
        vid, allele = c.rsplit("_", 1)
        vids.append(vid)
        alleles.append(allele)
    dos = df[geno_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    return iid, dos, vids, alleles


def load_csv_inputs(geno: pd.DataFrame, pheno: pd.DataFrame, phenotype_col: str):
    """Validate and align uploaded CSVs (used by the app and tests).

    Returns merged table, genotype matrix X (NaN preserved), phenotype y, feature names.
    Unlike v0.1 this does NOT mean-impute silently: missing values stay NaN and are handled by
    the association engine (complete-case per variant) and reported in QC.
    """
    if "sample_id" not in geno or "sample_id" not in pheno:
        raise ValueError("both files require a 'sample_id' column")
    if geno["sample_id"].duplicated().any() or pheno["sample_id"].duplicated().any():
        raise ValueError("duplicate sample_id values found")
    if phenotype_col not in pheno:
        raise ValueError(f"phenotype column {phenotype_col!r} not found")
    merged = geno.merge(pheno, on="sample_id", how="inner")
    if merged.empty:
        raise ValueError("no overlapping sample_id values between genotype and phenotype files")
    feats = [c for c in geno.columns if c != "sample_id"]
    X = merged[feats].apply(pd.to_numeric, errors="coerce")
    y = pd.to_numeric(merged[phenotype_col], errors="coerce")
    keep = y.notna()
    return (merged.loc[keep].reset_index(drop=True), X.loc[keep].reset_index(drop=True),
            y.loc[keep].reset_index(drop=True), feats)
