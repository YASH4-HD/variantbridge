"""SYNTHETIC data generator for deterministic tests and plumbing demos. EVERYTHING produced here is
simulated: it carries data_scope='synthetic' wherever it is persisted and must never be presented as
evidence about real resistance biology. Planted determinants prove only that the software recovers
what was planted."""
from __future__ import annotations

import numpy as np
import pandas as pd

DISC, REP, EXPL, EXPL_SUBSET = "38052776", "34485958", "38219757", "28720578"


def _mic_strings(log2_true: np.ndarray, cohort: str) -> list[str]:
    """Observed MIC strings on a two-fold grid with cohort-specific censoring (mimics the addendum's
    observed ranges: discovery '<=0.25' ... '>4'; replication '<=0.25' ... '>=4'). Synthetic only."""
    out = []
    lo_b, hi_b = -2, 2
    for v in log2_true:
        r = int(np.round(v))
        if r <= lo_b:
            out.append("<=0.25")
        elif r > hi_b and cohort == "discovery":
            out.append(">4")
        elif r >= hi_b and cohort == "replication":
            out.append(">=4")
        elif r > hi_b:
            out.append(">4")
        else:
            out.append(f"{2.0 ** r:g}")
    return out


def make_synthetic(seed: int = 0, n_disc: int = 500, n_rep: int = 250, n_expl: int = 150, n_expl_dup: int = 60,
                   n_lineages_disc: int = 25, n_lineages_rep: int = 12, n_genes: int = 100) -> dict:
    rng = np.random.default_rng(seed)
    feat_ids = [f"gene_{i:03d}" for i in range(n_genes)]
    planted = ["gyrA_S83L", "gyrA_D87N", "parC_S80I", "acrR_Q15*", "qnrS1", "lineage_marker_X"]
    all_ids = feat_ids + planted
    effects = {"gyrA_S83L": 3.0, "gyrA_D87N": 1.2, "parC_S80I": 1.6, "acrR_Q15*": 0.7, "qnrS1": 1.0, "lineage_marker_X": 0.0}

    def block(n, n_lin, cohort, prefix, lineage_prefix):
        sts = np.array([f"{lineage_prefix}{i}" for i in range(n_lin)])
        w = rng.dirichlet(np.ones(n_lin) * 0.6)
        st = rng.choice(sts, size=n, p=w)
        p_gene = {s: np.clip(rng.beta(0.3, 0.3, size=n_genes), 0.02, 0.98) for s in sts}
        G = np.zeros((n, len(all_ids)))
        for i, s in enumerate(st):
            G[i, :n_genes] = rng.random(n_genes) < p_gene[s]
        # resistance determinants: lineage-dependent prevalence (confounded) but variable within lineage
        prev = {k: {s: float(np.clip(rng.beta(0.6, 0.9), 0.03, 0.97)) for s in sts} for k in planted[:5]}
        for j, k in enumerate(planted[:5]):
            col = n_genes + j
            G[:, col] = [rng.random() < prev[k][s] for s in st]
        # confounder: marker present in ~90% of one lineage (which has +1 log2 baseline) and ~3% elsewhere
        hi_lin = sts[0]
        G[:, n_genes + 5] = [(rng.random() < (0.9 if s == hi_lin else 0.03)) for s in st]
        lin_eff = np.where(st == hi_lin, 1.0, 0.0) + rng.normal(0, 0.3, size=n_lin)[pd.Index(sts).get_indexer(st)]
        mu = -3.0 + lin_eff + sum(effects[k] * G[:, n_genes + j] for j, k in enumerate(planted))
        y = mu + rng.normal(0, 0.8, size=n)
        ids = [f"{prefix}.{i:05d}" for i in range(n)]
        return ids, st, G, y

    d_ids, d_st, d_G, d_y = block(n_disc, n_lineages_disc, "discovery", "562.1", "ST_D")
    r_ids, r_st, r_G, r_y = block(n_rep, n_lineages_rep, "replication", "562.2", "ST_R")
    e_ids, e_st, e_G, e_y = block(n_expl, 15, "exploratory", "562.3", "ST_E")
    G = np.vstack([d_G, r_G, e_G]); ids = d_ids + r_ids + e_ids; st = np.concatenate([d_st, r_st, e_st])
    cohort_of = ["discovery"] * n_disc + ["replication"] * n_rep + ["exploratory"] * n_expl
    ys = np.concatenate([d_y, r_y, e_y])
    recs = []
    for gid, y, c in zip(ids, ys, cohort_of):
        pm = {"discovery": DISC, "replication": REP, "exploratory": EXPL}[c]
        recs.append({"genome_id": gid, "study_pmid": pm, "laboratory_typing_method": "MIC", "_log2": y, "_c": c})
    rec = pd.DataFrame(recs)
    rec["measurement"] = ""
    for c in ("discovery", "replication", "exploratory"):
        m = rec["_c"] == c
        rec.loc[m, "measurement"] = _mic_strings(rec.loc[m, "_log2"].to_numpy(), c if c != "exploratory" else "exploratory")
    rec["resistant_phenotype"] = np.nan  # calls are NOT used (breakpoints unverified)
    # records that exercise duplicate logic (all synthetic):
    extra = []
    # (a) exploratory 28720578 records duplicating 38219757 genomes (superset rule: drop the 28720578 record)
    for gid in e_ids[:n_expl_dup]:
        row = rec[rec["genome_id"] == gid].iloc[0].copy(); row["study_pmid"] = EXPL_SUBSET; extra.append(row)
    # (b) a discovery genome also appearing in a lower-tier study (must stay discovery only)
    row = rec[rec["genome_id"] == d_ids[0]].iloc[0].copy(); row["study_pmid"] = EXPL; extra.append(row)
    # (c) exact duplicate record within discovery
    extra.append(rec[rec["genome_id"] == d_ids[1]].iloc[0].copy())
    # (d) conflicting duplicate within discovery (different MIC for same genome)
    row = rec[rec["genome_id"] == d_ids[2]].iloc[0].copy(); row["measurement"] = ">4" if rec.loc[rec["genome_id"] == d_ids[2], "measurement"].iloc[0] != ">4" else "<=0.25"; extra.append(row)
    # (e) a non-MIC (disk diffusion) record for an exploratory genome with no MIC
    row = rec[rec["genome_id"] == e_ids[-1]].iloc[0].copy(); row["genome_id"] = "562.3.99999"; row["laboratory_typing_method"] = "Disk diffusion"; row["measurement"] = "22"; extra.append(row)
    rec = pd.concat([rec, pd.DataFrame(extra)], ignore_index=True).drop(columns=["_log2", "_c"])
    # genome metadata (QC); includes a few deliberately failing genomes
    uids = sorted(set(ids) | {"562.3.99999"})
    meta = pd.DataFrame({"genome_id": uids, "genome_name": [f"Escherichia coli strain S{i}" for i in range(len(uids))],
                         "checkm_completeness": 99.0, "checkm_contamination": 0.5, "genome_length": 5_000_000, "contigs": 120})
    meta.loc[3, "checkm_completeness"] = 80.0
    meta.loc[4, "checkm_contamination"] = 7.0
    meta.loc[5, "genome_length"] = 3_000_000
    meta.loc[6, "genome_name"] = "Escherichia coli K-12 MG1655"
    meta.loc[7, "contigs"] = 900
    meta.loc[8, "checkm_completeness"] = np.nan
    mlst = dict(zip(ids, st)); meta["mlst"] = meta["genome_id"].map(mlst)
    fmeta = pd.DataFrame({"feature_id": all_ids,
                          "feature_type": ["gene_pa"] * n_genes + ["amino_acid_substitution"] * 3 + ["regulatory_mutation", "gene_pa", "gene_pa"],
                          "gene": [None] * n_genes + ["gyrA", "gyrA", "parC", "acrR", "qnrS1", None],
                          "aa_ref": [None] * n_genes + ["S", "D", "S", "Q", None, None],
                          "aa_pos": [np.nan] * n_genes + [83, 87, 80, 15, np.nan, np.nan],
                          "aa_alt": [None] * n_genes + ["L", "N", "I", "*", None, None]})
    return {"records": rec, "meta": meta, "features": pd.DataFrame(G, index=ids, columns=all_ids), "feature_meta": fmeta,
            "st": pd.Series(st, index=ids), "planted_effects": effects, "data_scope": "synthetic"}
