import copy

import numpy as np
import pandas as pd
import pytest

from variantbridge.amr import cohorts, config, data_io, structure, synthetic


@pytest.fixture(scope="session")
def cfg():
    return config.load_config()


@pytest.fixture(scope="session")
def answer_key():
    return config.load_answer_key()


@pytest.fixture(scope="session")
def syn(cfg):
    """SYNTHETIC cohorts through the same code path as real data (records -> parse -> tiers -> resolve)."""
    s = synthetic.make_synthetic(seed=0)
    rec = data_io.fetch_phenotypes(cfg, records=s["records"])
    kept, excl = cohorts.resolve_cohorts(rec, cfg)
    s.update(parsed=rec, kept=kept, exclusions=excl)
    return s


@pytest.fixture(scope="session")
def disc(syn, cfg):
    d = syn["kept"][syn["kept"]["cohort_tier"] == "discovery"].reset_index(drop=True)
    G = syn["features"].loc[d["genome_id"]].to_numpy()
    st = syn["st"].loc[d["genome_id"]].reset_index(drop=True)
    L, names = structure.lineage_design(st, cfg)
    return {"pheno": d, "G": G, "st": st, "L": L, "ids": list(syn["features"].columns)}


@pytest.fixture(scope="session")
def assoc(disc, cfg, answer_key):
    from variantbridge.amr import association as A
    return A.run_association(disc["G"], disc["ids"], "gene_pa", disc["pheno"], disc["L"], cfg, answer_key[1])


@pytest.fixture
def verified_cfg(cfg):
    """TEST-ONLY config with a VERIFIED breakpoint block holding ARBITRARY test numbers. These numbers are not
    a claim about any real breakpoint; they exercise the interval logic only."""
    c = copy.deepcopy(cfg)
    for k in ("ciprofloxacin_ecoli_discovery", "ciprofloxacin_ecoli_replication"):
        c["breakpoints"][k].update(status="VERIFIED", susceptible_max_mg_l=0.5, resistant_min_mg_l=2.0,
                                   verified_by="unit-test", source_excerpt="TEST FIXTURE ONLY")
    return c
