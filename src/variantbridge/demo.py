"""SIMULATED datasets. Everything produced here is synthetic and must be labelled SIMULATED in
any UI or report. Simulated data are never validation evidence (handoff rule 2)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def human_demo(seed=7, n=320, m=180):
    r = np.random.default_rng(seed)
    ids = [f"H{i:04d}" for i in range(n)]
    ch = r.integers(1, 23, m)
    pos = r.integers(1_000_000, 180_000_000, m)
    v = [f"chr{a}:{b}:A:G" for a, b in zip(ch, pos)]
    maf = r.uniform(.05, .48, m)
    anc = r.integers(0, 2, n)
    G = np.zeros((n, m), float)
    for j, p in enumerate(maf):
        pp = np.clip(p + (anc - .5) * (.18 if j < 25 else .02), .01, .95)
        G[:, j] = r.binomial(2, pp)
    G[r.random(G.shape) < .008] = np.nan
    F = np.where(np.isnan(G), np.nanmean(G, 0), G)
    age = r.normal(48, 13, n).clip(18, 85)
    sex = r.integers(0, 2, n)
    logit = -2.1 + .035 * (age - 45) + .25 * sex + .65 * F[:, 5] - .55 * F[:, 17] + .45 * F[:, 49]
    y = r.binomial(1, 1 / (1 + np.exp(-logit)))
    ge = pd.DataFrame(G, columns=v)
    ge.insert(0, "sample_id", ids)
    ph = pd.DataFrame({"sample_id": ids, "phenotype": y, "age": age.round(1), "sex": sex})
    me = pd.DataFrame({"variant": v, "chrom": ch, "pos": pos})
    return ge, ph, me


def amr_demo(seed=11, n=300, m=140):
    r = np.random.default_rng(seed)
    ids = [f"ISO{i:04d}" for i in range(n)]
    lin = r.integers(0, 3, n)
    g = [f"gene_{i:03d}" for i in range(m)]
    X = np.zeros((n, m))
    for j, p in enumerate(r.uniform(.05, .5, m)):
        X[:, j] = r.binomial(1, np.clip(p + (lin - 1) * r.uniform(-.10, .10), .01, .95))
    causal = [2, 12, 37, 75]
    z = -1.7 + 1.15 * X[:, 2] + .9 * X[:, 12] - .75 * X[:, 37] + .7 * X[:, 75] + .3 * (lin == 2)
    y = r.binomial(1, 1 / (1 + np.exp(-z)))
    ge = pd.DataFrame(X, columns=g)
    ge.insert(0, "sample_id", ids)
    ph = pd.DataFrame({"sample_id": ids, "phenotype": y, "lineage": lin})
    me = pd.DataFrame({"variant": g, "chrom": "chromosome", "pos": np.arange(1, m + 1) * 10000,
                       "known_amr": ["Known" if i in causal[:2] else "Uncertain" for i in range(m)]})
    return ge, ph, me


def simulate_structured_genotypes(n_per_pop=(60, 60, 60), m=2000, fst=0.08, seed=0):
    """Balding-Nichols simulation of populations with drift. Returns (G, labels). Test fixture only."""
    r = np.random.default_rng(seed)
    anc = r.uniform(0.1, 0.9, m)
    blocks, labs = [], []
    for k, n in enumerate(n_per_pop):
        a = anc * (1 - fst) / fst
        b = (1 - anc) * (1 - fst) / fst
        p = r.beta(a, b)
        blocks.append(r.binomial(2, p, (n, m)).astype(float))
        labs += [f"POP{k + 1}"] * n
    return np.vstack(blocks), np.array(labs)


def simulate_discovery_target(n_disc=1500, n_target=500, m=300, h2_snps=12, seed=3):
    """Simulated quantitative-trait cohorts with disjoint samples for the PRS method demonstration."""
    r = np.random.default_rng(seed)
    maf = r.uniform(0.1, 0.5, m)
    causal = r.choice(m, h2_snps, replace=False)
    beta = np.zeros(m)
    beta[causal] = r.normal(0, 0.25, h2_snps)

    def make(n, prefix):
        G = r.binomial(2, maf, (n, m)).astype(float)
        y = G @ beta + r.normal(0, 1.0, n)
        return pd.DataFrame(G, columns=[f"snp{j}" for j in range(m)],
                            index=[f"{prefix}{i:05d}" for i in range(n)]), pd.Series(y, index=[f"{prefix}{i:05d}" for i in range(n)])

    Gd, yd = make(n_disc, "D")
    Gt, yt = make(n_target, "T")
    return Gd, yd, Gt, yt, beta
