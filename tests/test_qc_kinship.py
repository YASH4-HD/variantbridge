import numpy as np
import pandas as pd
import pytest
from math import lgamma, log, exp

from variantbridge import qc


def brute_hwe(n_aa, n_ab, n_bb):
    """Independent exact HWE p-value by enumerating the conditional distribution."""
    n = n_aa + n_ab + n_bb
    n_a = 2 * n_aa + n_ab                    # count of allele A
    def logp(h):
        a = (n_a - h) // 2
        b = n - h - a
        return (lgamma(n + 1) - lgamma(a + 1) - lgamma(h + 1) - lgamma(b + 1) + h * log(2)
                + lgamma(n_a + 1) + lgamma(2 * n - n_a + 1) - lgamma(2 * n + 1))
    hs = [h for h in range(n_a % 2, min(n_a, 2 * n - n_a) + 1, 2)]
    ps = np.array([exp(logp(h)) for h in hs])
    obs = exp(logp(n_ab))
    return float(ps[ps <= obs * (1 + 1e-9)].sum())


@pytest.mark.parametrize("counts", [(10, 20, 10), (25, 0, 25), (50, 10, 1), (3, 4, 2), (100, 50, 8), (0, 0, 5), (7, 1, 0)])
def test_hwe_exact_matches_bruteforce(counts):
    assert qc.hwe_exact_pvalue(*counts) == pytest.approx(brute_hwe(*counts), rel=1e-6, abs=1e-12)


def test_hwe_equilibrium_not_rejected_and_deviation_rejected():
    assert qc.hwe_exact_pvalue(250, 500, 250) > 0.5
    assert qc.hwe_exact_pvalue(500, 0, 500) == pytest.approx(brute_hwe(500, 0, 500), rel=1e-6, abs=1e-250)
    assert qc.hwe_exact_pvalue(500, 0, 500) < 1e-9


def test_sample_and_variant_metrics_and_thresholds():
    r = np.random.default_rng(0)
    G = r.binomial(2, 0.4, (100, 200)).astype(float)
    G[0, :60] = np.nan                        # bad call rate sample
    sm_ = qc.sample_metrics(G)
    thr = qc.QCThresholds()
    mask = qc.sample_qc_mask(sm_, thr)
    assert not mask.iloc[0] and mask.iloc[1:].mean() > 0.9
    vm = qc.variant_metrics(G)
    assert (vm["maf"] <= 0.5).all()
    assert qc.variant_qc_mask(vm, thr).dtype == bool


def test_sex_check():
    x_het = np.array([0.001, 0.30, 0.002, 0.28])
    inferred = qc.infer_sex_from_x_het(x_het, 0.1)
    assert inferred.tolist() == [1, 2, 1, 2]
    c = qc.sex_concordance(inferred, np.array([1, 2, 1, 1]))
    assert c["n_compared"] == 4 and c["n_concordant"] == 3 and c["concordance"] == 0.75


def simulate_family(m=6000, seed=1):
    r = np.random.default_rng(seed)
    p = r.uniform(0.1, 0.5, m)
    hap = lambda: (r.random(m) < p).astype(int)
    mom = hap() + hap(); dad = hap() + hap()
    def transmit(g):  # one allele from a diploid genotype (0/1/2)
        return np.where(g == 2, 1, np.where(g == 0, 0, (r.random(m) < 0.5).astype(int)))
    child = transmit(mom) + transmit(dad)
    unrel = hap() + hap()
    unrel2 = hap() + hap()
    return np.vstack([mom, dad, child, unrel, unrel2]).astype(float)


def test_king_robust_known_relationships():
    G = simulate_family()
    G = np.vstack([G, G[2]])            # exact duplicate of child (index 5)
    phi = qc.king_robust(G)
    assert phi[2, 5] == pytest.approx(0.5, abs=1e-9)             # duplicate
    assert 0.2 < phi[0, 2] < 0.3 and 0.2 < phi[1, 2] < 0.3        # parent-offspring ~0.25
    assert abs(phi[0, 1]) < 0.03                                   # unrelated parents
    assert abs(phi[3, 4]) < 0.03
    assert qc.classify_kinship(phi[2, 5]) == "duplicate/MZ twin"
    assert qc.classify_kinship(phi[0, 2]) == "1st-degree"
    pairs = qc.related_pairs(phi, [f"s{i}" for i in range(6)])
    assert set(pairs["degree"]) <= {"duplicate/MZ twin", "1st-degree"}
