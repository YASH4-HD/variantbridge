import numpy as np
import pandas as pd
import pytest

from variantbridge import population_structure as ps, validation, demo, capability


@pytest.fixture(scope="module")
def sim():
    return demo.simulate_structured_genotypes((60, 60, 60), m=3000, fst=0.08, seed=0)


def test_pca_is_deterministic_and_orthogonal(sim):
    G, lab = sim
    p1, r1, _ = ps.genotype_pca(G, 5)
    p2, r2, _ = ps.genotype_pca(G, 5)
    assert np.allclose(p1, p2) and np.allclose(r1, r2)
    assert (np.diff(r1) <= 1e-12).all() and r1.sum() <= 1.0
    gram = p1.T @ p1
    off = gram - np.diag(np.diag(gram))
    assert np.abs(off).max() < 1e-6 * np.abs(np.diag(gram)).max()


def test_v1_recovers_simulated_population_structure_implementation(sim):
    """Implementation test on a SIMULATED fixture: shows the V1 function works. NOT validation evidence."""
    G, lab = sim
    pcs, _, _ = ps.genotype_pca(G, 5)
    ids = [f"s{i}" for i in range(len(lab))]
    res = validation.v1_population_structure(ps.pcs_frame(pcs, ids), pd.Series(lab, index=ids))
    assert res["knn_cv_balanced_accuracy"] > 0.9 and res["status"] == capability.PASSED
    assert res["eta_squared_by_pc"]["PC1"] > 0.5


def test_v1_fails_when_labels_unrelated_to_pcs(sim):
    G, lab = sim
    pcs, _, _ = ps.genotype_pca(G, 5)
    ids = [f"s{i}" for i in range(len(lab))]
    shuffled = pd.Series(np.random.default_rng(0).permutation(lab), index=ids)
    res = validation.v1_population_structure(ps.pcs_frame(pcs, ids), shuffled)
    assert res["knn_cv_balanced_accuracy"] < 0.6 and res["status"] == capability.FAILED


def test_ld_prune_and_region_mask():
    r = np.random.default_rng(0)
    base = r.binomial(2, 0.4, (200, 1)).astype(float)
    G = np.hstack([base, base, r.binomial(2, 0.4, (200, 3)).astype(float)])   # col 1 duplicates col 0
    keep = ps.ld_prune_greedy(G, 50, 0.2)
    assert keep[0] and not keep[1] and keep[2:].all()
    regions = pd.DataFrame({"chrom": ["6"], "start": [25_000_000], "end": [34_000_000]})
    inside = ps.regions_to_mask(["6", "6", "7", "chr6"], [30_000_000, 40_000_000, 30_000_000, 26_000_000], regions)
    assert inside.tolist() == [True, False, False, True]
