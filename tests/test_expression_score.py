import numpy as np
import pandas as pd
import pytest

from variantbridge import expression_score as es, capability, validation
from variantbridge.io import GenotypeBlock


def build(n=300, n_genes=60, k=8, signal=True, seed=0):
    r = np.random.default_rng(seed)
    samples = np.array([f"S{i:03d}" for i in range(n)])
    rows, G_all, pos, vids, ca, oa = [], [], [], [], [], []
    expr, tss = {}, {}
    for g in range(n_genes):
        gid = f"ENSG{g:05d}"
        t = 1_000_000 + g * 5_000_000
        G = r.binomial(2, r.uniform(0.2, 0.5, k), (n, k)).astype(np.float32)
        z = r.uniform(3, 8, k) * r.choice([-1, 1], k)
        # discovery orientation: half of the weights are given for the OTHER allele (need flipping)
        flip = r.random(k) < 0.5
        for j in range(k):
            vid = f"rs{g}_{j}"
            vids.append(vid); pos.append(t - 200_000 + j * 20_000); ca.append("G"); oa.append("A")
            rows.append(dict(gene_id=gid, variant_id=vid, pvalue=1e-8, zscore=-z[j] if flip[j] else z[j],
                             assessed_allele="A" if flip[j] else "G", other_allele="G" if flip[j] else "A"))
        G_all.append(G)
        truth = G @ (z / 20.0)
        expr[gid] = (truth if signal else 0 * truth) + r.normal(0, 1.0, n)
        tss[gid] = t
    pos = np.array(pos); order = np.argsort(pos)
    block = GenotypeBlock("1", pos[order], np.array(vids)[order], np.array(ca)[order], np.array(oa)[order], np.hstack(G_all)[:, order], samples)
    return pd.DataFrame(rows), block, pd.Series(tss), pd.DataFrame(expr, index=samples)


def test_alignment_flips_and_drops():
    W, block, _, _ = build(n_genes=2, k=4)
    extra = pd.DataFrame([
        dict(gene_id="ENSG00000", variant_id="rs0_0", pvalue=1e-9, zscore=1.0, assessed_allele="C", other_allele="T"),   # mismatch
        dict(gene_id="ENSG00000", variant_id="nope", pvalue=1e-9, zscore=1.0, assessed_allele="G", other_allele="A"),     # absent
    ])
    w, rep = es.align_weights(pd.concat([W, extra], ignore_index=True), block)
    assert rep["dropped_not_in_target"] == 1 and rep["dropped_allele_mismatch"] == 1
    assert rep["n_aligned"] == 8
    # weights given on A (the other allele) must be sign-flipped into the counted-allele (G) frame
    src = W.set_index("variant_id")
    for r in w.itertuples():
        expected = src.loc[r.variant_id, "zscore"] * (1 if src.loc[r.variant_id, "assessed_allele"] == "G" else -1)
        assert r.w == pytest.approx(expected)


def test_palindromic_snps_are_dropped():
    W, block, _, _ = build(n_genes=1, k=3)
    block.counted_allele[:] = "A"; block.other_allele[:] = "T"
    W["assessed_allele"] = "A"; W["other_allele"] = "T"
    w, rep = es.align_weights(W, block)
    assert rep["dropped_palindromic"] == 3 and len(w) == 0


def test_pipeline_recovers_known_signal_implementation():
    """Implementation test on SIMULATED data with a planted relationship: shows the scoring pipeline works."""
    W, block, tss, expr = build(signal=True)
    pred, info, rep = es.predict_expression(W, block, tss, es.ExpressionScoreParams(p_threshold=1e-4, r2_clump=0.9))
    ev = es.evaluate_expression_score(pred, expr, n_permutations=200, n_bootstrap=200, seed=1)
    assert ev["n_genes_scored"] == 60 and ev["observed_mean_r"] > 0.3
    assert ev["empirical_p_one_sided"] == pytest.approx(1 / 201)
    assert ev["bootstrap_ci95_mean_r"][0] > 0.2 and ev["fraction_genes_r_positive"] > 0.9
    assert ev["null_ci95"][1] < 0.1


def test_null_data_gives_null_biological_result_but_implementation_passes():
    """A null biological result must NOT fail the build (amendment 3)."""
    W, block, tss, expr = build(signal=False, seed=3)
    pred, info, rep = es.predict_expression(W, block, tss, es.ExpressionScoreParams(r2_clump=0.9))
    ev = es.evaluate_expression_score(pred, expr, n_permutations=300, n_bootstrap=200, seed=1)
    assert ev["empirical_p_one_sided"] > 0.05
    checks = {"genes_scored_gt_0": ev["n_genes_scored"] > 0, "null_generated": ev["n_permutations"] >= 100}
    res = validation.v6_expression_score({k: v for k, v in ev.items() if k not in ("per_gene", "null_distribution")}, rep, checks)
    assert res["implementation"]["status"] == capability.PASSED and res["status"] == capability.PASSED
    assert res["biological_evidence"]["evidence_above_permutation_null"] is False
    assert "valid null biological result" in res["biological_evidence"]["interpretation"]


def test_implementation_failure_is_reported_separately():
    ev = dict(empirical_p_one_sided=0.001, bootstrap_ci95_mean_r=[0.1, 0.2], observed_mean_r=0.15, observed_mean_r2=0.03,
              n_genes_scored=10, n_samples=100, null_ci95=[-0.01, 0.01], n_permutations=1000)
    res = validation.v6_expression_score(ev, {}, {"cohort_disjointness_attested": False, "genes": True})
    assert res["status"] == capability.FAILED and res["biological_evidence"]["evidence_above_permutation_null"] is True


def test_terminology_never_prs_in_expression_module():
    import inspect
    src = inspect.getsource(es)
    head = src.split('"""')[1]
    assert "NOT a polygenic risk score" in head
    assert not hasattr(es, "prs_score")
