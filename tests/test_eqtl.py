import numpy as np
import pandas as pd
import pytest

from variantbridge import eqtl
from variantbridge.io import GenotypeBlock
from variantbridge.stats_utils import genomic_inflation_lambda


def make_dataset(n=150, n_genes=40, n_var_per_gene=60, planted=(0, 1, 2, 3, 4), effect=0.8, seed=0):
    r = np.random.default_rng(seed)
    samples = np.array([f"S{i:03d}" for i in range(n)])
    pos, vids, gene_meta, y = [], [], {}, {}
    G_all = []
    for g in range(n_genes):
        tss = 5_000_000 + g * 3_000_000          # genes are far apart: independent cis windows
        p0 = tss - 400_000 + np.arange(n_var_per_gene) * 13_000
        maf = r.uniform(0.1, 0.45, n_var_per_gene)
        G = r.binomial(2, maf, (n, n_var_per_gene)).astype(np.float32)
        pos.append(p0); vids += [f"rs{g}_{j}" for j in range(n_var_per_gene)]; G_all.append(G)
        gid = f"ENSG{g:05d}"
        gene_meta[gid] = {"chrom": "22", "tss": tss}
        yy = r.normal(size=n)
        if g in planted:
            yy = yy + effect * (G[:, 10] - G[:, 10].mean())
        y[gid] = yy
    pos = np.concatenate(pos); order = np.argsort(pos, kind="stable")
    G_all = np.hstack(G_all)[:, order]
    block = GenotypeBlock("22", pos[order], np.array(vids)[order], np.array(["G"] * len(pos)), np.array(["A"] * len(pos)), G_all, samples)
    pheno = pd.DataFrame(y, index=samples)
    gm = pd.DataFrame(gene_meta).T
    gm["tss"] = gm["tss"].astype(int)
    cov = pd.DataFrame({"PC1": r.normal(size=n)}, index=samples)
    return pheno, gm, {"22": block}, cov


def test_cis_slice():
    pos = np.array([1, 5, 10, 20, 30])
    assert eqtl.cis_slice(pos, 15, 5) == slice(2, 4)
    assert eqtl.cis_slice(pos, 100, 5) == slice(5, 5)


def test_planted_eqtls_recovered_with_permutation_fdr():
    pheno, gm, blocks, cov = make_dataset(effect=1.3)
    best, pairs = eqtl.run_cis_scan(pheno, gm, blocks, cov, "TEST", n_perm=500, seed=1)
    assert len(best) == 40
    hits = set(best.loc[best.q_bh < 0.05, "gene_id"])
    planted = {f"ENSG{g:05d}" for g in range(5)}
    assert planted <= hits
    assert len(hits - planted) <= 3
    top = best.set_index("gene_id").loc["ENSG00000"]
    assert top.best_variant == "rs0_10" and top.effect > 0.4
    assert (best["empirical_p"] >= 1 / 501).all()   # resolution floor of the +1 permutation p


def test_pair_retention_and_sink_and_consistency():
    from variantbridge.association import check_statistical_consistency
    pheno, gm, blocks, cov = make_dataset(n_genes=3, planted=(0,))
    sunk = []
    best, pairs = eqtl.run_cis_scan(pheno, gm, blocks, cov, "T", keep_pairs={("ENSG00000", "rs0_10")}, pair_sink=sunk.append)
    assert len(pairs) == 1 and pairs.iloc[0]["variant"] == "rs0_10"
    assert len(sunk) == 3
    assert check_statistical_consistency(pd.concat(sunk)) == []


def test_missing_genotypes_not_silently_imputed():
    pheno, gm, blocks, cov = make_dataset(n_genes=2, planted=())
    blocks["22"].dosage[0, 5] = np.nan
    with pytest.raises(ValueError, match="allow_mean_impute"):
        eqtl.run_cis_scan(pheno, gm, blocks, cov, "T")
    best, _ = eqtl.run_cis_scan(pheno, gm, blocks, cov, "T", allow_mean_impute=True)
    assert len(best) == 2


def test_sample_order_mismatch_raises():
    pheno, gm, blocks, cov = make_dataset(n_genes=2, planted=())
    with pytest.raises(ValueError, match="sample order"):
        eqtl.run_cis_scan(pheno.iloc[::-1], gm, blocks, cov, "T")


def test_permuted_phenotype_lambda_near_one_implementation():
    pheno, gm, blocks, cov = make_dataset(n_genes=30, planted=(0, 1, 2, 3, 4))
    lam = eqtl.permuted_scan_lambda(pheno, gm, blocks, cov, n_permutations=3, seed=2)
    assert len(lam) == 3 and 0.8 < np.median(lam) < 1.2
    real = []
    eqtl.run_cis_scan(pheno, gm, blocks, cov, "T", pair_sink=lambda d: real.append(d["p"].to_numpy()))
    assert genomic_inflation_lambda(np.concatenate(real)) > np.median(lam)   # true effects inflate the real scan


def test_maf_either_population_and_covariates():
    r = np.random.default_rng(0)
    def blk(freqs, n=200):
        d = np.column_stack([r.binomial(2, f, n) for f in freqs]).astype(np.float32)
        m = d.shape[1]
        return GenotypeBlock("1", np.arange(m), np.array([f"v{i}" for i in range(m)]), np.array(["G"] * m), np.array(["A"] * m), d, np.array([f"s{i}" for i in range(n)]))
    eur = blk([0.30, 0.01, 0.02]); yri = blk([0.30, 0.20, 0.01])
    mask = eqtl.maf_either_population([eur, yri], 0.05)
    assert mask.tolist() == [True, True, False]
    pcs = pd.DataFrame({"sample_id": ["a", "b", "c"], "PC1": [1., 2, 3], "PC2": [0., 1, 0], "PC3": [5., 5, 6]})
    cov = eqtl.build_covariates(["c", "a"], pcs, 2, pd.Series({"a": 0, "b": 1, "c": 1}))
    assert list(cov.columns) == ["PC1", "PC2", "imputed"] and cov.loc["c", "imputed"] == 1.0
