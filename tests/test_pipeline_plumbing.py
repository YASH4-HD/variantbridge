"""End-to-end PLUMBING tests: scripts 03 -> 05 -> 06 on SIMULATED GEUVADIS-like files using a real
plink2 binary. These prove the scripts run and that results agree with an independent statsmodels
re-computation under the ASSUMED file formats. They are NOT V4 (no published answer key is used) and
produce no biological or replication claim."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from scipy import stats

from synthetic import write_geuvadis_like

ROOT = Path(__file__).resolve().parents[1]


def run(script, *args, env_plink2):
    cp = subprocess.run([sys.executable, str(ROOT / "scripts" / script), *args], capture_output=True, text=True,
                        env={**__import__("os").environ, "PLINK2": env_plink2}, cwd=ROOT)
    assert cp.returncode == 0, cp.stdout[-2000:] + cp.stderr[-3000:]
    return cp.stdout


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, plink2_bin):
    d = tmp_path_factory.mktemp("geu")
    info = write_geuvadis_like(d)
    wd, ad = d / "work", d / "artifacts"
    out = run("03_prepare_geuvadis.py", "--vcf-template", str(d / "geu.chr{chrom}.vcf"), "--chroms", "22",
              "--panel", str(d / "panel.txt"), "--expression", str(d / "expression.txt"),
              "--imputation-status-file", str(d / "imputed.tsv"), "--id-converter", str(d / "idconvert.txt"),
              "--no-long-range-exclusion", "--work-dir", str(wd), env_plink2=plink2_bin)
    return d, wd, ad, info, out


def test_prepare_outputs_and_deviations_recorded(prepared):
    d, wd, ad, info, out = prepared
    man = json.loads((wd / "geuvadis_prep.manifest.json").read_text())
    assert any("long-range" in x for x in man["inclusion_exclusion"]["deviations"])
    assert any("90%" in x for x in man["inclusion_exclusion"]["deviations"])
    cov = pd.read_csv(wd / "covariates_EUR.csv", index_col=0)
    assert list(cov.columns) == ["PC1", "PC2", "PC3", "imputed"]
    assert list(pd.read_csv(wd / "covariates_YRI.csv", index_col=0).columns) == ["PC1", "PC2", "imputed"]
    expr = pd.read_csv(wd / "expression_EUR.csv.gz", index_col=0)
    assert expr.columns[0] == "ENSG00000"                       # version suffix stripped
    from variantbridge.io import GenotypeBlock
    blk = GenotypeBlock.load(wd / "blocks" / "chr22.npz")
    assert blk.variant_ids[0].startswith("rs")                  # dbSNP-style IDs applied via the converter
    assert set(blk.counted_allele) == {"G"} and set(blk.other_allele) == {"A"}   # ALT counted explicitly (plink2 default would count REF)
    assert np.all(np.diff(blk.positions) >= 0)


def test_scan_matches_independent_statsmodels_and_recovers_planted(prepared, plink2_bin):
    d, wd, ad, info, _ = prepared
    run("05_run_eqtl_scan.py", "--work-dir", str(wd), "--artifact-dir", str(ad), "--n-perm", "300", env_plink2=plink2_bin)
    best = pd.read_csv(ad / "geuvadis_eqtl_EUR.csv.gz")
    cov = pd.read_csv(wd / "covariates_EUR.csv", index_col=0)
    expr = pd.read_csv(wd / "expression_EUR.csv.gz", index_col=0)
    from variantbridge.io import GenotypeBlock
    blk = GenotypeBlock.load(wd / "blocks" / "chr22.npz").subset_samples(cov.index.to_numpy())
    row = best.set_index("gene_id").loc["ENSG00000"]
    j = np.flatnonzero(blk.variant_ids == row["best_variant"])[0]
    f = sm.OLS(expr["ENSG00000"].to_numpy(), sm.add_constant(np.column_stack([cov.to_numpy(), blk.dosage[:, j].astype(float)]))).fit()
    assert row["effect"] == pytest.approx(f.params[-1], rel=1e-6) and row["p_nominal"] == pytest.approx(f.pvalues[-1], rel=1e-5)
    assert row["effect"] > 0.5      # planted positive ALT-dosage effect has the right sign (orientation test)
    hits = set(best.loc[best.q_bh < 0.05, "gene_id"])
    assert {f"ENSG{g:05d}" for g in info["planted"][:4]} <= hits
    man = json.loads((ad / "geuvadis_eqtl_EUR.manifest.json").read_text())
    from variantbridge import manifest as mf
    assert mf.validate_manifest(man) == [] and man["output_sha256"] == mf.file_sha256(ad / "geuvadis_eqtl_EUR.csv.gz")


def test_v4_script_plumbing_with_independently_built_key(prepared, plink2_bin):
    d, wd, ad, info, _ = prepared
    cov = pd.read_csv(wd / "covariates_EUR.csv", index_col=0)
    expr = pd.read_csv(wd / "expression_EUR.csv.gz", index_col=0)
    from variantbridge.io import GenotypeBlock
    blk = GenotypeBlock.load(wd / "blocks" / "chr22.npz").subset_samples(cov.index.to_numpy())
    gmeta = pd.read_csv(wd / "gene_meta.csv", index_col=0)
    rows = []
    for gene in expr.columns:
        tss = gmeta.at[gene, "tss"]
        idx = np.flatnonzero(np.abs(blk.positions - tss) <= 1_000_000)
        best = None
        for j in idx:
            x = blk.dosage[:, j].astype(float)
            if x.std() == 0:
                continue
            f = sm.OLS(expr[gene].to_numpy(), sm.add_constant(np.column_stack([cov.to_numpy(), x]))).fit()
            if best is None or f.pvalues[-1] < best[1]:
                best = (blk.variant_ids[j], f.pvalues[-1], stats.spearmanr(x, expr[gene].to_numpy())[0])
        rows.append((best[0], gene, best[1], best[2], 0))
    key = pd.DataFrame(rows, columns=["SNP_ID", "GENE_ID", "pvalue", "rho", "distance"])
    key = key[key["pvalue"] < 1e-3]            # published-like: only significant pairs
    key.to_csv(d / "key_best.txt.gz", sep=" ", index=False)
    key.to_csv(d / "key_all.txt.gz", sep=" ", index=False)
    run("05_run_eqtl_scan.py", "--work-dir", str(wd), "--artifact-dir", str(ad), "--n-perm", "100",
        "--key-files", f"EUR={d / 'key_best.txt.gz'}", env_plink2=plink2_bin)
    out = run("06_validate_answer_key_v4.py", "--population", "EUR", "--key-best", str(d / "key_best.txt.gz"),
              "--key-all", str(d / "key_all.txt.gz"), "--work-dir", str(wd), "--artifact-dir", str(ad), env_plink2=plink2_bin)
    res = json.loads((ad / "v4_EUR.json").read_text())
    assert res["same_direction_rate"] == 1.0 and res["n_published_best_pairs_compared"] == len(key)
    assert res["spearman_neglog10p"] > 0.99
    # The fixture has far fewer than the pre-registered minimum of 100 compared pairs, so the verdict must be
    # FAILED (insufficient evidence) rather than PASSED, and YRI was never run: V4 can never be PASSED here.
    st = json.loads((ad / "VALIDATION_STATUS.json").read_text())["V4"]
    assert res["status"] == "FAILED" and res["n_published_best_pairs_compared"] < 100
    assert st["status"] != "PASSED" and set(st["populations"]) == {"EUR"}


def test_script_refuses_silent_deviations(tmp_path, plink2_bin):
    cp = subprocess.run([sys.executable, str(ROOT / "scripts" / "03_prepare_geuvadis.py"), "--vcf-template", "x{chrom}",
                         "--panel", "p", "--expression", "e", "--work-dir", str(tmp_path)], capture_output=True, text=True,
                        env={**__import__("os").environ, "PLINK2": plink2_bin})
    assert cp.returncode != 0 and "imputation" in (cp.stdout + cp.stderr)


def test_v6_script_reports_implementation_and_biology_separately(prepared, plink2_bin):
    d, wd, ad, info, _ = prepared
    from variantbridge.io import GenotypeBlock
    blk = GenotypeBlock.load(wd / "blocks" / "chr22.npz")
    rng = np.random.default_rng(0)
    rows = []
    for g in range(24):
        for j in range(70):
            if j == 10:
                z = 6.0 if g in info["planted"] else rng.normal(0, 1)
                p = 1e-9 if g in info["planted"] else 0.5
            else:
                z, p = rng.normal(0, 1), 0.5
            rows.append(dict(Pvalue=p, SNP=f"rs{g}_{j}", AssessedAllele=blk.counted_allele[np.flatnonzero(blk.variant_ids == f"rs{g}_{j}")[0]],
                             OtherAllele="A" if blk.counted_allele[np.flatnonzero(blk.variant_ids == f"rs{g}_{j}")[0]] == "G" else "G",
                             Zscore=z, Gene=f"ENSG{g:05d}"))
    pd.DataFrame(rows).to_csv(d / "weights.txt", sep="\t", index=False)
    common = ["--weights", str(d / "weights.txt"), "--work-dir", str(wd), "--artifact-dir", str(ad), "--n-permutations", "200",
              "--n-bootstrap", "100", "--p-threshold", "1e-4"]
    run("08_run_expression_score_v6.py", *common, "--cohort-disjointness-attested", env_plink2=plink2_bin)
    st = json.loads((ad / "VALIDATION_STATUS.json").read_text())["V6"]
    assert st["implementation"]["status"] == "PASSED" and st["status"] == "PASSED"
    assert st["biological_evidence"]["evidence_above_permutation_null"] is True
    assert "Genetic Expression Score" in st["terminology"] and "PRS" in st["terminology"]   # only to say it is NOT a PRS
    assert "not PRS" in st["terminology"]
    # without the attestation the implementation (and hence build) status is FAILED, independent of the biology
    run("08_run_expression_score_v6.py", *common, "--no-record", env_plink2=plink2_bin)   # no-record keeps previous status
    cp = subprocess.run([sys.executable, str(ROOT / "scripts" / "08_run_expression_score_v6.py"), *common],
                        capture_output=True, text=True, env={**__import__("os").environ, "PLINK2": plink2_bin}, cwd=ROOT)
    st2 = json.loads((ad / "VALIDATION_STATUS.json").read_text())["V6"]
    assert st2["implementation"]["checks"]["cohort_disjointness_attested"] is False and st2["status"] == "FAILED"
    assert st2["biological_evidence"]["evidence_above_permutation_null"] is True       # biology reported independently
