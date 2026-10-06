"""V5 (RELEASE-BLOCKING): VariantBridge association vs `plink2 --glm` on an identical SIMULATED dataset/model.
Skipped (V5 stays NOT RUN) when no plink2 binary is available."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_script():
    spec = importlib.util.spec_from_file_location("v5script", ROOT / "scripts" / "07_benchmark_plink2_v5.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.plink2
def test_v5_plink2_glm_benchmark(tmp_path, plink2_bin):
    mod = load_script()
    rec = mod.run_benchmark(tmp_path / "work", plink2_bin, record_root=str(tmp_path / "art"))
    assert rec["status"] == "PASSED", json.dumps({k: v for k, v in rec["results"].items()}, indent=1, default=str)
    for model, r in rec["results"].items():
        assert r["n_variants_matched"] == 400 and r["n_unmatched"] == 0
        assert r["max_rel_diff_effect"] <= 1e-4 and r["max_rel_diff_se"] <= 1e-4 and r["obs_count_identical"]
        assert r["internal_consistency_violations"] == 0
    assert "SIMULATED" in rec["dataset"]
    status = json.loads((tmp_path / "art" / "VALIDATION_STATUS.json").read_text())
    assert status["V5"]["status"] == "PASSED" and "plink2_version" in status["V5"]


@pytest.mark.plink2
def test_v5_real_fileset_mode_plumbing_uses_same_code_path(tmp_path, plink2_bin):
    """Exercises --pfile mode with a SIMULATED fileset standing in for a real chr22 subset (plumbing only)."""
    mod = load_script()
    wd = tmp_path / "w"; wd.mkdir()
    mod.write_simulated_plink_inputs(wd, seed=7)
    from variantbridge.prep import run_plink2
    run_plink2(plink2_bin, ["--vcf", str(wd / "sim.vcf"), "--make-pgen", "--out", str(wd / "sim")])
    rec = mod.run_benchmark(tmp_path / "w2", plink2_bin, record_root=None, real=dict(
        pfile=str(wd / "sim"), pheno_file=str(wd / "pheno.tsv"), pheno_col="y_lin", covar_file=str(wd / "covar.tsv"),
        covar_cols=["age", "sex", "pc1"], model="linear"))
    assert rec["scope"] == "real_fileset" and rec["status"] == "PASSED" and rec["results"]["linear"]["n_unmatched"] == 0


@pytest.mark.plink2
def test_v5_comparator_would_catch_a_wrong_model(tmp_path, plink2_bin):
    """Sensitivity check: omitting the covariates in VariantBridge must make the benchmark FAIL (the test can fail)."""
    mod = load_script()
    wd = tmp_path / "w"; wd.mkdir()
    mod.write_simulated_plink_inputs(wd, seed=7)
    from variantbridge.prep import run_plink2
    run_plink2(plink2_bin, ["--vcf", str(wd / "sim.vcf"), "--make-pgen", "--out", str(wd / "sim")])
    # plink2 with covariates vs VariantBridge without: use pheno col but pass covar to plink2 only
    from variantbridge import association, validation
    import pandas as pd
    lg = []
    args = ["--pfile", str(wd / "sim"), "--pheno", str(wd / "pheno.tsv"), "--pheno-name", "y_lin", "--covar", str(wd / "covar.tsv"),
            "--covar-name", "age", "sex", "pc1", "--glm", "hide-covar", "--out", str(wd / "ref")]
    run_plink2(plink2_bin, args, lg)
    pl = mod.parse_glm(next(wd.glob("ref.y_lin.glm.*")), "linear")
    from variantbridge.io import read_plink_raw
    from variantbridge.prep import read_pvar
    pv = read_pvar(str(wd / "sim") + ".pvar")
    pd.DataFrame({"id": pv["ID"], "allele": pv["ALT"]}).to_csv(wd / "alt.txt", sep="\t", header=False, index=False)
    run_plink2(plink2_bin, ["--pfile", str(wd / "sim"), "--export-allele", str(wd / "alt.txt"), "--export", "A", "--out", str(wd / "raw")], lg)
    iid, dos, vids, counted = read_plink_raw(wd / "raw.raw")
    assert (pv.set_index("ID").loc[vids, "ALT"].to_numpy() == counted).all()   # orientation is right, so only the model differs
    y = pd.read_csv(wd / "pheno.tsv", sep="\t").set_index("#IID").loc[iid, "y_lin"].to_numpy()
    wrong = association.scan(dos.astype(float), y, None, variant_ids=vids, model="linear")
    rep = validation.v5_compare_to_plink2(wrong[["variant", "effect", "se", "p"]], pl[["variant", "effect", "se", "p"]])
    assert rep["status"] == "FAILED"
