#!/usr/bin/env python3
"""V5 (RELEASE-BLOCKING): benchmark VariantBridge association statistics against `plink2 --glm`.

Modes
  --simulate (default): deterministic SIMULATED dataset written to PLINK format. This validates the
      ARITHMETIC EQUIVALENCE of the engines on an identical dataset/model; it says nothing about
      biology and is labelled SIMULATED in the record.
  --pfile PREFIX --pheno-file F --pheno-col C [--covar-file F --covar-cols a,b] --model linear|logistic:
      run on a real (small, e.g. chr22 subset) PLINK2 fileset (pheno/covar TSVs keyed by '#IID';
      binary phenotypes coded 1/2). Recorded with scope=real_fileset.

Writes artifacts/VALIDATION_STATUS.json (key V5) only after the comparison actually ran.
Status is PASSED/FAILED against tolerances pre-registered in variantbridge.validation.
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import _common  # noqa: F401
from variantbridge import association, capability, validation
from variantbridge.io import read_plink_raw
from variantbridge.prep import PrepError, find_plink2, plink2_version, read_pvar, run_plink2


def write_simulated_plink_inputs(outdir: Path, seed: int = 20261005, n: int = 300, m: int = 400):
    r = np.random.default_rng(seed)
    maf = r.uniform(0.05, 0.5, m)
    G = r.binomial(2, maf, (n, m)).astype(float)
    G[r.random(G.shape) < 0.01] = np.nan
    age = r.normal(50, 12, n)
    sex = r.integers(0, 2, n).astype(float)
    pc1 = r.normal(0, 1, n)
    Gf = np.nan_to_num(G, nan=0.0)
    beta = np.zeros(m)
    beta[r.choice(m, 20, replace=False)] = r.normal(0, 0.35, 20)
    lin = Gf @ beta + 0.03 * (age - 50) + 0.4 * sex + 0.5 * pc1 + r.normal(0, 1.0, n)
    p = 1 / (1 + np.exp(-(-0.5 + 0.5 * (Gf @ beta) + 0.02 * (age - 50) + 0.3 * pc1)))
    binv = r.binomial(1, p)
    sids = [f"S{i:04d}" for i in range(n)]
    vids = [f"snp{j:04d}" for j in range(m)]
    lines = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sids)]
    code = {0.0: "0/0", 1.0: "0/1", 2.0: "1/1"}
    for j in range(m):
        gts = ["./." if np.isnan(g) else code[g] for g in G[:, j]]
        lines.append("\t".join(["22", str(16_000_000 + 1000 * j), vids[j], "A", "G", ".", ".", ".", "GT"] + gts))
    (outdir / "sim.vcf").write_text("\n".join(lines) + "\n")
    pd.DataFrame({"#IID": sids, "y_lin": lin, "y_bin": binv + 1}).to_csv(outdir / "pheno.tsv", sep="\t", index=False)
    pd.DataFrame({"#IID": sids, "age": age, "sex": sex, "pc1": pc1}).to_csv(outdir / "covar.tsv", sep="\t", index=False)
    return G, np.column_stack([age, sex, pc1]), lin, binv.astype(float), vids, sids


def parse_glm(path: Path, model: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[df["TEST"] == "ADD"]
    if "ERRCODE" in df and df["ERRCODE"].astype(str).str.strip().replace({".": ""}).ne("").any():
        bad = df[df["ERRCODE"].astype(str).str.strip().replace({".": ""}).ne("")]
        print(f"WARNING: {len(bad)} variants with PLINK2 ERRCODE excluded from comparison", file=sys.stderr)
        df = df.drop(bad.index)
    if model == "linear":
        out = pd.DataFrame({"variant": df["ID"], "effect": df["BETA"], "se": df["SE"], "p": df["P"], "n": df["OBS_CT"],
                            "a1": df["A1"], "alt": df["ALT"], "ref": df["REF"]})
    else:
        out = pd.DataFrame({"variant": df["ID"], "effect": df["OR"], "se": df["LOG(OR)_SE"], "p": df["P"],
                            "n": df["OBS_CT"], "a1": df["A1"], "alt": df["ALT"], "ref": df["REF"]})
    # PLINK2 tests A1, which is not always ALT (observed: A1 = REF for some variants). VariantBridge
    # models the ALT dosage, so re-express PLINK2 results on the ALT allele: beta -> -beta, OR -> 1/OR
    # (SE of beta/log-OR and P are unchanged). The number flipped is recorded in the report.
    neither = ~((out["a1"] == out["alt"]) | (out["a1"] == out["ref"]))
    if neither.any():
        raise PrepError("PLINK2 tested allele is neither REF nor ALT for some variants (multi-allelic input?)")
    flip = (out["a1"] == out["ref"]) & (out["a1"] != out["alt"])
    out.loc[flip, "effect"] = (-out.loc[flip, "effect"]) if model == "linear" else (1.0 / out.loc[flip, "effect"])
    out.attrs["n_a1_is_ref_reoriented"] = int(flip.sum())
    return out.drop(columns=["a1", "alt", "ref"])


def benchmark_pfile(workdir: Path, plink2: str, pfile: str, pheno_file: Path, pheno_col: str, covar_file: Path | None,
                    covar_cols: list[str], model: str, log: list) -> dict:
    """Run `plink2 --glm` and VariantBridge on the SAME PLINK2 fileset, phenotype and covariates; compare.

    VariantBridge genotypes are exported from the same fileset with the ALT allele counted explicitly.
    Binary phenotypes in the pheno file must use PLINK coding (1 = control, 2 = case).
    """
    tag = f"{model}_{pheno_col}"
    args = ["--pfile", pfile, "--pheno", str(pheno_file), "--pheno-name", pheno_col]
    if covar_file is not None:
        args += ["--covar", str(covar_file), "--covar-name", *covar_cols]
    args += ["--glm", "hide-covar", *(["no-firth"] if model == "logistic" else []), "--out", str(workdir / f"res_{tag}")]
    run_plink2(plink2, args, log)
    pl = parse_glm(next(workdir.glob(f"res_{tag}.{pheno_col}.glm.*")), model)
    # export the SAME genotypes with the ALT allele counted
    pv = read_pvar(f"{pfile}.pvar")
    pd.DataFrame({"id": pv["ID"], "allele": pv["ALT"]}).to_csv(workdir / f"alt_{tag}.txt", sep="\t", header=False, index=False)
    run_plink2(plink2, ["--pfile", pfile, "--export-allele", str(workdir / f"alt_{tag}.txt"), "--export", "A",
                        "--out", str(workdir / f"raw_{tag}")], log)
    iid, dos, vids, counted = read_plink_raw(workdir / f"raw_{tag}.raw")
    alt = pv.set_index("ID").loc[vids, "ALT"].to_numpy()
    if not (np.asarray(counted) == alt).all():
        raise PrepError("exported counted allele != ALT")
    ph = pd.read_csv(pheno_file, sep="\t").set_index("#IID").loc[iid]
    y = ph[pheno_col].to_numpy(dtype=float)
    if model == "logistic":
        y = y - 1.0
    cov = None
    if covar_file is not None:
        cov = pd.read_csv(covar_file, sep="\t").set_index("#IID").loc[iid, covar_cols].to_numpy(dtype=float)
    mine = association.scan(dos.astype(float), y, cov, variant_ids=vids, model=model)
    problems = association.check_statistical_consistency(mine)
    m = mine[["variant", "effect", "se", "p", "n", "odds_ratio"]].copy()
    if model == "logistic":
        m["effect"] = m["odds_ratio"]
    rep = validation.v5_compare_to_plink2(m[["variant", "effect", "se", "p"]], pl[["variant", "effect", "se", "p"]], model=model)
    pl_n = pl.set_index("variant")["n"]
    n_ok = bool((m.set_index("variant")["n"].reindex(pl_n.index) == pl_n).all())
    rep["obs_count_identical"] = n_ok
    rep["n_plink2_a1_is_ref_reoriented_to_alt"] = int(pl.attrs.get("n_a1_is_ref_reoriented", 0))
    rep["internal_consistency_violations"] = len(problems)
    rep["effect_compared_on"] = "BETA" if model == "linear" else "OR (se compared on log-odds scale)"
    if not n_ok or problems:
        rep["status"] = capability.FAILED
    return rep


def run_benchmark(workdir: Path, plink2: str, simulate: bool = True, record_root: str | None = "artifacts",
                  seed: int = 20261005, real: dict | None = None) -> dict:
    """Simulated mode (default) or real-fileset mode (``real`` = dict with pfile/pheno_file/pheno_col/covar_file/covar_cols/model)."""
    workdir.mkdir(parents=True, exist_ok=True)
    log: list = []
    results = {}
    if real is None:
        write_simulated_plink_inputs(workdir, seed)
        base = str(workdir / "sim")
        run_plink2(plink2, ["--vcf", str(workdir / "sim.vcf"), "--make-pgen", "--out", base], log)
        for model, col in (("linear", "y_lin"), ("logistic", "y_bin")):
            results[model] = benchmark_pfile(workdir, plink2, base, workdir / "pheno.tsv", col, workdir / "covar.tsv",
                                             ["age", "sex", "pc1"], model, log)
        dataset = ("SIMULATED (deterministic, seed=%d): software arithmetic-equivalence benchmark only; not biological validation" % seed)
        scope = "simulated"
    else:
        results[real["model"]] = benchmark_pfile(workdir, plink2, real["pfile"], Path(real["pheno_file"]), real["pheno_col"],
                                                 Path(real["covar_file"]) if real.get("covar_file") else None,
                                                 real.get("covar_cols", []), real["model"], log)
        dataset = f"real PLINK2 fileset supplied by user: {real['pfile']}"
        scope = "real_fileset"
    overall_ok = all(r["status"] == capability.PASSED for r in results.values())
    record = {
        "status": capability.PASSED if overall_ok else capability.FAILED,
        "executed_at": validation.now(), "scope": scope, "dataset": dataset,
        "plink2_version": plink2_version(plink2), "results": results,
        "logistic_note": "plink2 --glm ... no-firth (Wald); VariantBridge IRLS Wald. No Firth implemented.",
        "plink2_commands": log,
    }
    if scope == "simulated":
        record["remaining"] = "Repeat on a real chr22 subset with --pfile (scope=real_fileset) before claiming real-data benchmark."
    if record_root:
        capability.record_status("V5", record, record_root)
    return record


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plink2", default=None)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--artifact-dir", default="artifacts")
    ap.add_argument("--no-record", action="store_true", help="do not write VALIDATION_STATUS.json")
    ap.add_argument("--pfile", help="REAL mode: PLINK2 fileset prefix (small subset, e.g. chr22 slice)")
    ap.add_argument("--pheno-file"); ap.add_argument("--pheno-col")
    ap.add_argument("--covar-file"); ap.add_argument("--covar-cols", default="")
    ap.add_argument("--model", choices=["linear", "logistic"], default="linear")
    a = ap.parse_args()
    plink2 = find_plink2(a.plink2)
    wd = Path(a.workdir) if a.workdir else Path(tempfile.mkdtemp(prefix="vb_v5_"))
    real = None
    if a.pfile:
        if not (a.pheno_file and a.pheno_col):
            sys.exit("--pfile mode needs --pheno-file and --pheno-col")
        real = dict(pfile=a.pfile, pheno_file=a.pheno_file, pheno_col=a.pheno_col, covar_file=a.covar_file,
                    covar_cols=[c for c in a.covar_cols.split(",") if c], model=a.model)
    rec = run_benchmark(wd, plink2, record_root=None if a.no_record else a.artifact_dir, real=real)
    print(json.dumps({k: v for k, v in rec.items() if k != "plink2_commands"}, indent=2, default=str))
    sys.exit(0 if rec["status"] == capability.PASSED else 1)


if __name__ == "__main__":
    main()
