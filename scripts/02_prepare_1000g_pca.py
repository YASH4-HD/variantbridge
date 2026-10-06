#!/usr/bin/env python3
"""Offline: 1000 Genomes -> QC -> LD-pruned PCA -> V1 (+ optional KING relatedness, V2 sex check).

Outputs compact artifacts consumed by Streamlit:
    artifacts/1000g_pca.csv.gz  + artifacts/1000g_pca.manifest.json
    artifacts/1000g_king_pairs.csv.gz (+ manifest) when --king
and records V1 (and V2 if --x-vcf) in artifacts/VALIDATION_STATUS.json.

Methodology (docs/METHODS.md; audit stages 1-4): sample QC (call rate, heterozygosity) ->
variant QC (MAF, missingness; NO pooled-multi-ancestry HWE filter, Wahlund effect) -> LD pruning
(PLINK2 --indep-pairwise, window/step/r2 = project design choices) with long-range-LD exclusion
(region table must be supplied; otherwise an explicit manifest-recorded deviation) -> PCA.
Never silently substitutes: if --long-range-ld-file is absent you MUST pass --no-long-range-exclusion.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import _common  # noqa: F401
from variantbridge import capability, manifest as mf, population_structure as ps, qc, validation
from variantbridge.io import read_plink_raw
from variantbridge.prep import PrepError, find_plink2, plink2_version, read_panel, run_plink2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--panel", required=True, help="1000G sample panel (sample pop super_pop gender)")
    ap.add_argument("--work-dir", default="data/work/1000g")
    ap.add_argument("--artifact-dir", default="artifacts")
    ap.add_argument("--plink2", default=None)
    ap.add_argument("--long-range-ld-file", help="TSV chrom,start,end (1-based inclusive, b37), no header")
    ap.add_argument("--no-long-range-exclusion", action="store_true", help="explicit deviation; recorded in manifest")
    ap.add_argument("--pca-engine", choices=["plink2", "python"], default="plink2")
    ap.add_argument("--n-pcs", type=int, default=10)
    ap.add_argument("--king", action="store_true", help="also compute KING kinship pairs (plink2 --make-king-table)")
    ap.add_argument("--x-vcf", help="chrX VCF for the V2 sex check (optional)")
    ap.add_argument("--source-url", default="NOT RECORDED")
    ap.add_argument("--access-date", default="NOT RECORDED")
    a = ap.parse_args()
    if not a.long_range_ld_file and not a.no_long_range_exclusion:
        sys.exit("Provide --long-range-ld-file (Price 2008 regions) or explicitly pass --no-long-range-exclusion.")
    plink2 = find_plink2(a.plink2)
    wd, ad = Path(a.work_dir), Path(a.artifact_dir)
    wd.mkdir(parents=True, exist_ok=True)
    ad.mkdir(parents=True, exist_ok=True)
    log: list = []
    thr = qc.QCThresholds()
    deviations = []
    base = str(wd / "qc")
    run_plink2(plink2, ["--vcf", a.vcf, "--max-alleles", "2", "--snps-only", "just-acgt",
                        "--set-missing-var-ids", "@:#:$r:$a", "--rm-dup", "force-first",
                        "--maf", str(thr.maf_min), "--geno", str(thr.variant_missing_max),
                        "--make-pgen", "--out", base], log)
    # ---- sample QC
    run_plink2(plink2, ["--pfile", base, "--missing", "--het", "--out", str(wd / "sq")], log)
    smiss = pd.read_csv(wd / "sq.smiss", sep="\t")
    het = pd.read_csv(wd / "sq.het", sep="\t")
    iid = "#IID" if "#IID" in smiss else smiss.columns[0]
    smet = pd.DataFrame({"call_rate": 1 - smiss["F_MISS"].to_numpy()}, index=smiss[iid].astype(str))
    ho = het.set_index(het.columns[0] if "#IID" not in het else "#IID")
    smet["het_rate"] = (1 - ho["O(HOM)"] / ho["OBS_CT"]).reindex(smet.index).to_numpy()
    ok = qc.sample_qc_mask(smet, thr)
    (wd / "remove.txt").write_text("\n".join(["#IID"] + list(smet.index[~ok])) + "\n")
    # ---- LD pruning
    prune_args = ["--pfile", base, "--remove", str(wd / "remove.txt"), "--indep-pairwise", "50", "5", "0.2", "--out", str(wd / "prune")]
    if a.long_range_ld_file:
        lr = pd.read_csv(a.long_range_ld_file, sep="\t", header=None, names=["chrom", "start", "end"])
        lr.to_csv(wd / "longrange.bed1", sep="\t", header=False, index=False)
        prune_args = ["--pfile", base, "--remove", str(wd / "remove.txt"), "--exclude", "bed1", str(wd / "longrange.bed1"),
                      "--indep-pairwise", "50", "5", "0.2", "--out", str(wd / "prune")]
    else:
        deviations.append("long-range LD regions NOT excluded before PCA (explicit --no-long-range-exclusion)")
    run_plink2(plink2, prune_args, log)
    # ---- PCA
    if a.pca_engine == "plink2":
        run_plink2(plink2, ["--pfile", base, "--remove", str(wd / "remove.txt"), "--extract", str(wd / "prune.prune.in"),
                            "--pca", str(a.n_pcs), "--out", str(wd / "pca")], log)
        ev = pd.read_csv(wd / "pca.eigenvec", sep="\t")
        ev.columns = [c.lstrip("#") for c in ev.columns]
        pcs = ev.rename(columns={"IID": "sample_id"}).drop(columns=[c for c in ("FID",) if c in ev])
        pcs["sample_id"] = pcs["sample_id"].astype(str)
        engine = "plink2 --pca"
    else:
        run_plink2(plink2, ["--pfile", base, "--remove", str(wd / "remove.txt"), "--extract", str(wd / "prune.prune.in"),
                            "--export", "A", "--out", str(wd / "pr")], log)
        s, dos, _, _ = read_plink_raw(wd / "pr.raw")
        p, _ratio, _ = ps.genotype_pca(dos.astype(float), a.n_pcs)
        pcs = ps.pcs_frame(p, s)
        engine = "variantbridge.population_structure.genotype_pca"
    panel = read_panel(a.panel).rename(columns={"sample": "sample_id"})
    out = pcs.merge(panel, on="sample_id", how="left")
    if out["super_pop"].isna().any():
        deviations.append(f"{int(out['super_pop'].isna().sum())} samples absent from panel (kept without label)")
    outp = ad / "1000g_pca.csv.gz"
    out.to_csv(outp, index=False)
    lab = out.dropna(subset=["super_pop"]).set_index("sample_id")["super_pop"]
    v1 = validation.v1_population_structure(out[out["sample_id"].isin(lab.index)][["sample_id"] + [c for c in out if c.startswith("PC")]], lab)
    v1["engine"] = engine
    capability.record_status("V1", v1, ad)
    n_vars = sum(1 for _ in open(wd / "prune.prune.in"))
    man = mf.build_manifest(
        dataset_name="1000 Genomes Project (phase 3 / IGSR), GRCh37", source=a.source_url,
        citation=["doi:10.1038/nature15393", "doi:10.1093/nar/gkz836"],
        license_access="Open (Fort Lauderdale data-use framework); no registration", access_date=a.access_date,
        phenotype_definition=mf.NOT_APPLICABLE + " (genotype-only; superpopulation labels from the sample panel)",
        inclusion_exclusion={"samples_in_panel_labelled": int(lab.shape[0]), "samples_removed_qc": int((~ok).sum()),
                             "variants_after_prune": n_vars, "deviations": deviations},
        qc_thresholds={**thr.provenance(), "hwe": "NOT APPLIED (pooled multi-ancestry; Wahlund)",
                       "ld_prune": "plink2 --indep-pairwise 50 5 0.2 [project design choice]"},
        covariates=mf.NOT_APPLICABLE, model_specification={"pca_engine": engine, "n_pcs": a.n_pcs,
                                                          "long_range_ld_excluded": bool(a.long_range_ld_file)},
        random_seed=0, limitations=["genotype-only resource; no disease phenotypes", "single chromosome unless full genome supplied"],
        software_versions=mf.software_versions({"plink2": plink2_version(plink2)}),
        output_file=str(outp), output_sha256=mf.file_sha256(outp), plink2_commands=log)
    mf.write_manifest(ad / "1000g_pca.manifest.json", man)
    print(json.dumps({"V1": v1["status"], "balanced_accuracy": v1["knn_cv_balanced_accuracy"]}, indent=2))
    if a.king:
        run_plink2(plink2, ["--pfile", base, "--remove", str(wd / "remove.txt"), "--extract", str(wd / "prune.prune.in"),
                            "--make-king-table", "--king-table-filter", "0.0442", "--out", str(wd / "king")], log)
        k = pd.read_csv(wd / "king.kin0", sep="\t")
        k.columns = [c.lstrip("#") for c in k.columns]
        k["degree"] = [qc.classify_kinship(x) for x in k["KINSHIP"]]
        kp = ad / "1000g_king_pairs.csv.gz"
        k.to_csv(kp, index=False)
        man2 = dict(man, output_file=str(kp), output_sha256=mf.file_sha256(kp),
                    model_specification={"method": "KING-robust (plink2 --make-king-table)", "filter": 0.0442})
        mf.write_manifest(ad / "1000g_king_pairs.manifest.json", man2)
    if a.x_vcf:
        xb = str(wd / "x")
        run_plink2(plink2, ["--vcf", a.x_vcf, "--split-par", "b37", "--max-alleles", "2", "--make-pgen", "--out", xb], log)
        sx = panel[["sample_id", "gender"]].assign(SEX=lambda d: d["gender"].str.lower().map({"male": 1, "female": 2}))
        sx[["sample_id", "SEX"]].rename(columns={"sample_id": "#IID"}).to_csv(wd / "sex.txt", sep="\t", index=False)
        run_plink2(plink2, ["--pfile", xb, "--update-sex", str(wd / "sex.txt"), "--check-sex", "--out", str(wd / "sexchk")], log)
        sc = pd.read_csv(wd / "sexchk.sexcheck", sep="\t")
        sc.columns = [c.lstrip("#") for c in sc.columns]
        conc = qc.sex_concordance(sc["SNPSEX"].to_numpy(), sc["PEDSEX"].to_numpy())
        capability.record_status("V2", validation.v2_sex_check(conc), ad)
        print("V2 recorded", conc)


if __name__ == "__main__":
    main()
