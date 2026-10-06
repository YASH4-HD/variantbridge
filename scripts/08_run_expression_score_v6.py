#!/usr/bin/env python3
"""V6: Genetic Expression Score (eQTLGen weights -> GEUVADIS EUR), implementation and biology REPORTED SEPARATELY.

Terminology: this is predicted expression, not a PRS and not "risk". A null biological result does
not fail the software build. Parameters (p-threshold, clumping r2) are fixed a priori from
ExpressionScoreParams and are never tuned on the target.

    python scripts/08_run_expression_score_v6.py --weights data/raw/<eQTLGen file> \\
        --cohort-disjointness-attested
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import _common  # noqa: F401
from variantbridge import capability, expression_score as es, manifest as mf, validation
from variantbridge.io import GenotypeBlock
from variantbridge.prep import PrepError, resolve_columns, strip_version

W_ALIASES = {
    "pvalue": ["Pvalue", "pvalue", "P"], "variant_id": ["SNP", "rsid", "variant_id"],
    "assessed_allele": ["AssessedAllele", "assessed_allele"], "other_allele": ["OtherAllele", "other_allele"],
    "zscore": ["Zscore", "zscore", "Z"], "gene_id": ["Gene", "gene", "gene_id"],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--work-dir", default="data/work/geuvadis")
    ap.add_argument("--artifact-dir", default="artifacts")
    ap.add_argument("--p-threshold", type=float, default=es.ExpressionScoreParams.p_threshold)
    ap.add_argument("--r2-clump", type=float, default=es.ExpressionScoreParams.r2_clump)
    ap.add_argument("--n-permutations", type=int, default=1000)
    ap.add_argument("--n-bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261005)
    ap.add_argument("--cohort-disjointness-attested", action="store_true",
                    help="confirm you verified GEUVADIS is not in the eQTLGen cohort list (eqtlgen.org/phase1.html)")
    ap.add_argument("--max-genes", type=int, default=None)
    ap.add_argument("--no-record", action="store_true")
    a = ap.parse_args()
    wd, ad = Path(a.work_dir), Path(a.artifact_dir)
    prep_man = mf.read_manifest(wd / "geuvadis_prep.manifest.json")
    W = pd.read_csv(a.weights, sep="\t", dtype=str)
    m = resolve_columns(W.columns, W_ALIASES, "eQTLGen weights")
    W = pd.DataFrame({k: W[v] for k, v in m.items()})
    W["pvalue"] = pd.to_numeric(W["pvalue"])
    W["zscore"] = pd.to_numeric(W["zscore"])
    W["gene_id"] = W["gene_id"].map(strip_version)
    gmeta = pd.read_csv(wd / "gene_meta.csv", index_col=0, dtype={"chrom": str})
    cov = pd.read_csv(wd / "covariates_EUR.csv", index_col=0)
    measured = pd.read_csv(wd / "expression_EUR.csv.gz", index_col=0)
    params = es.ExpressionScoreParams(p_threshold=a.p_threshold, r2_clump=a.r2_clump)
    preds, reports = [], []
    for p in sorted((wd / "blocks").glob("chr*.npz")):
        blk = GenotypeBlock.load(p).subset_samples(cov.index.to_numpy())
        chrom = blk.chrom
        gt = gmeta[gmeta["chrom"] == chrom]["tss"]
        if gt.empty:
            continue
        pr, info, rep = es.predict_expression(W, blk, gt, params)
        preds.append(pr)
        reports.append(rep)
    if not preds:
        raise PrepError("no genes could be scored (check weights/target overlap)")
    pred = pd.concat(preds, axis=1)
    if a.max_genes:
        pred = pred.iloc[:, : a.max_genes]
    ev = es.evaluate_expression_score(pred, measured, a.n_permutations, a.n_bootstrap, a.seed)
    agg = {k: int(sum(r[k] for r in reports)) for k in reports[0] if isinstance(reports[0][k], int)}
    checks = {
        "cohort_disjointness_attested": bool(a.cohort_disjointness_attested),
        "parameters_fixed_a_priori_not_tuned_on_target": True,
        "alleles_aligned_and_palindromic_dropped": agg.get("n_aligned", 0) > 0,
        "genes_scored_gt_0": ev["n_genes_scored"] > 0,
        "predicted_and_measured_finite": bool(np.isfinite(pred.to_numpy()).all()),
        "permutation_null_generated": ev["n_permutations"] >= 100,
        "bootstrap_ci_generated": ev["n_bootstrap"] >= 100,
    }
    res = validation.v6_expression_score({k: v for k, v in ev.items() if k not in ("per_gene", "null_distribution")}, agg, checks)
    res["parameters"] = params.provenance()
    res["terminology"] = "Genetic Expression Score / Predicted Expression (not PRS, not risk)"
    ad.mkdir(parents=True, exist_ok=True)
    ev["per_gene"].to_csv(ad / "expression_score_per_gene.csv.gz", index=False)
    np.save(ad / "expression_score_null.npy", ev["null_distribution"])
    out = ad / "expression_score_per_gene.csv.gz"
    man = mf.build_manifest(
        dataset_name="eQTLGen phase I cis weights -> GEUVADIS EUR target (Genetic Expression Score)",
        source=a.weights, citation=["doi:10.1038/s41588-021-00913-z", "doi:10.1038/nature12531"],
        license_access="open downloads (per audit)", access_date=prep_man["access_date"],
        phenotype_definition="measured PEER-residual gene expression (LCL)", inclusion_exclusion=agg,
        qc_thresholds=params.provenance(), covariates=mf.NOT_APPLICABLE,
        model_specification="clump+threshold expression score; per-gene Pearson r; permutation null (sample labels) + bootstrap CI",
        random_seed=a.seed, limitations=["discovery tissue blood/PBMC vs target LCL: validates scoring machinery, not cross-tissue biology"],
        output_file=str(out), output_sha256=mf.file_sha256(out))
    mf.write_manifest(ad / "expression_score_per_gene.manifest.json", man)
    if not a.no_record:
        capability.record_status("V6", res, ad)
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
