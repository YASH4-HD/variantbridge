#!/usr/bin/env python3
"""Offline: GEUVADIS-native cis-eQTL scan (EUR and YRI separately) from scripts/03 outputs.

Outputs (compact, consumed by Streamlit):
    artifacts/geuvadis_eqtl_{EUR,YRI}.csv.gz   one row per gene: best cis variant, effect, SE, p,
                                                permutation p, BH q  (+ manifest)
    <work-dir>/pairs_{pop}.csv.gz              full statistics for published answer-key pairs (for V4)
Optionally runs V3 (permuted-phenotype lambda) and records it.

Model: see variantbridge.eqtl module docstring. No GWAS thresholds, no pooled analysis.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import _common  # noqa: F401
from variantbridge import capability, eqtl, manifest as mf, validation
from variantbridge.io import GenotypeBlock
from variantbridge.prep import parse_answer_key


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", default="data/work/geuvadis")
    ap.add_argument("--artifact-dir", default="artifacts")
    ap.add_argument("--populations", default="EUR,YRI")
    ap.add_argument("--n-perm", type=int, default=1000, help="per-gene permutations for the empirical gene p-value")
    ap.add_argument("--seed", type=int, default=20261005)
    ap.add_argument("--key-best", help="published *.best.* answer key template with {pop}, e.g. data/raw/{pop}.gene.cis.FDR5.best.rs137.txt.gz (EUR373/YRI89 names: pass full path per run)")
    ap.add_argument("--key-files", nargs="*", default=[], help="pairs POP=path of published 'all' and 'best' files to retain pair stats for")
    ap.add_argument("--v3-permutations", type=int, default=0, help="if >0 run V3 on EUR with this many permuted scans")
    ap.add_argument("--v3-max-genes", type=int, default=2000)
    ap.add_argument("--allow-mean-impute", action="store_true")
    ap.add_argument("--max-genes", type=int, default=None, help="restrict to the first N genes (smoke tests only; recorded)")
    a = ap.parse_args()
    wd, ad = Path(a.work_dir), Path(a.artifact_dir)
    ad.mkdir(parents=True, exist_ok=True)
    prep_man = mf.read_manifest(wd / "geuvadis_prep.manifest.json")
    gmeta = pd.read_csv(wd / "gene_meta.csv", index_col=0, dtype={"chrom": str})
    blocks = {p.stem.replace("chr", ""): GenotypeBlock.load(p) for p in sorted((wd / "blocks").glob("chr*.npz"))}
    pops = a.populations.split(",")
    # keep pairs from published files (key=POP=path)
    keep = {p: set() for p in pops}
    for kv in a.key_files:
        pop, path = kv.split("=", 1)
        k = parse_answer_key(path)
        keep[pop] |= set(zip(k["gene_id"], k["snp_id"]))
    pop_blocks, pop_cov, pop_expr = {}, {}, {}
    for pop in pops:
        cov = pd.read_csv(wd / f"covariates_{pop}.csv", index_col=0)
        expr = pd.read_csv(wd / f"expression_{pop}.csv.gz", index_col=0)
        pop_cov[pop], pop_expr[pop] = cov, expr
        pop_blocks[pop] = {c: b.subset_samples(cov.index.to_numpy()) for c, b in blocks.items()}
    # variant mask: MAF > threshold in EITHER population (computed from both populations' samples)
    masks = {}
    for c in blocks:
        masks[c] = eqtl.maf_either_population([pop_blocks[p][c] for p in pops], eqtl.DEFAULT_MAF)
    for pop in pops:
        expr = pop_expr[pop]
        if a.max_genes:
            expr = expr.iloc[:, : a.max_genes]
        best, pairs = eqtl.run_cis_scan(expr, gmeta, pop_blocks[pop], pop_cov[pop], pop, masks, n_perm=a.n_perm,
                                        seed=a.seed, keep_pairs=keep[pop] or None, allow_mean_impute=a.allow_mean_impute)
        out = ad / f"geuvadis_eqtl_{pop}.csv.gz"
        best.to_csv(out, index=False)
        if len(pairs):
            pairs.to_csv(wd / f"pairs_{pop}.csv.gz", index=False)
        man = mf.build_manifest(
            dataset_name=f"GEUVADIS cis-eQTL scan, {pop}", source=prep_man["source"], citation=prep_man["citation"],
            license_access=prep_man["license_access"], access_date=prep_man["access_date"],
            phenotype_definition=prep_man["phenotype_definition"],
            inclusion_exclusion={"n_samples": len(expr), "n_genes_input": int(expr.shape[1]), "n_genes_tested": int(len(best)),
                                 "max_genes_restriction": a.max_genes, "prep_deviations": prep_man["inclusion_exclusion"].get("deviations")},
            qc_thresholds={"maf_gt_either_population": eqtl.DEFAULT_MAF, "cis_window_bp": eqtl.DEFAULT_WINDOW, "autosomes_only": True},
            covariates=list(pop_cov[pop].columns),
            model_specification={"model": "OLS linear, t-test; FWL-vectorised", "multiple_testing": f"direct per-gene permutation (n={a.n_perm}) + BH across genes [project design choice, see METHODS.md]"},
            random_seed=a.seed, limitations=prep_man["limitations"] + ["permutation scheme differs from the 2013 paper's"],
            output_file=str(out), output_sha256=mf.file_sha256(out))
        mf.write_manifest(ad / f"geuvadis_eqtl_{pop}.manifest.json", man)
        print(pop, "genes tested:", len(best), "eGenes (BH q<0.05):", int((best.get("q_bh", pd.Series(dtype=float)) < 0.05).sum()))
    if a.v3_permutations > 0:
        pop = "EUR"
        lam = eqtl.permuted_scan_lambda(pop_expr[pop], gmeta, pop_blocks[pop], pop_cov[pop], a.v3_permutations, a.seed,
                                        variant_mask=masks, max_genes=a.v3_max_genes)
        rec = validation.v3_null_calibration(lam)
        rec["scope"] = f"EUR, first {a.v3_max_genes} genes, {a.v3_permutations} permuted scans"
        capability.record_status("V3", rec, ad)
        print("V3:", rec["status"], rec["median_lambda"])


if __name__ == "__main__":
    main()
