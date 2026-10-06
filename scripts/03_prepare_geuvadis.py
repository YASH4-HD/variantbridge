#!/usr/bin/env python3
"""Offline: GEUVADIS genotypes + expression -> compact per-chromosome genotype blocks, expression
matrices per population (EUR/YRI), and covariates (population-specific PCs + imputation status).

Everything is written under --work-dir (git-ignored). No raw data enter Git; Streamlit never
reads these files. Methodology follows audit Stage 0-5 (source: Lappalainen 2013 SI via audit):
 * EUR (CEU+FIN+GBR+TSI) and YRI analysed separately; variants kept when MAF > 5% in EITHER population;
 * covariates: PCs 1-3 (EUR) / 1-2 (YRI) computed within population from LD-pruned genotypes
   [PCA engine: variantbridge.population_structure.genotype_pca] + imputation status (0|1);
 * expression: the supplied PEER-residual matrix is used as-is (resk10); the ">0 in >90% of
   individuals" filter is applied only if --rpkm-file is supplied, else recorded as a deviation;
 * dbSNP137 ID conversion applied when --id-converter is supplied (column order UNVERIFIED, see --converter-cols).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import _common  # noqa: F401
from variantbridge import eqtl, manifest as mf, population_structure as ps
from variantbridge.io import GenotypeBlock, read_plink_raw
from variantbridge.prep import (PrepError, find_plink2, geuvadis_population, parse_expression_matrix, plink2_version,
                                read_id_converter, read_panel, read_pvar, run_plink2)


def parse_chroms(s: str) -> list[str]:
    out = []
    for part in s.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            out += [str(i) for i in range(int(lo), int(hi) + 1)]
        else:
            out.append(part)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf-template", required=True, help="path with {chrom}, e.g. data/raw/GEUVADIS.chr{chrom}.PH1PH2_465....vcf.gz")
    ap.add_argument("--chroms", default="22", help="e.g. 22 or 1-22 (autosomes only)")
    ap.add_argument("--panel", required=True)
    ap.add_argument("--expression", required=True)
    ap.add_argument("--expr-gene-col"); ap.add_argument("--expr-chrom-col"); ap.add_argument("--expr-tss-col")
    ap.add_argument("--imputation-status-file", help="TSV: sample<TAB>imputed(0/1)")
    ap.add_argument("--no-imputation-covariate", action="store_true", help="explicit deviation; recorded")
    ap.add_argument("--rpkm-file", help="non-PEER RPKM matrix for the >90%%-expressed filter (optional)")
    ap.add_argument("--id-converter"); ap.add_argument("--converter-cols", default="0,1")
    ap.add_argument("--long-range-ld-file"); ap.add_argument("--no-long-range-exclusion", action="store_true")
    ap.add_argument("--maf", type=float, default=eqtl.DEFAULT_MAF)
    ap.add_argument("--work-dir", default="data/work/geuvadis")
    ap.add_argument("--plink2", default=None)
    ap.add_argument("--source-url", default="NOT RECORDED"); ap.add_argument("--access-date", default="NOT RECORDED")
    a = ap.parse_args()
    if not a.imputation_status_file and not a.no_imputation_covariate:
        sys.exit("Provide --imputation-status-file or explicitly pass --no-imputation-covariate (recorded deviation).")
    if not a.long_range_ld_file and not a.no_long_range_exclusion:
        sys.exit("Provide --long-range-ld-file or explicitly pass --no-long-range-exclusion (recorded deviation).")
    plink2 = find_plink2(a.plink2)
    wd = Path(a.work_dir)
    (wd / "blocks").mkdir(parents=True, exist_ok=True)
    log, deviations = [], []
    panel = read_panel(a.panel)
    panel["group"] = panel["pop"].map(geuvadis_population)
    panel = panel.dropna(subset=["group"])
    expr_all, gmeta = parse_expression_matrix(a.expression, None, a.expr_gene_col, a.expr_chrom_col, a.expr_tss_col)
    samples = [s for s in panel["sample"] if s in expr_all.index]
    group = panel.set_index("sample").loc[samples, "group"]
    if a.rpkm_file:
        rp, _ = parse_expression_matrix(a.rpkm_file, None, a.expr_gene_col, a.expr_chrom_col, a.expr_tss_col)
        expressed = {}
        for g in ("EUR", "YRI"):
            sub = rp.loc[[s for s in samples if group[s] == g]]
            expressed[g] = set(sub.columns[(sub > 0).mean() > 0.9])
    else:
        expressed = None
        deviations.append("expressed (>0) in >90% of individuals filter NOT applied (no --rpkm-file supplied)")
    if a.no_imputation_covariate:
        deviations.append("imputation-status covariate NOT used (explicit deviation from published model)")
    if a.no_long_range_exclusion:
        deviations.append("long-range LD regions NOT excluded before PCA (explicit deviation)")
    for g in ("EUR", "YRI"):
        (wd / f"keep_{g}.txt").write_text("\n".join(["#IID"] + [s for s in samples if group[s] == g]) + "\n")
    conv = None
    if a.id_converter:
        c0, c1 = (int(x) for x in a.converter_cols.split(","))
        conv = read_id_converter(a.id_converter, c0, c1)
        conv.to_csv(wd / "idconvert.txt", sep="\t", header=False, index=False)
    lr_args = []
    if a.long_range_ld_file:
        lr = pd.read_csv(a.long_range_ld_file, sep="\t", header=None, names=["chrom", "start", "end"])
        lr.to_csv(wd / "longrange.bed1", sep="\t", header=False, index=False)
        lr_args = ["--exclude", "bed1", str(wd / "longrange.bed1")]
    pca_cols = {g: [] for g in ("EUR", "YRI")}
    pca_ids = {g: [] for g in ("EUR", "YRI")}
    chroms = parse_chroms(a.chroms)
    for ch in chroms:
        base = str(wd / f"chr{ch}")
        args = ["--vcf", a.vcf_template.format(chrom=ch), "--max-alleles", "2", "--snps-only", "just-acgt"]
        if conv is not None:
            args += ["--update-name", str(wd / "idconvert.txt")]
        run_plink2(plink2, args + ["--make-pgen", "--out", base], log)
        for g in ("EUR", "YRI"):
            run_plink2(plink2, ["--pfile", base, "--keep", str(wd / f"keep_{g}.txt"), "--maf", str(a.maf),
                                "--write-snplist", "--out", str(wd / f"snps_{g}_chr{ch}")], log)
        union = sorted(set(open(wd / f"snps_EUR_chr{ch}.snplist").read().split()) |
                       set(open(wd / f"snps_YRI_chr{ch}.snplist").read().split()))
        (wd / f"union_chr{ch}.txt").write_text("\n".join(union) + "\n")
        keep_all = wd / "keep_all.txt"
        keep_all.write_text("\n".join(["#IID"] + samples) + "\n")
        # plink2 --export A counts REF by default (observed); count the ALT allele EXPLICITLY so effect signs
        # have a fixed, documented orientation (effect = change per ALT allele).
        pv0 = read_pvar(wd / f"chr{ch}.pvar").drop_duplicates("ID").set_index("ID")
        pd.DataFrame({"id": union, "allele": pv0.loc[union, "ALT"].to_numpy()}).to_csv(
            wd / f"alt_chr{ch}.txt", sep="\t", header=False, index=False)
        run_plink2(plink2, ["--pfile", base, "--keep", str(keep_all), "--extract", str(wd / f"union_chr{ch}.txt"),
                            "--export-allele", str(wd / f"alt_chr{ch}.txt"),
                            "--export", "A", "--out", str(wd / f"raw_chr{ch}")], log)
        sid, dos, vids, counted = read_plink_raw(wd / f"raw_chr{ch}.raw")
        pv = read_pvar(wd / f"chr{ch}.pvar")
        pv = pv.drop_duplicates("ID").set_index("ID")
        pos = pv.loc[vids, "POS"].astype(np.int64).to_numpy()
        ref, alt = pv.loc[vids, "REF"].to_numpy(), pv.loc[vids, "ALT"].to_numpy()
        counted = np.array(counted)
        other = np.where(counted == alt, ref, alt)
        bad = counted != alt
        if bad.any():
            raise PrepError(f"counted allele != ALT for {int(bad.sum())} variants on chr{ch}: --export-allele not honoured")
        order = np.argsort(pos, kind="stable")
        blk = GenotypeBlock(ch, pos[order], np.array(vids)[order], counted[order], other[order], dos[:, order],
                            np.array(sid))
        blk = blk.subset_samples(samples)
        blk.save(wd / "blocks" / f"chr{ch}.npz")
        # LD-pruned variants per population for PCA (plink2, within population)
        for g in ("EUR", "YRI"):
            run_plink2(plink2, ["--pfile", base, "--keep", str(wd / f"keep_{g}.txt"), "--extract", str(wd / f"snps_{g}_chr{ch}.snplist"),
                                *lr_args, "--indep-pairwise", "50", "5", "0.2", "--out", str(wd / f"prune_{g}_chr{ch}")], log)
            keep_ids = set(open(wd / f"prune_{g}_chr{ch}.prune.in").read().split())
            m = np.array([v in keep_ids for v in blk.variant_ids])
            gi = np.flatnonzero(group.loc[blk.sample_ids].to_numpy() == g)
            pca_cols[g].append(blk.dosage[np.ix_(gi, np.flatnonzero(m))])
        print(f"chr{ch}: {len(blk.variant_ids)} variants, {len(samples)} samples")
    # PCs per population
    cov_info = {}
    for g in ("EUR", "YRI"):
        gs = [s for s in samples if group[s] == g]
        D = np.hstack(pca_cols[g]).astype(float)
        pcs, ratio, used = ps.genotype_pca(D, 10)
        pcdf = ps.pcs_frame(pcs, gs)
        pcdf.to_csv(wd / f"pcs_{g}.csv", index=False)
        imp = None
        if a.imputation_status_file:
            imp = pd.read_csv(a.imputation_status_file, sep=None, engine="python", header=None, names=["sample", "imputed"])
            imp = imp.set_index("sample")["imputed"]
            missing_imp = [s for s in gs if s not in imp.index]
            if missing_imp:
                raise PrepError(f"imputation status missing for {len(missing_imp)} samples, e.g. {missing_imp[:3]}")
        cov = eqtl.build_covariates(gs, pcdf, eqtl.PCS_BY_POPULATION[g], imp)
        cov.to_csv(wd / f"covariates_{g}.csv")
        ex = expr_all.loc[gs]
        if expressed is not None:
            ex = ex[[c for c in ex.columns if c in expressed[g]]]
        ex.to_csv(wd / f"expression_{g}.csv.gz")
        cov_info[g] = {"n_samples": len(gs), "n_pca_variants": int(used.sum()), "n_genes": int(ex.shape[1]),
                       "variance_ratio_pc1_3": [float(x) for x in ratio[:3]]}
    gmeta.to_csv(wd / "gene_meta.csv")
    man = mf.build_manifest(
        dataset_name="GEUVADIS (E-GEUV-1) mRNA-seq + 1000G phase1 genotypes, GRCh37", source=a.source_url,
        citation=["doi:10.1038/nature12531"], license_access="freely and openly available with no restrictions (per audit)",
        access_date=a.access_date, phenotype_definition="PEER-residual (K=10), standard-normal-transformed gene expression, LCL",
        inclusion_exclusion={"chromosomes": chroms, "populations": cov_info, "deviations": deviations,
                             "snps_only_biallelic": True, "maf_either_population": a.maf},
        qc_thresholds={"maf_gt_in_either_population": a.maf, "hwe": "NOT APPLIED (published pipeline; imputed/phased)"},
        covariates={"EUR": "PC1-3 (+imputation status)", "YRI": "PC1-2 (+imputation status)"},
        model_specification="data preparation only; association model in scripts/05_run_eqtl_scan.py", random_seed=0,
        limitations=["LCL tissue", "EUR+YRI only", "Gencode v12 / 2013 quantification", "format assumptions UNVERIFIED: see docs/DATASETS.md"],
        software_versions=mf.software_versions({"plink2": plink2_version(plink2)}), plink2_commands=log)
    mf.write_manifest(wd / "geuvadis_prep.manifest.json", man)
    print("prepared; deviations:", deviations)


if __name__ == "__main__":
    main()
