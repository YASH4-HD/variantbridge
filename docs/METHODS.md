# METHODS

Every method is labelled **[source]** (taken from the CohortOmics dataset audit, which cites the
primary literature) or **[design choice]** (decided in this project and *not* claimed to come from
a paper). Nothing here has been run on real data inside this repository's delivery; see
`TEST_RESULTS.md` and `artifacts/VALIDATION_STATUS.json` (absent ⇒ NOT RUN).

## Architecture
`src/variantbridge/` holds all science with no Streamlit import. `app.py` is a thin viewer.
`scripts/` are offline, resumable and produce compact artifacts with manifests. The app never
downloads or processes large data (amendment 4).

## Statistical consistency (release-blocking, amendment 5) — `association.py`
Effect, SE, CI and p all come from **one fitted model** per variant.
- Quantitative phenotype: OLS with a t-test, df = n − k − 2 (k covariates), vectorised by
  Frisch–Waugh–Lovell residualisation. [source: standard theory]
- Binary phenotype: IRLS logistic regression with a Wald z-test; no Firth correction.
  Non-convergence, separation and degenerate columns are flagged in `status`, never silently
  reported. [design choice: no Firth, to match PLINK2 `no-firth`]
- v0.1 reported a logistic coefficient with a point-biserial p-value. That is removed.
- `check_statistical_consistency` recomputes stat/p/CI from (effect, SE, df) and fails on mismatch.
- **V5**: same subset, same model, compared with `plink2 --glm hide-covar no-firth`. A1=REF rows
  are re-oriented to ALT (β→−β, OR→1/OR) and counted in the report. Pre-registered tolerances:
  relative 1e-4 on β/SE, absolute 5e-3 on −log10 p (p > 1e-30). [design choice]

## QC — `qc.py`
Sample call rate/heterozygosity, variant call rate/MAF, exact HWE (Wigginton 2005), X-heterozygosity
sex inference (V2), KING-robust kinship. Thresholds carry provenance labels in `QCThresholds`.
Long-range LD regions are **not hard-coded** (the audit flags region coordinates as unverified):
scripts require `--long-range-ld-file` or an explicit `--no-long-range-exclusion`, which is
recorded in the manifest as a deviation.

## Population structure — `population_structure.py`
Patterson-normalised genotype PCA on LD-pruned autosomal variants (Gram-matrix eigendecomposition,
deterministic sign convention). **V1 (amendment 2)** is "leading PCs recover broad 1000G
superpopulation structure", quantified against labels: kNN 5-fold CV balanced accuracy on PCs 1–5
(≥ 0.90), per-PC η², silhouette. Thresholds pre-registered. [design choice] The old claim that
"PC1 separates superpopulations" is not used.
`feature_pca_v01` preserves the v0.1 scaler+PCA for the EvoResist-AI branch only.

## GEUVADIS-native cis-eQTL — `eqtl.py`
- Populations analysed **separately** (EUR, YRI); no pooled analysis. [source: Lappalainen 2013 SI via audit]
- Covariates: imputation status + PCs 1–3 (EUR) / 1–2 (YRI). [source]
- ±1 Mb of the gene coordinate, MAF > 5 % in either population, autosomes. [source]
- PEER (resk10) normalised expression as provided by GEUVADIS. [source]
- Significance: direct per-gene permutation (best |r| over cis variants vs. permuted
  residualised phenotype), then Benjamini–Hochberg across genes. [design choice — differs from the
  2013 paper's scheme; replication therefore tests agreement, not identity of procedure]
- No 5×10⁻⁸ threshold; no LD pruning of tested variants. [source]
- The permutation p-value has a resolution floor 1/(n_perm+1) (empirical p uses the +1 correction).
- **V3**: permuted-phenotype λ (|median λ − 1| ≤ 0.10). **V4 (release-blocking, amendment 6)**:
  EUR and YRI each compared with the published answer key: same-direction fraction ≥ 0.90 over
  ≥ 100 matched pairs. Allele orientation is never auto-flipped; a near-zero agreement rate is
  flagged as a probable orientation error rather than corrected. Replication is not claimed until
  `VALIDATION_STATUS.json` records an executed PASS for both populations.

## Genetic Expression Score — `expression_score.py` (amendment 1)
Named **Genetic Expression Score / Predicted Expression**; never PRS, never risk.
eQTLGen (blood cis-eQTL z-scores) weights → clump + p-threshold fixed a priori (p ≤ 1e-4,
r² ≤ 0.1) [design choice] → palindromic SNPs dropped, alleles aligned → predicted expression in
GEUVADIS. Primary statistic: mean per-gene Pearson r; bootstrap CI; **sample-label permutation null
(permuted jointly across genes)** and empirical p; per-gene Fisher-z CIs.
**V6 (amendment 3)** reports observed performance, uncertainty, null and empirical p.
*Implementation success* (cohorts disjoint, alleles aligned, genes covered) is recorded separately
from *biological evidence*. A null biological result is a valid outcome and does not fail the build.
GEUVADIS LCLs vs eQTLGen blood is a tissue mismatch; eQTLGen includes GEUVADIS-derived samples
unless excluded — `--cohort-disjointness-attested` forces the operator to state this explicitly.

## PRS — `prs.py`
Reserved for genuine complex traits, with `assert_disjoint_cohorts` leakage guard. With open data
there is no valid disease-trait PRS (needs controlled-access target); the UI shows only a
SIMULATED PRS demo on disjoint simulated cohorts.

## Provenance — `manifest.py`, `artifacts.py`, `capability.py`
Manifests record the handoff §12 / audit E fields plus file SHA-256. The artifact loader refuses
missing manifests and checksum mismatches. Capability tiers: TECHNICALLY POSSIBLE, STATISTICALLY
DEFENSIBLE, PORTFOLIO DEMONSTRATION ONLY, SIMULATED. A tier is only a *target*; the validation
status (NOT RUN / PASSED / FAILED) is displayed next to it.

## EvoResist-AI — `prediction.py`
Held-out random forest + permutation importance, **methodology unchanged from v0.1** (not
redesigned). The only code change: prediction after fitting uses `n_jobs=1` (identical output,
faster). Exploratory only; predictive importance is not mechanism.
