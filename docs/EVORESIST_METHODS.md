# EvoResist-AI v0.2 — METHODS

Frozen scientific basis: `docs/source_reports/report_evoresist_ai_audit.md` + `report_evoresist_cohort_addendum.md`
(addendum wins). Items marked **[user-directed]** were required by the implementation brief and are
*additions/clarifications*, not silent redesigns; each is listed in CHANGELOG. Nothing here has been run on real data.

## Question and system
E. coli × ciprofloxacin. Discovery PMID 38052776 · replication PMID 34485958 · exploratory = remaining eligible
BV-BRC data (descriptive only). A genome in more than one tier goes to the highest tier only; the 38219757 ⊃ 28720578
overlap keeps the 38219757 record. Tier mapping lives in `config/amr/config.yaml`.

## MIC censoring (user-critical) — `amr/phenotypes.py`
Raw string, operator, boundary, censoring direction (left/right/exact/invalid), inclusivity, `log2_mic_exact` (exact
records only; NaN for censored), `log2_boundary`, the interval `[log2_lower, log2_upper]` and an ordinal code are
**separate columns**. A boundary is never written into the exact-measurement column. Conflicting
`sign`/`value` vs string representations are *invalid*, not guessed. Off-grid dilutions are flagged, never altered.

## Primary association model [user-directed] — `amr/association.py`
Interval-censored Gaussian ("Tobit-type") regression of log2 MIC per feature, lineage fixed effects (ST groups):
effect (log2-MIC difference per feature present), SE, CI, z and **Wald p from the same fitted model**. Newton
iterations with analytic gradient/Hessian (verified against numerical derivatives in tests). A likelihood-ratio p
is reported only as a diagnostic column. BH-FDR (primary) and Bonferroni over the features actually tested; an
unadjusted fit is reported beside the adjusted one. The model takes **no annotation, pathway or answer-key input**;
feature frequency filters (`min_isolates_with_feature`, `maf_min`) are thresholds from config and every filtered
feature is listed with its reason. Only the discovery tier is accepted. Replication re-tests pre-specified discovery
hits with the same model; hit replicated = same direction and p < 0.05 / (#hits).
*Why not pyseer:* audit §5 names pyseer, which treats the phenotype as continuous (naive for censored MIC) and could
not be installed in the build environment. `association.run_pyseer` wraps it for unitig/gene-PA screening with a binary (verified-breakpoint) phenotype, or as an
explicitly acknowledged *censoring-naive sensitivity analysis only* (output labelled; cards may not use it). Its command flags and
output columns are UNVERIFIED (tested with an injected fake runner and a format fixture). Unitig GWAS is therefore
**NOT RUN**; the native model runs on gene presence/absence, SNP and amino-acid-substitution features.

## Breakpoints (never inferred)
`config.yaml → breakpoints` is `UNVERIFIED` for both cohorts. `assign_binary`, the binary logistic association and
any S/I/R statistic raise `BreakpointUnverified` unless status is `VERIFIED` **and** `verified_by`, `source_excerpt`
and both numeric values are filled. What was found for the discovery cohort (Kayama 2023 Methods, retrieved
2026-10-05): broth microdilution (MicroScan WalkAway panels) and "MIC cutoffs … according to the Clinical and
Laboratory Standards Institute 2021". Numeric ciprofloxacin values and the concentration range/censoring convention were
not in the retrieved text; the 2024 erratum was not readable. Intermediate calls whose interval straddles a breakpoint
are `indeterminate`.

## Predictive baseline — `amr/prediction.py`
Breakpoint-free targets: (1) censored log2 MIC with L2-penalised interval-censored regression (metric: interval-aware
MAE); (2) an *extreme-class contrast* (left-censored vs right-censored with strict separation of bounds) with
elastic-net logistic (metric: AUC) — **not** an S/R classifier [user-directed, because breakpoints are unverified].
Nested CV, outer *and* inner folds are lineage-blocked (whole STs held out); feature screening and hyper-parameters
are chosen inside training folds. Random-split CV exists only as a labelled leakage demonstration. Primary external
validation is discovery → replication; exploratory data are rejected by the metric functions.

## Interpretation layer — `amr/annotation.py`, `amr/evidence.py`
`answer_key.yaml` is pre-registered, hashed before association and used only for validation labels and card labels.
Candidate ranking is purely statistical (BH q, |z|); annotation columns are appended unchanged and never reorder it;
there is no composite score. Evidence Cards (schema v1.1 = supplied v1 + `causal_status`, `replication_evidence`,
cohort/censoring provenance; additive only) carry four independent dimensions and are rejected if they use causal or
mechanistic language without experimental evidence, claim "novel determinants", merge dimensions in one statement,
claim lineage-adjusted evidence that was not lineage-adjusted, cite random-split prediction, or (synthetic) make any
biological claim.

## Pre-registered answer-key criteria (thresholds in config with provenance)
A1: ≥1 gyrA/parC QRDR cluster FDR-significant and within the top 20 clusters (linked features, r² ≥ 0.95, count once),
positive direction; and lineage-blocked extreme-contrast AUC ≥ 0.70 (audit says "materially above chance"; the
number 0.70 was fixed here before any real-data run). Secondary determinants (qnr, aac(6')-Ib-cr, acrR, marR, soxR)
are reported, not scored; an honest null is a result. CRyPTIC rpoB/katG benchmark: **not implemented**.

## Compute / external tools (UNTESTED here)
unitig-caller, panaroo/Roary, AMRFinderPlus, mash/ANI, Gubbins and pyseer were not available/run. Feature matrices
must be supplied by those tools in the documented npz format (`amr/features_io.py`).
