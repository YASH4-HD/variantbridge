# CHANGELOG

## 0.2.1 - EvoResist-AI v0.2 (CohortOmics modules byte-identical to 0.2.0; enforced by a hash test)
### Added
- `src/variantbridge/amr/`: config (threshold provenance), phenotypes (explicit MIC censoring), cohorts (tiers, duplicate
  exclusion, separation guards), data_io (offline resumable BV-BRC fetch, injectable network), qc, structure, annotation
  (interpretation only), association (interval-censored model, pyseer wrapper), prediction (lineage-blocked nested CV),
  evidence (schema-validated cards + claim restrictions), validation (pre-registered answer key, status registry), reporting.
- `config/amr/`: `config.yaml`, pre-registered `answer_key.yaml`, supplied Evidence Card schema v1 and additive v1.1
  (`causal_status`, `replication_evidence`, cohort/censoring provenance), `dataset_provenance.csv`.
- Scripts `amr_00` (synthetic inputs) and `amr_01`-`amr_06`; `amr_view.py` viewer; `docs/EVORESIST_*.md`.
- Tests: `tests/amr/` (parsing, censoring, cohorts, duplicates, QC, association consistency, lineage-blocked CV, cards, claim guards, scripts, viewer, CohortOmics-unchanged).
### Decisions needing review (user-directed additions, not silent redesigns)
- Primary association model is interval-censored Gaussian (native), not pyseer-on-continuous-MIC; pyseer is wrapped but only for
  binary-verified-breakpoint or acknowledged censoring-naive sensitivity use, and could not be run here.
- Predictive targets are breakpoint-free (censored log2 MIC; extreme-class contrast) because the breakpoint is unverified.
- `answer_key_auc_min = 0.70` and `top_n_hits = 20` were pre-registered (audit gave no numbers).
- `app.py`: one added sidebar option and a 4-line dispatch to `amr_view.py`; the CohortOmics path is unchanged.
### Fixed during build
- A threshold missing `provenance` was silently skipped by the config walker (now rejected); BV-BRC paging limit/offset could disagree (now consistent); script `--phenotypes` path ignored its directory.
### NOT RUN / UNVERIFIED
- Every real-data stage; breakpoint for both cohorts; pyseer/unitig GWAS; AMRFinderPlus; species ANI; BV-BRC field names/paging; CRyPTIC benchmark.

## 0.2.0
### Fixed (scientific)
- v0.1 paired logistic coefficients with a point-biserial p-value (inconsistent). Association now
  uses one fitted model (OLS t-test / IRLS logistic Wald) for effect, SE, CI, p. Release-blocking
  consistency check + V5 PLINK2 benchmark harness.
- PCA/QC: genotype PCA on LD-pruned variants with deterministic signs; exact HWE; KING-robust
  kinship; no silent mean-imputation of genotypes in the human branch (v0.1 AMR path preserved).
### Added
- Modular package `src/variantbridge/`; Streamlit app is a thin UI.
- Offline scripts 01–03, 05–08; manifests with SHA-256; validated-artifact loader.
- GEUVADIS-native cis-eQTL scan with per-gene permutation and BH; V3/V4 validations.
- Genetic Expression Score (eQTLGen → GEUVADIS), never called PRS; V6 with permutation null.
- Validation status registry (absent ⇒ NOT RUN); capability tiers; evidence cards.
- Tests, docs (METHODS, DATASETS, LIMITATIONS, INTERVIEW_NOTES).
### Changed
- V1 redefined: leading PCs recover broad superpopulation structure, quantified.
- EvoResist-AI prediction moved unchanged to `prediction.py` (only post-fit `n_jobs=1`; same output).
- Requirement: streamlit>=1.50 (`width="stretch"`).
### Not done / NOT RUN
- No real-data validation executed (V1–V4, V6); V5 only on simulated data. Script `04` intentionally
  absent (eQTLGen weight prep is inside script 08).
