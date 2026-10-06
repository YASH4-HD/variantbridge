# TEST RESULTS (executed 2026-10-05, Linux, Python 3.13, streamlit 1.65)

## EvoResist-AI v0.2 build - full suite
| Run | Result |
|---|---|
| With PLINK2 (built from source, `PLINK v2.0.0-b.1-devNL`, commit a2b85291) | **221 passed** in 124 s |
| Without PLINK2 | **211 passed, 10 skipped** in 105 s (all skips are plink2-dependent CohortOmics tests; V5 NOT RUN in that mode) |

Of these, 71 are the CohortOmics v0.2 tests (unchanged) and 150 are new (`tests/amr/`): MIC parsing/censoring, config provenance,
cohort tiers and duplicate exclusion, QC gates, censored-model derivative/OLS-equivalence/bias tests, association-output consistency,
discovery-only guards, pyseer wrapper (fake runner), lineage-blocked CV and the leakage lesson, Evidence Card schema + claim
restrictions, validation registry, data_io (no network), script plumbing (amr_00..06), Streamlit viewer smoke tests, and a hash test
proving the CohortOmics modules/scripts are byte-identical to the delivered v0.2 archive.

**All EvoResist tests use SYNTHETIC or fixture data.** They verify software behaviour (including recovery of PLANTED effects), not
biology. The only real-source check is the discovery-breakpoint literature read (A6 = UNVERIFIED).

## Real vs synthetic
| Item | Data |
|---|---|
| Discovery breakpoint/AST standard (A6) | REAL source text (partial): CLSI 2021 stated; numeric values NOT found; erratum NOT read -> UNVERIFIED |
| Everything else in EvoResist (parsing, cohorts, association, CV, cards, scripts, viewer) | SYNTHETIC / fixtures |
| PLINK2 --glm V5 benchmark (CohortOmics) | SIMULATED inputs vs real plink2 binary (see below) |

## NOT RUN / UNVERIFIED (EvoResist)
BV-BRC download and field names; real phenotype build (A5); real association (A2), answer-key recovery (A1), predictive
validation (A3), replication (A4); any breakpoint-dependent analysis; pyseer/unitig GWAS; AMRFinderPlus annotation; species ANI;
Gubbins/mash; CRyPTIC benchmark; literature review for candidates; all wet-lab validation.

---
## Earlier CohortOmics results (unchanged)

## pytest
| Run | Result |
|---|---|
| With PLINK2 (built from source, `PLINK v2.0.0-b.1-devNL`, commit a2b85291, see tools/) | **71 passed** in 68 s |
| Without PLINK2 | **61 passed, 10 skipped** (all skips are plink2-dependent; V5 NOT RUN in that mode) |

Tests use SIMULATED data only. They verify code correctness and plumbing, not biology.

## V5 (PLINK2 benchmark) — SIMULATED inputs only
Executed against real `plink2 --glm` (400 simulated variants, linear and logistic, identical subset/model).
Max relative difference: ≈4.4e-6–4.8e-6 on effect/SE; ≈1.9e-6 on −log10 p; observation counts identical;
3 variants re-oriented (A1=REF→ALT). Pre-registered tolerances (1e-4 rel; 5e-3 abs) were met.
This was run with `--no-record`: **`artifacts/VALIDATION_STATUS.json` is not shipped**, so V5 shows
NOT RUN in the app until you run it on a real subset. It shows arithmetic equivalence only.

## NOT RUN (require real data / unavailable in build environment)
| ID | Status |
|---|---|
| V1 population structure (real 1000G) | NOT RUN |
| V2 X-het sex check (real chrX; `--x-vcf` path untested) | NOT RUN |
| V3 permuted-phenotype λ on real GEUVADIS | NOT RUN |
| V4 EUR/YRI answer-key replication (release-blocking) | NOT RUN — no replication claim |
| V5 on a real chr22 subset | NOT RUN |
| V6 Genetic Expression Score on real eQTLGen/GEUVADIS | NOT RUN |
| `plink2 --pca` engine | UNTESTED (build used NO_LAPACK=1) |
| Real-file header/format assumptions | UNVERIFIED (docs/DATASETS.md) |
