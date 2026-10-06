# VariantBridge v0.2 (+ EvoResist-AI v0.2)

Two branches in one Streamlit app, with the science in an importable, tested package:

- **CohortOmics** — human cohort/population genomics: QC, relatedness, PCA, association,
  GEUVADIS cis-eQTL, Genetic Expression Score.
- **EvoResist-AI** — pathogen/AMR genomics (exploratory; methodology unchanged from v0.1).

> **Status: research/portfolio code. No real-data validation has been executed in this delivery.**
> V1–V4 and V6 are **NOT RUN**; V5 has only been run on *simulated* inputs. Built-in demo data are
> **SIMULATED**. See `TEST_RESULTS.md`, `docs/LIMITATIONS.md`. Replication of GEUVADIS is not claimed.

## EvoResist-AI v0.2 (E. coli × ciprofloxacin) — status
Implemented in `src/variantbridge/amr/` with offline scripts `scripts/amr_00…06` and a viewer (`amr_view.py`, sidebar:
EvoResist-AI → *EvoResist v0.2 artifacts*). Basis: `docs/source_reports/` (audit + cohort addendum). See
`docs/EVORESIST_METHODS.md`, `EVORESIST_DATASETS.md`, `EVORESIST_LIMITATIONS.md`, `EVORESIST_NEXT_STEPS.md`.

* **No real-data analysis has been run.** Executed so far: unit/plumbing tests on SYNTHETIC data and one literature check
  (discovery breakpoint: **UNVERIFIED**). Passing software tests is not evidence that EvoResist-AI recovers real resistance biology.
* MIC censoring is explicit (raw / boundary / direction / transforms are separate columns); the primary model is
  interval-censored; no naive linear model on `<=`/`>` values is fitted. Breakpoint-dependent analyses refuse to run until verified.
* Cohort tiers: discovery PMID 38052776 · replication PMID 34485958 · exploratory = remainder (descriptive only).
* Four evidence dimensions (`association_evidence`, `predictive_contribution`, `biological_plausibility`, `causal_status`) are
  never merged; `config/amr/answer_key.yaml` is pre-registered and only labels/validates.
* Real-data commands: `docs/EVORESIST_NEXT_STEPS.md`.

## Layout
```
app.py                     thin Streamlit UI (never downloads/processes large data)
src/variantbridge/         science (no Streamlit import)
scripts/                   offline pipeline (01 download, 02 1000G PCA, 03 GEUVADIS prep,
                           05 eQTL scan, 06 V4, 07 V5 PLINK2 benchmark, 08 V6 expression score)
config/datasets.json       audit-verified URLs (UNVERIFIED items flagged)
artifacts/                 compact outputs + manifests + VALIDATION_STATUS.json (absent => NOT RUN)
docs/                      METHODS, DATASETS, LIMITATIONS, INTERVIEW_NOTES
tests/                     pytest suite
tools/                     plink2 build notes
```
(`04` is intentionally absent: eQTLGen weight handling lives in script 08.)

## Quick start (demo, SIMULATED)
```
pip install -r requirements.txt
streamlit run app.py
```
Upload mode: genotype/feature CSV (`sample_id` + numeric columns), phenotype CSV
(`sample_id`, `phenotype`), optional variant metadata (`variant`, `chrom`, `pos`, annotations such as `known_amr`).
Missing genotypes are not silently imputed in the human branch.

## Real-data workflow (offline; needs `requirements-offline.txt` and PLINK2)
```
python scripts/01_download_data.py --list
python scripts/01_download_data.py --dry-run          # then --yes to download
python scripts/02_prepare_1000g_pca.py --vcf <1000G chr22 vcf.gz> --panel <panel> \
       (--long-range-ld-file <tsv> | --no-long-range-exclusion) [--king] [--x-vcf <chrX vcf>]   # V1 (+V2)
python scripts/03_prepare_geuvadis.py --vcf-template <path with {chrom}> --panel <panel> \
       --expression <resk10 matrix> --chroms 22 (--long-range-ld-file ... | --no-long-range-exclusion)
python scripts/05_run_eqtl_scan.py --key-files EUR=<EUR best> ... --v3-permutations 20      # V3
python scripts/06_validate_answer_key_v4.py --population EUR --key-best ... --key-all ...    # V4 (and YRI)
python scripts/07_benchmark_plink2_v5.py --pfile <real subset prefix>                         # V5 real subset
python scripts/08_run_expression_score_v6.py --weights <eQTLGen file> --cohort-disjointness-attested  # V6
```
Then select **Validated artifacts (real data)** in the app. The app refuses artifacts with a missing
manifest or SHA-256 mismatch. Check each script's `--help`; several format assumptions are UNVERIFIED
(docs/DATASETS.md) and the scripts fail loudly rather than guess.

## Tests
```
PLINK2=/path/to/plink2 python -m pytest tests -q   # without PLINK2, plink2-dependent tests SKIP (V5 NOT RUN)
```

## Claim control
Genetic Expression Score is **not** a PRS and not "risk". No disease-trait PRS is possible with open
data. Capability tiers are targets; the validation status shown beside them is what has actually run.
Not a clinical tool.
