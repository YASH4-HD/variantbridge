# Addendum: Cohort Selection for the EvoResist-AI E. coli–Ciprofloxacin Pilot

**Date:** 5 October 2026. Companion to `report_evoresist_ai_audit.md`. All BV-BRC numbers queried live from `genome_amr`/`genome` APIs (laboratory-evidence records only) on 5 Oct 2026. Items not verifiable from API or primary-paper abstracts are marked **UNVERIFIED**.

---

## Part 1 — Cohort-Selection Addendum

### 1.1 Per-cohort characterisation

| Dimension | **PMID 38052776** — Kayama et al., Nat Commun 2023 | **PMID 34485958** — Cambridge Univ. Hospitals study, Lancet Microbe 2021 | **PMID 38219757** — Lancet Microbe 2024 | **PMID 28720578** — Kallonen et al., Genome Res 2017 |
|---|---|---|---|---|
| Isolates (cipro records / unique genomes) | 2,429 / 2,429 | 1,092 / 1,092 | 4,286 / 4,286 | 1,504 / 1,504 |
| Records with quantitative MIC | **2,429 (100%)** | 1,092 (100%, heavily censored) | 1,553 (36%) | **0** (methods are dilution-based, but no measurement values in BV-BRC — categorical R/S only in the aggregate) |
| AST method | Uniform quantitative MIC (single national standard) | Broth dilution | MIC 1,553 + disk diffusion 2,733 | Agar dilution 1,094 + broth dilution 410 |
| Ciprofloxacin units / values | mg/L; censored scale `<=0.25` (764), 1 (253), 2 (23), 4 (14), `>4` (1,375) | mg/L; `<=0.25` (606), `>=4` (446), 0.5–2 sparse (34) | mg/L; continuous range 0.008–256 in MIC subset | not available (R/S only: 1,226 S / 12 I / 266 R) |
| Breakpoint standard | **UNVERIFIED** at abstract level (protocol in paper methods; an erratum exists, Nat Commun 2024;15:782) | **UNVERIFIED** | **UNVERIFIED** | **UNVERIFIED** |
| Sequencing availability | WGS for all (open BV-BRC/NCBI); long-read subset in study | WGS for all (open) | WGS for all (open) | WGS for all (open) |
| Geography | Japan (national) | England (two haematology wards, Cambridge) | Norway 3,192 + UK 1,094 | UK 1,504 |
| Collection period | 2019 (dated per isolate) | 2015 (dated per isolate) | **UNVERIFIED** (longitudinal; dates not in BV-BRC) | **UNVERIFIED** (dates not in BV-BRC) |
| Lineage diversity (unique MLST) | 235 STs | 96 STs | 467 STs | 227 STs |
| Genome↔phenotype join key | BV-BRC `genome_id` (verified in pilot) | same | same | same |
| Duplicate isolates across studies | none (disjoint from all others) | none (disjoint) | **shares exactly 1,094 genomes with PMID 28720578** (its UK arm = the 2017 UK collection) | contained within 38219757's UK arm |
| Accessibility | Open (BV-BRC/NCBI); open paper | Open | Open | Open |
| Sampling frame caveat | Enriched for 3GC-resistant / reduced-carbapenem-susceptibility isolates (study design) — not a random population sample | Hospital-specific (single centre, two wards) | Mixed methods within one study | — |

### 1.2 Overlap matrix (shared `genome_id`s)

| | 38052776 | 34485958 | 38219757 | 28720578 |
|---|---|---|---|---|
| **38052776** | — | 0 | 0 | 0 |
| **34485958** | 0 | — | 0 | 0 |
| **38219757** | 0 | 0 | — | **1,094** |
| **28720578** | 0 | 0 | **1,094** | — |

Total unique genomes across the four studies: 8,217. The 38219757↔28720578 overlap means these two must never be split across discovery/replication; treat 28720578 as a subset of 38219757.

### 1.3 Recommended arrangement

1. **Primary discovery cohort — PMID 38052776 (Japan NIID national surveillance), n = 2,429.**
   The only cohort with 100% standardized quantitative ciprofloxacin MIC from a single laboratory standard, per-isolate dates, 235 STs of lineage diversity, and full open WGS. This is the cohort that satisfies "standardized quantitative MIC methodology".
   *Declared limitations:* (a) the sampling frame enriches for 3GC-resistant/reduced-carbapenem-susceptibility isolates, so prevalence-based claims do not generalise — association claims remain internal to the cohort; (b) MIC is censored at `<=0.25` and `>4` — analyse as censored log2 MIC or ordinal, not naive linear; (c) exact breakpoint/AST standard to be extracted from the paper methods (**UNVERIFIED** here); (d) an erratum exists — record it in the provenance manifest.
2. **Independent replication cohort — PMID 34485958 (Cambridge, England, 2015), n = 1,092.**
   Disjoint from discovery by construction (0 shared genomes), different country, period, and study; broth-dilution MIC with censoring (`<=0.25`/`>=4`) supports replication of association direction and effect presence, though with reduced resolution. Secondary replication option: the 459-genome Norway MIC subset of 38219757 (disjoint from both), too small alone but usable as a third point.
3. **Exploratory / generalisation set — the full BV-BRC aggregate (10,904 unique genomes) minus all discovery and replication genomes.**
   Includes the mixed-methods remainder (disk diffusion, agar dilution, uncalled phenotypes). Used only for: frequency spectra, structure description, and qualitative consistency checks. **Never** pooled into discovery or used for headline association or predictive metrics. The 1,094 duplicated genomes are counted once (keep the 38219757 record; drop the 28720578 duplicate).

**Rule for all three tiers:** a genome that appears in more than one tier is assigned to the highest tier only; tier membership is recorded in the isolate manifest.

---

## Part 2 — Implementation Specification delta for Claude

These amendments modify the Implementation Specification in `report_evoresist_ai_audit.md` §10; everything else there stands.

1. **QC/frequency thresholds are configuration, not constants.** All thresholds (CheckM completeness/contamination, genome-length range, MAF/rare-variant frequency cutoff, MIC censoring handling) move from hard-coded defaults into `config.yaml` under a `thresholds:` block. Each entry must carry a `provenance` field with one of: `source_derived` (justified from the source study's own QC, with citation) or `pre_registered` (a project design choice, recorded with the analysis date and rationale). The pipeline must print the provenance class of every threshold in the methods report. No threshold may be presented as immutable.
2. **Association is unbiased and genome-wide.** The association module (`association.py`) runs over the complete unitig/gene-PA/SNP feature space with no biological pre-filtering. Efflux-pathway membership, AMR-database annotation, and any biological classification live only in `annotation.py`/`evidence.py` as *interpretation* layers and must not enter feature selection, the association model, or the statistical priority score. The candidate-prioritisation score must be decomposable in output (statistical component reported separately from any annotation component), so a reviewer can see the pure association ranking.
3. **Known efflux determinants are pre-registered controls, not a search space.** The answer-key set (gyrA/parC QRDR, acrR/marR, qnr, etc.) is fixed in a version-controlled `answer_key.yaml` before the association run and used only for (a) pipeline validation (§6 known-answer pass criteria) and (b) labelling of Evidence Cards. Candidate discovery must not be restricted to efflux genes or any gene set; hits outside the answer key are legitimate discoveries and must not be down-weighted for lacking efflux annotation.
4. **Four separate evidence dimensions, never merged.** Evidence Cards and the prioritisation output must carry four independent fields — `association_evidence` (statistics from the fitted model), `predictive_contribution` (outer-fold, lineage-blocked CV only), `biological_plausibility` (annotation/literature interpretation), and `causal_status` (empty/"not established" unless allele-replacement or equivalent experimental evidence exists). No composite score may sum or average across these dimensions; the UI must render them as separate, labelled sections. Any statement combining them (e.g., "predictive and plausible, therefore causal") is a schema violation.

**Cohort-tier interface change:** `data_io.fetch_phenotypes(config)` gains a `cohort_tier` column (`discovery` | `replication` | `exploratory`) derived from the PMID→tier mapping in config; `association.run_pyseer` accepts only `discovery`; `prediction.py` reports discovery→replication as the primary validation and treats exploratory data as descriptive only.
