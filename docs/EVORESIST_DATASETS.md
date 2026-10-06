# EvoResist-AI — DATASETS

Source: `config/amr/dataset_provenance.csv` (supplied, audit 2026-10-05) and the cohort addendum. **No data is in this
repository.** Acquisition is offline/resumable (`scripts/amr_01_fetch_bvbrc.py`); the BV-BRC API was unreachable from the
build environment, so every count below is *as reported by the addendum*, not re-queried.

| Tier | Study | n (addendum) | Notes |
|---|---|---|---|
| discovery | PMID 38052776, Kayama 2023 | 2,429 | 100 % quantitative MIC; censored `<=0.25`(764) `1`(253) `2`(23) `4`(14) `>4`(1,375); sampling frame enriched for 3GC-resistant isolates; erratum Nat Commun 2024;15:782 not read |
| replication | PMID 34485958, Lancet Microbe 2021 | 1,092 | broth dilution, `<=0.25`(606) `>=4`(446), sparse middle |
| exploratory | everything else eligible (incl. 38219757; 28720578 collapsed into it) | up to ~10,904 unique genomes in total | descriptive only |

Reconciliation flags (UNVERIFIED): the discovery paper text reports 4195 isolates tested / 5143 sequenced collected 2019–2020,
the addendum 2,429 records collected 2019 — reconcile with real data (script amr_02 reports A5 against the addendum counts).
BV-BRC field names (`measurement`, `measurement_sign`, `measurement_value`, `testing_standard`, `pmid`, MLST/CheckM names), RQL
paging and the metadata query are **UNVERIFIED** against the live API; loaders fail loudly (`FIELD_ALIASES`).
Licence: BV-BRC public genomes; cite BV-BRC + source PMIDs; do not re-host the aggregate; CARD data files are non-commercial
(not bundled); AMRFinderPlus licence text UNVERIFIED.
