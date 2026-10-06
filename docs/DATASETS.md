# DATASETS

Source of truth: `config/datasets.json` (URLs copied from the CohortOmics audit, 2026-10-05).
Items the audit did not verify carry `UNVERIFIED`. **No raw data is in this repository**; nothing is
downloaded by the app. Downloads: `python scripts/01_download_data.py --list|--dry-run|--yes`.

| Dataset | Role | Notes |
|---|---|---|
| 1000 Genomes phase 3 (IGSR), GRCh37 | genotype backbone: QC, relatedness, PCA (V1, V2) | chr22 VCF ≈206 MB size verified by audit. Panel filename UNVERIFIED. |
| NYGC 30x 1000G (GRCh38) | optional | per-chromosome filename UNVERIFIED; pass `--url`. |
| GEUVADIS (E-GEUV-1), GRCh37, Gencode v12 | association engine (V3, V4, V6 target) | 462 samples; genotypes per chromosome (chr22 ≈651 MB), PEER resk10 expression ≈91 MB. Open, no restrictions. |
| GEUVADIS answer keys EUR373/YRI89 gene cis FDR5 | V4 | file names verified; **column headers UNVERIFIED** |
| eQTLGen | V6 weights | column names UNVERIFIED |
| PGP-UK | optional height score demo (portfolio only) | small; not implemented as validated |
| HipSci | partially open | not used |
| openSNP | rejected by audit | not used |
| Controlled-access (UK Biobank, dbGaP, …) | out of scope | no disease PRS possible without them |

## UNVERIFIED format assumptions (must be checked on first real run)
1. Answer-key and expression-matrix header names (`prep.resolve_columns` fails loudly if absent).
2. Whether the answer-key `Coord`/position is the gene TSS.
3. ID-converter column order (Ensembl ↔ sample IDs).
4. eQTLGen column names.
5. 1000G panel filename.
6. Genome build agreement between GEUVADIS (GRCh37) and eQTLGen (GRCh37 per audit; confirm).

Scripts raise `PrepError` rather than guessing when a column cannot be resolved.

## Licence / citation
1000G: Fort Lauderdale; cite doi:10.1038/nature15393, doi:10.1093/nar/gkz836. GEUVADIS:
Lappalainen et al. 2013, Nature. See audit references.
