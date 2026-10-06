# Next real-data validation stage (commands; none executed in this build)

0. **Verify the discovery breakpoint first.** Read Kayama et al. 2023 (Nat Commun, PMID 38052776) Methods + Supplementary and the
   2024 erratum (Nat Commun 15:782). Only if the S/I/R values are stated there, fill `config/amr/config.yaml →
   breakpoints.ciprofloxacin_ecoli_discovery` (`status: VERIFIED`, both numbers, `verified_by`, `source_excerpt`) and repeat for replication.
1. `python scripts/amr_01_fetch_bvbrc.py --dry-run` then
   `python scripts/amr_01_fetch_bvbrc.py --allow-network --raw-dir data/raw/bvbrc --with-metadata` (resumable; needs BV-BRC access; check field names).
2. Genomes → features (EXTERNAL, UNTESTED here): download FASTA from BV-BRC `genome_sequence` (`variantbridge.amr.data_io.contigs_json_to_fasta`),
   call gene presence/absence (panaroo/Roary), point mutations (AMRFinderPlus), optionally unitigs (unitig-caller); assemble
   `features.npz` with `variantbridge.amr.features_io.save_feature_matrix(..., data_scope="real")` (samples = genome_ids, `feature_meta` with gene/aa_ref/aa_pos/aa_alt).
3. `python scripts/amr_02_build_phenotypes.py --raw-dir data/raw/bvbrc --meta-csv data/raw/bvbrc/genome_metadata.csv` → records A5.
4. `python scripts/amr_03_run_association.py --features data/work/amr_features.npz` → A2, A4 (answer-key hash recorded first).
5. `python scripts/amr_04_run_prediction.py --features data/work/amr_features.npz`
6. `python scripts/amr_05_validate_answer_key.py --features data/work/amr_features.npz` → A1, A3.
7. `python scripts/amr_06_build_evidence_cards.py --features data/work/amr_features.npz --top 25`
8. `streamlit run app.py` → EvoResist-AI → *EvoResist v0.2 artifacts*.
Do not edit `answer_key.yaml` after step 4 (script 05/06 refuse a changed hash). Report failures and nulls as results.
