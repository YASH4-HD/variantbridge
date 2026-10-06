# LIMITATIONS

- **Nothing has been validated on real data in this delivery.** V1–V4 and V6 are NOT RUN; the
  downloads were blocked in the build environment. V5 was run only against *simulated* inputs
  (see TEST_RESULTS.md). Replication of GEUVADIS answer keys must not be claimed until V4 executes.
- Built-in demo data are SIMULATED and prove only that code runs.
- The cis-eQTL significance procedure (direct permutation + BH) is a design choice and differs
  from the 2013 paper; V4 measures agreement, not procedural identity.
- Permutation p-values have resolution 1/(n_perm+1).
- eQTL results are statistical association with expression, not causal mechanism.
- Genetic Expression Score: LCL (GEUVADIS) vs whole blood (eQTLGen) tissue mismatch; possible
  sample overlap; European-ancestry weight bias; a null result is informative, not a bug.
- No open-data disease PRS; the PRS module is a leakage-guarded tool, demonstrated only on
  simulated cohorts.
- KING-robust kinship and X-heterozygosity sex-check on real 1000G data: the X path (`--x-vcf`)
  is untested; sex inference depends on X pseudo-autosomal handling.
- Long-range LD exclusion regions are not hard-coded (unverified); runs without them are
  flagged as a manifest deviation.
- `plink2 --pca` engine path is untested (the build used `NO_LAPACK=1`); the default python PCA
  engine is tested on simulated data only.
- Logistic regression has no Firth correction; rare variants/small groups are flagged
  (`separation`), not repaired.
- EvoResist-AI random forest is exploratory; lineage/population structure confounds AMR
  prediction and is not corrected.
- Not a clinical or diagnostic tool.
