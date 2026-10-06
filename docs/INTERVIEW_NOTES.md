# INTERVIEW NOTES (truthful scope)

What this project is: a transparent, test-backed research/portfolio codebase for open-data human
cohort genomics, with explicit validation gates and claim control.

Say:
- "I built the QC/PCA/association/eQTL modules and a validation framework with pre-registered
  criteria; real-data validations are specified and scripted, and I report which have been run."
- "The v0.1 app paired a logistic coefficient with a point-biserial p-value; I replaced it with a
  single-model engine and benchmarked it against PLINK2 on identical inputs."
- "The expression score is a Genetic Expression Score, not a PRS; a null result is reported
  as a result."

Do not say:
- That GEUVADIS has been replicated, or that V1–V4/V6 passed, unless `VALIDATION_STATUS.json`
  shows an executed pass from your own run.
- That you have biobank, controlled-access or clinical PRS experience.
- That EvoResist-AI is a validated AMR predictor: it is a developing, exploratory direction.

Likely questions: why separate EUR/YRI; why permutation over Bonferroni; what leakage would
inflate a score; why not hard-code long-range LD regions; what V4 near-zero agreement means
(allele orientation).
