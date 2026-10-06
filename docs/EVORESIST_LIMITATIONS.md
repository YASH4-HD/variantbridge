# EvoResist-AI — LIMITATIONS

- **No real-data run has been executed.** All executed tests use SYNTHETIC or fixture data. A planted determinant being
  recovered shows the software works, not that EvoResist-AI recovers real resistance biology.
- Discovery breakpoint **UNVERIFIED**; no S/I/R, binary-logistic or clinical-category analysis exists until it is verified.
- Gaussian assumption on log2 MIC; exact dilution values are treated as exact (not ±½-dilution intervals).
- Extreme-contrast AUC is not a clinical classifier; class definitions differ between cohorts (`>4` vs `>=4`).
- Fixed-effect ST adjustment cannot remove all structure (recombination, accessory-genome linkage); linked features are not separable.
- Unitig GWAS, pyseer (LMM), AMRFinderPlus/CARD annotation, species ANI check, Gubbins, mash tree, homoplasy analysis and the
  CRyPTIC TB benchmark are **NOT RUN / NOT IMPLEMENTED**.
- Discovery sampling frame is enriched for 3GC-resistant isolates; single-cohort associations are internal to it.
- No clinical, treatment or causal claim is supported; allele replacement is the only route to causal language.
- Literature support in cards is empty by design (no fabricated citations); a literature review is a listed missing evidence item.
