# VariantBridge v1

Two Streamlit branches:
- CohortOmics — human cohort/population genomics
- EvoResist-AI — pathogen/AMR genomics

Run:
`pip install -r requirements.txt`
`streamlit run app.py`

Input genotype/features CSV: `sample_id` + numeric feature columns.
Phenotype CSV: `sample_id` + `phenotype`.
Optional metadata: `variant`, `chrom`, `pos`, annotations such as `known_amr`.

This is a transparent research/portfolio prototype. Built-in data are simulated. Lightweight association and score modules are exploratory and should be reproduced with validated domain-specific tools for research publication.
