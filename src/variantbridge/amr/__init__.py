"""EvoResist-AI v0.2: microbial GWAS / AMR pilot (E. coli - ciprofloxacin).

Science modules only; no Streamlit import, no network access at import time. Long-running
real-data operations are offline scripts (scripts/amr_*.py) that write compact manifest-validated
artifacts under artifacts/amr/. Methodology: report_evoresist_ai_audit.md + addendum (frozen).
"""
AMR_VERSION = "0.2.0"
