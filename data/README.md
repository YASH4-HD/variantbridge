# data/

Raw and intermediate genomic data are **never committed** (see `.gitignore`).

- `data/raw/`  : downloads from `scripts/01_download_data.py` (+ `DOWNLOAD_LOG.json` with SHA-256).
- `data/work/` : intermediate PLINK2 filesets / per-chromosome blocks.

The Streamlit app never reads these directories and never downloads anything.
Only compact, manifest-validated files in `artifacts/` are loaded by the app.
