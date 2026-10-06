# artifacts/

Compact outputs of the offline scripts, each with a `*.manifest.json` (SHA-256, inputs, versions,
parameters, deviations). The app refuses an artifact whose manifest is missing/invalid or whose
checksum does not match.

`VALIDATION_STATUS.json` is written **only** by code that actually executed a validation
(`capability.record_status` requires `executed_at`). If a validation id is absent, its status is
**NOT RUN**. Nothing in this directory is shipped pre-filled: this repository ships with no
real-data results.
