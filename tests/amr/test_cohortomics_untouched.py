"""EvoResist-AI v0.2 must not modify or break the CohortOmics modules delivered in VariantBridge v0.2.
Hashes were taken from the delivered archive (variantbridge_v0.2.zip) before any EvoResist change."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FROZEN = json.loads((Path(__file__).parent / "cohortomics_v0_2_hashes.json").read_text())


def test_cohortomics_modules_and_scripts_are_byte_identical():
    changed = [p for p, h in FROZEN.items() if hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != h]
    assert changed == [], f"CohortOmics files modified: {changed}"


def test_frozen_list_covers_the_core_modules():
    for m in ("association.py", "qc.py", "eqtl.py", "expression_score.py", "population_structure.py", "validation.py", "manifest.py", "prediction.py"):
        assert f"src/variantbridge/{m}" in FROZEN
