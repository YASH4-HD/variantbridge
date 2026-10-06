import os
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "scripts", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


@pytest.fixture(scope="session")
def plink2_bin():
    cand = os.environ.get("PLINK2") or shutil.which("plink2")
    if not cand:
        pytest.skip("plink2 not available (set PLINK2 env var or put plink2 on PATH); V5 benchmark NOT RUN")
    return cand
