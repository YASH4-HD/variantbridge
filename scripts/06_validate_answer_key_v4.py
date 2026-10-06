#!/usr/bin/env python3
"""V4 (RELEASE-BLOCKING): compare re-derived GEUVADIS cis-eQTL results with the published answer key.

Run scripts/05 first with --key-files so pair statistics for published gene-SNP pairs are retained.

    python scripts/06_validate_answer_key_v4.py --population EUR \\
        --key-best data/raw/EUR373.gene.cis.FDR5.best.rs137.txt.gz \\
        --key-all  data/raw/EUR373.gene.cis.FDR5.all.rs137.txt.gz

Nothing is recorded as PASSED unless this script actually executes on real data. Pre-registered
criterion: same-direction rate >= 0.90 over >= 100 compared published best pairs. The allele-
orientation convention of the published rho is not assumed (see validation.v4 flags).
"""
import argparse
import json
from pathlib import Path

import pandas as pd

import _common  # noqa: F401
from variantbridge import capability, validation
from variantbridge.prep import parse_answer_key


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--population", required=True, choices=["EUR", "YRI"])
    ap.add_argument("--key-best", required=True)
    ap.add_argument("--key-all", required=True)
    ap.add_argument("--work-dir", default="data/work/geuvadis")
    ap.add_argument("--artifact-dir", default="artifacts")
    ap.add_argument("--no-record", action="store_true")
    a = ap.parse_args()
    ad = Path(a.artifact_dir)
    mine_best = pd.read_csv(ad / f"geuvadis_eqtl_{a.population}.csv.gz")
    pairs_path = Path(a.work_dir) / f"pairs_{a.population}.csv.gz"
    if not pairs_path.exists():
        raise SystemExit("pairs file missing: run scripts/05 with --key-files POP=<answer key> first")
    mine_pairs = pd.read_csv(pairs_path)
    res = validation.v4_answer_key_replication(a.population, mine_best, mine_pairs, parse_answer_key(a.key_best),
                                               parse_answer_key(a.key_all))
    out = ad / f"v4_{a.population}.json"
    out.write_text(json.dumps(res, indent=2, default=str))
    if not a.no_record:
        prev = capability.load_status(ad)["V4"]
        per = prev.get("populations", {}) if prev.get("status") != capability.NOT_RUN else {}
        per[a.population] = res
        statuses = [r["status"] for r in per.values()]
        overall = capability.PASSED if len(per) == 2 and all(s == capability.PASSED for s in statuses) else (
            capability.FAILED if any(s == capability.FAILED for s in statuses) else capability.NOT_RUN)
        capability.record_status("V4", {"status": overall, "executed_at": validation.now(), "populations": per,
                                        "note": "PASSED requires BOTH EUR and YRI runs to pass"}, ad)
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
