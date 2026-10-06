#!/usr/bin/env python
"""OFFLINE, RESUMABLE BV-BRC acquisition (genome_amr records, then genome metadata). Never run by Streamlit.

    python scripts/amr_01_fetch_bvbrc.py --dry-run
    python scripts/amr_01_fetch_bvbrc.py --allow-network --raw-dir data/raw/bvbrc
    python scripts/amr_01_fetch_bvbrc.py --allow-network --with-metadata

Pages already on disk are never re-downloaded. If the network is blocked the script exits with code 2 and the
stage is NOT RUN. BV-BRC field names / RQL paging / metadata query syntax are UNVERIFIED in the build environment."""
import argparse
import sys
import urllib.request

import _amr_common as K
from variantbridge.amr import data_io as D


def http_get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--raw-dir", default="data/raw/bvbrc")
    ap.add_argument("--allow-network", action="store_true")
    ap.add_argument("--with-metadata", action="store_true", help="also fetch genome metadata for all genome_ids in the cached AMR pages")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    cfg = K.C.load_config(a.config)
    if a.dry_run:
        print("would query:", D.API_AMR + "?" + D.build_rql(cfg, 0))
        print("raw dir:", a.raw_dir, "| page size:", D.PAGE)
        return
    try:
        pages = D.fetch_genome_amr(cfg, a.raw_dir, http_get, a.allow_network)
        print(f"{len(pages)} genome_amr page(s) cached in {a.raw_dir}")
        if a.with_metadata:
            rec = D.load_cached_records(pages)
            meta = D.fetch_genome_metadata(rec["genome_id"].astype(str).tolist(), a.raw_dir, http_get, a.allow_network)
            meta.to_csv(f"{a.raw_dir}/genome_metadata.csv", index=False)
            print(f"metadata for {len(meta)} genomes -> {a.raw_dir}/genome_metadata.csv")
    except (D.DataNotAvailable, OSError) as e:
        print(f"NOT RUN: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
