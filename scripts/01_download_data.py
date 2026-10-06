#!/usr/bin/env python3
"""Offline, resumable dataset download with SHA-256 logging. NEVER run by the Streamlit app.

    python scripts/01_download_data.py --list
    python scripts/01_download_data.py --item geuvadis_expression_resk10 --dry-run
    python scripts/01_download_data.py --item 1000g_phase3_chr22_vcf --yes

Large items require --yes after you have read the size. Each download appends an entry
(url, size, sha256, access date) to data/raw/DOWNLOAD_LOG.json for the provenance manifests.
Raw data are git-ignored and must never be committed.
"""
import argparse
import datetime as dt
import json
import sys
import urllib.request
from pathlib import Path

import _common  # noqa: F401
from variantbridge.manifest import file_sha256

CFG = json.loads((_common.ROOT / "config" / "datasets.json").read_text())


def resolve_url(name: str, item: dict, chrom: str | None, override: str | None) -> str:
    if override:
        return override
    if "url" in item and item["url"].startswith("http"):
        return item["url"]
    if "path" in item:
        return CFG["geuvadis_base"] + item["path"]
    if "path_template" in item:
        if not chrom:
            raise SystemExit(f"{name} needs --chrom")
        return CFG["geuvadis_base"] + item["path_template"].format(chrom=chrom)
    raise SystemExit(f"{name}: no downloadable URL configured ({item.get('verified_by')}). Provide --url.")


def download(url: str, dest: Path, chunk: int = 1 << 20) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    start = part.stat().st_size if part.exists() else 0
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-"} if start else {})
    with urllib.request.urlopen(req) as r, open(part, "ab" if start and r.status == 206 else "wb") as fh:
        while True:
            b = r.read(chunk)
            if not b:
                break
            fh.write(b)
    part.rename(dest)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--item")
    ap.add_argument("--chrom")
    ap.add_argument("--url", help="override URL (required for items marked UNVERIFIED)")
    ap.add_argument("--dest", default="data/raw")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--yes", action="store_true", help="confirm download of large files")
    a = ap.parse_args()
    items = {k: v for k, v in CFG["items"].items()}
    if a.list or not a.item:
        for k, v in items.items():
            print(f"{k:42s} ~{v.get('approx_size_mb', '?')} MB  [{v['verified_by']}]")
        return
    if a.item not in items:
        sys.exit(f"unknown item {a.item}")
    item = items[a.item]
    url = resolve_url(a.item, item, a.chrom, a.url)
    size = item.get("approx_size_mb")
    if isinstance(size, dict):
        size = size.get(a.chrom, "?")
    print(f"item={a.item}\nurl={url}\napprox_size_mb={size}\nlicense={item.get('license', 'see audit')}\nverified_by={item['verified_by']}")
    if a.dry_run:
        print("dry run: nothing downloaded")
        return
    if not a.yes and (not isinstance(size, (int, float)) or size > 50):
        sys.exit("size >50 MB or unknown: re-run with --yes to confirm")
    dest = Path(a.dest) / url.rsplit("/", 1)[-1]
    download(url, dest)
    entry = {"item": a.item, "url": url, "file": str(dest), "bytes": dest.stat().st_size,
             "sha256": file_sha256(dest), "access_date": dt.date.today().isoformat(),
             "license": item.get("license", "see audit"), "verified_by": item["verified_by"]}
    logp = Path(a.dest) / "DOWNLOAD_LOG.json"
    log = json.loads(logp.read_text()) if logp.exists() else []
    log.append(entry)
    logp.write_text(json.dumps(log, indent=2))
    print("downloaded", dest)


if __name__ == "__main__":
    main()
