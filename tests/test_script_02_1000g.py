"""Plumbing test for scripts/02 on a SIMULATED 1000G-like VCF using the python PCA engine (the plink2 --pca
path needs a LAPACK-enabled plink2 and is NOT exercised here). Not V1 evidence."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def write_vcf_and_panel(d, seed=0):
    from variantbridge.demo import simulate_structured_genotypes
    G, lab = simulate_structured_genotypes((50, 50, 50), m=1500, fst=0.1, seed=seed)
    n, m = G.shape
    samples = [f"NA{i:05d}" for i in range(n)]
    code = {0: "0|0", 1: "0|1", 2: "1|1"}
    lines = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for j in range(m):
        lines.append("\t".join(["22", str(20_000_000 + 5000 * j), ".", "A", "G", ".", ".", ".", "GT"] + [code[int(x)] for x in G[:, j]]))
    (d / "t.vcf").write_text("\n".join(lines) + "\n")
    sp = {"POP1": "EUR", "POP2": "AFR", "POP3": "EAS"}
    pd.DataFrame({"sample": samples, "pop": lab, "super_pop": [sp[x] for x in lab], "gender": ["male"] * n}).to_csv(d / "panel.txt", sep="\t", index=False)


def test_script_02_records_v1_and_manifest(tmp_path, plink2_bin):
    write_vcf_and_panel(tmp_path)
    cp = subprocess.run([sys.executable, str(ROOT / "scripts" / "02_prepare_1000g_pca.py"), "--vcf", str(tmp_path / "t.vcf"),
                         "--panel", str(tmp_path / "panel.txt"), "--work-dir", str(tmp_path / "w"), "--artifact-dir", str(tmp_path / "a"),
                         "--no-long-range-exclusion", "--pca-engine", "python", "--king", "--n-pcs", "5"],
                        capture_output=True, text=True, cwd=ROOT, env={**__import__("os").environ, "PLINK2": plink2_bin})
    assert cp.returncode == 0, cp.stdout[-1500:] + cp.stderr[-2500:]
    st = json.loads((tmp_path / "a" / "VALIDATION_STATUS.json").read_text())
    assert st["V1"]["knn_cv_balanced_accuracy"] > 0.9 and "V4" not in st and "V5" not in st
    man = json.loads((tmp_path / "a" / "1000g_pca.manifest.json").read_text())
    assert any("long-range" in x for x in man["inclusion_exclusion"]["deviations"])
    assert man["qc_thresholds"]["hwe"].startswith("NOT APPLIED")
    pca = pd.read_csv(tmp_path / "a" / "1000g_pca.csv.gz")
    assert {"sample_id", "PC1", "super_pop"} <= set(pca.columns) and len(pca) == 150


def test_script_02_refuses_silent_long_range_omission(tmp_path, plink2_bin):
    cp = subprocess.run([sys.executable, str(ROOT / "scripts" / "02_prepare_1000g_pca.py"), "--vcf", "x", "--panel", "y"],
                        capture_output=True, text=True, cwd=ROOT, env={**__import__("os").environ, "PLINK2": plink2_bin})
    assert cp.returncode != 0 and "long-range" in (cp.stdout + cp.stderr)
