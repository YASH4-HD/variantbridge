"""SIMULATED GEUVADIS-like inputs in the format ASSUMED by scripts/03 (real-file format UNVERIFIED).

Used only for pipeline-plumbing tests. Nothing generated here is validation evidence."""
from pathlib import Path

import numpy as np
import pandas as pd


def write_geuvadis_like(outdir: Path, seed=0, n_eur=100, n_yri=60, n_genes=24, n_var=70, planted=(0, 1, 2, 3, 4, 5), effect=1.1):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    r = np.random.default_rng(seed)
    pops = ["CEU"] * (n_eur // 4) + ["FIN"] * (n_eur // 4) + ["GBR"] * (n_eur // 4) + ["TSI"] * (n_eur // 4) + ["YRI"] * n_yri
    samples = [f"HG{i:05d}" for i in range(len(pops))]
    group = np.array(["EUR" if p != "YRI" else "YRI" for p in pops])
    # panel
    pd.DataFrame({"sample": samples, "pop": pops, "super_pop": np.where(group == "EUR", "EUR", "AFR"),
                  "gender": r.choice(["male", "female"], len(samples))}).to_csv(outdir / "panel.txt", sep="\t", index=False)
    pd.DataFrame({"s": samples, "imputed": (r.random(len(samples)) < 0.1).astype(int)}).to_csv(outdir / "imputed.tsv", sep="\t", header=False, index=False)
    n = len(samples)
    lines = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    code = {0: "0|0", 1: "0|1", 2: "1|1"}
    expr_rows, G_by_gene, conv = [], {}, []
    # population stratification in allele frequencies
    for g in range(n_genes):
        tss = 20_000_000 + g * 4_000_000
        base = r.uniform(0.15, 0.45, n_var)
        shift = r.uniform(-0.1, 0.1, n_var)
        freq = np.clip(np.where(group[:, None] == "EUR", base, base + shift), 0.02, 0.98)
        G = r.binomial(2, freq).astype(int)
        for j in range(n_var):
            pos = tss - 300_000 + j * 8_000
            old = f"chr22_{pos}_A_G"; rs = f"rs{g}_{j}"
            conv.append((old, rs))
            lines.append("\t".join(["22", str(pos), old, "A", "G", ".", ".", ".", "GT"] + [code[x] for x in G[:, j]]))
        G_by_gene[g] = G
        y = r.normal(size=n)
        if g in planted:
            y = y + effect * (G[:, 10] - G[:, 10].mean())
        expr_rows.append((f"ENSG{g:05d}.4", "22", tss, y))
    # VCF must be position-sorted overall
    head, body = lines[:3], sorted(lines[3:], key=lambda l: int(l.split("\t")[1]))
    (outdir / "geu.chr22.vcf").write_text("\n".join(head + body) + "\n")
    pd.DataFrame(conv).to_csv(outdir / "idconvert.txt", sep="\t", header=False, index=False)
    ex = pd.concat([pd.DataFrame({"Gene_Symbol": [e[0] for e in expr_rows], "Chr": [e[1] for e in expr_rows], "Coord": [e[2] for e in expr_rows]}),
                    pd.DataFrame({s: [e[3][i] for e in expr_rows] for i, s in enumerate(samples)})], axis=1)
    ex.to_csv(outdir / "expression.txt", sep="\t", index=False)
    return dict(samples=samples, group=group, G_by_gene=G_by_gene, planted=planted)
