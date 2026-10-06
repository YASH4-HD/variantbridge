"""VariantBridge v0.2 - Streamlit UI (presentation only).

All scientific computation lives in src/variantbridge and is independent of Streamlit. The app
NEVER downloads or processes raw genomic datasets: real-data results come from compact validated
artifacts produced offline by scripts/ (artifacts/). Anything not yet executed shows NOT RUN.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from variantbridge import (association, artifacts, capability, demo, evidence, eqtl, io as vio, plots, population_structure as ps,
                           prediction, prs, qc, reporting)
from variantbridge.stats_utils import genomic_inflation_lambda

ART = ROOT / "artifacts"
COHORT, AMR = "CohortOmics - Human cohorts", "EvoResist-AI - Pathogen/AMR"

st.set_page_config(page_title="VariantBridge", page_icon="🧬", layout="wide")
st.title("🧬 VariantBridge v0.2")
st.caption("CohortOmics + EvoResist-AI • Interpretable genomics-to-evidence workbench")
st.info("Research/education prototype. Outputs are exploratory: not clinical, diagnostic, therapeutic or causal claims. "
        "Association is not causation; prediction is not mechanism.")


def badge(tier: str) -> None:
    st.markdown(f"**Capability tier:** `{tier}`")


def status_chip(s: str) -> str:
    return {"PASSED": "✅ PASSED", "FAILED": "❌ FAILED"}.get(s, f"⚪ {s}")


@st.cache_data(show_spinner=False)
def cached_demo(branch: str):
    return demo.human_demo() if branch == COHORT else demo.amr_demo()


@st.cache_data(show_spinner=False)
def cached_scan(X: pd.DataFrame, y: pd.Series, cov, ids: tuple):
    return association.scan(X.to_numpy(), y.to_numpy(), cov, variant_ids=list(ids))


@st.cache_data(show_spinner="Training held-out Random Forest (v0.1 method, unchanged)...")
def cached_rf(X: pd.DataFrame, y: pd.Series):
    return prediction.holdout_random_forest(X, y)


def learning(title: str, biology: str, statistics: str, implementation: str, validation: str, failure: str, interview: str):
    with st.expander(f"Learning mode: {title}"):
        for h, t in (("Biology", biology), ("Statistics", statistics), ("Implementation", implementation),
                     ("Validation", validation), ("Common failure", failure), ("Interview question", interview)):
            st.markdown(f"**{h}.** {t}")


# ------------------------------------------------------------------------------------ sidebar
with st.sidebar:
    branch = st.radio("Branch", [COHORT, AMR])
    options = ["Demo data (SIMULATED)", "Upload CSV"] + (["Validated artifacts (real data)"] if branch == COHORT else ["EvoResist v0.2 artifacts (offline results)"])
    source = st.radio("Data", options)
    gf = pf = mf_ = None
    if source == "Upload CSV":
        gf = st.file_uploader("Genotype/features CSV", type="csv")
        pf = st.file_uploader("Phenotype CSV", type="csv")
        mf_ = st.file_uploader("Metadata CSV (optional)", type="csv")
        st.caption("Uploads: sample_id + numeric features; phenotype file: sample_id + phenotype. Missing values are kept, not silently imputed.")
    st.caption("Streamlit never downloads or processes raw genomic data. Real results come from offline scripts → artifacts/.")

is_cohort = branch == COHORT

# --------------------------------------------------------------------- tabs available for all
tab_names = ["Overview", "Provenance & Status", "QC", "Population Structure", "Association", "Scores / Prediction", "Evidence Cards", "Methods & Limits"]
tabs = st.tabs(tab_names)

if branch == AMR and source.startswith("EvoResist v0.2"):
    import amr_view          # EvoResist-AI v0.2 artifact viewer (presentation only); CohortOmics code path untouched
    amr_view.render(tabs)
    st.stop()

# ------------------------------------------------------------------------------ data loading
ge = ph = meta = None
data_tier = capability.SIMULATED
phen_col = "phenotype"
M = X = y = None
features: list = []
if source != "Validated artifacts (real data)":
    if source.startswith("Demo"):
        ge, ph, meta = cached_demo(branch)
        data_tier = capability.SIMULATED
    else:
        data_tier = "USER-UPLOADED DATA (tier depends on your data; app makes no validation claim)"
        if not gf or not pf:
            with tabs[0]:
                st.warning("Upload genotype/features and phenotype CSVs to continue.")
            st.stop()
        ge, ph = pd.read_csv(gf), pd.read_csv(pf)
        meta = pd.read_csv(mf_) if mf_ else None
    if "sample_id" not in ge or "sample_id" not in ph:
        st.error("Both files require a sample_id column.")
        st.stop()
    with st.sidebar:
        cands = [c for c in ph if c != "sample_id"]
        phen_col = "phenotype" if "phenotype" in ph else st.selectbox("Phenotype column", cands)
    try:
        M, X, y, features = vio.load_csv_inputs(ge, ph, phen_col)
    except ValueError as e:
        st.error(str(e))
        st.stop()

# ----------------------------------------------------------------------------------- Overview
with tabs[0]:
    st.subheader("What this app does")
    st.write("VariantBridge separates **association, prediction, biological interpretation and causal validation** instead of treating a "
             "model output as biological proof. Every analysis carries a capability tier; every validation shows its real status.")
    if X is not None:
        a, b, c, d = st.columns(4)
        a.metric("Samples", len(M)); b.metric("Features", len(features)); c.metric("Mode", "Human" if is_cohort else "AMR")
        d.metric("Data", "SIMULATED" if source.startswith("Demo") else "Uploaded")
        badge(data_tier)
        with st.expander("Preview"):
            st.dataframe(ge.head(), width="stretch"); st.dataframe(ph.head(), width="stretch")
    else:
        st.write("Showing validated artifacts only (offline-prepared real data).")
    st.markdown("#### Capability tiers")
    stt = capability.load_status(ART)
    st.dataframe(pd.DataFrame([{"Analysis": k, "Tier": t, "Validation": (f"{v} - {status_chip(stt[v]['status'])}" if v else "")}
                               for k, (t, v) in capability.CAPABILITIES.items()]), hide_index=True, width="stretch")

# ----------------------------------------------------------------------- Provenance & Status
with tabs[1]:
    st.subheader("Validation status (never fabricated)")
    st.write("Statuses are read from `artifacts/VALIDATION_STATUS.json`, written only by scripts after real execution. Absent = **NOT RUN**.")
    stt = capability.load_status(ART)
    for vid, desc in capability.VALIDATION_IDS.items():
        rec = stt[vid]
        st.markdown(f"**{vid}** - {desc}: {status_chip(rec['status'])}" + (f"  _(executed {rec['executed_at']})_" if "executed_at" in rec else ""))
        if vid == "V5" and rec["status"] != capability.NOT_RUN:
            st.caption(f"Scope: {rec.get('scope')} - {rec.get('dataset')}. {rec.get('remaining', '')}")
        if vid == "V6" and rec["status"] != capability.NOT_RUN:
            be = rec["biological_evidence"]
            st.caption(f"Implementation: {rec['implementation']['status']} (software). Biological evidence reported separately: "
                       f"mean r = {be['observed_mean_r']:.3f}, 95% bootstrap CI [{be['bootstrap_ci95_mean_r'][0]:.3f}, {be['bootstrap_ci95_mean_r'][1]:.3f}], "
                       f"empirical p = {be['empirical_p_one_sided']:.4f} ({be['n_permutations']} permutations). {be['interpretation']}")
    st.markdown("#### Dataset manifests")
    found = False
    for mp in sorted(ART.glob("*.manifest.json")):
        found = True
        with st.expander(mp.name):
            st.json(artifacts.read_manifest(mp) if hasattr(artifacts, "read_manifest") else __import__("json").loads(mp.read_text()))
    if not found:
        st.info("No dataset manifests yet: offline preparation has NOT RUN. See README for scripts/01-08.")

# ------------------------------------------------------------------------------------ QC
with tabs[2]:
    if X is None:
        st.info("QC of real data is performed offline (PLINK2 + scripts/02-03); see manifests under Provenance & Status.")
    else:
        badge(data_tier)
        G = X.to_numpy(float)
        sm_, vm_ = qc.sample_metrics(G), qc.variant_metrics(G)
        thr = qc.QCThresholds()
        st.dataframe(pd.DataFrame({"Metric": ["Samples", "Features", "Mean sample missingness", "Mean feature missingness", "Median MAF/minor frequency"],
                                   "Value": [len(M), len(features), 1 - sm_["call_rate"].mean(), vm_["missing_rate"].mean(), float(np.nanmedian(vm_["maf"]))]}),
                     hide_index=True, width="stretch")
        a, b = st.columns(2)
        a.plotly_chart(px.histogram(1 - sm_["call_rate"], nbins=30, title="Sample missingness"), width="stretch")
        b.plotly_chart(px.histogram(vm_["missing_rate"], nbins=30, title="Feature missingness"), width="stretch")
        st.plotly_chart(px.histogram(vm_["maf"].dropna(), nbins=30, title="Minor allele / feature frequency"), width="stretch")
        st.caption("Reference thresholds (project design choices within conventional ranges) are displayed, not applied automatically: "
                   f"{thr.provenance()}. HWE is never applied to pooled multi-ancestry data (Wahlund effect).")
        learning("QC", "Poor-quality samples and variants create false associations.", "Call rate, heterozygosity (±3 SD), MAF and missingness filters; exact HWE only within populations/controls.",
                 "qc.py computes metrics; thresholds are explicit arguments.", "Known-answer sex check (V2) and unit tests incl. an exact HWE brute-force check.",
                 "A single pooled HWE filter removes real variants under population structure.", "Why not filter HWE on the combined 1000 Genomes data?")

# ------------------------------------------------------------------------ Population structure
with tabs[3]:
    if source == "Validated artifacts (real data)":
        df, man, msg = artifacts.load_table_artifact("1000g_pca.csv.gz", ART)
        if df is None:
            st.warning(f"NOT RUN: {msg}")
        else:
            badge(capability.STATISTICALLY_DEFENSIBLE)
            st.plotly_chart(px.scatter(df, x="PC1", y="PC2", color="super_pop", hover_name="sample_id", title="1000 Genomes PCA by superpopulation"), width="stretch")
            v1 = capability.load_status(ART)["V1"]
            st.write(f"V1: {status_chip(v1['status'])}")
            if v1["status"] != capability.NOT_RUN:
                st.write(f"kNN cross-validated balanced accuracy of superpopulation from the first {v1['n_pcs_used']} PCs: "
                         f"**{v1['knn_cv_balanced_accuracy']:.3f}**; silhouette {v1['silhouette']:.3f}.")
                st.json(v1["eta_squared_by_pc"])
    else:
        badge(data_tier)
        if is_cohort:
            pcs_, ratio = ps.genotype_pca(X.to_numpy(float), 5)[:2]
        else:
            pcs_, ratio = ps.feature_pca_v01(X.to_numpy(float), 5)
        P = ps.pcs_frame(pcs_, M["sample_id"])
        st.plotly_chart(plots.pca_figure(P, y, ratio, "PCA"), width="stretch")
        st.dataframe(pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(ratio))], "Explained variance": ratio}), hide_index=True)
        st.caption("Colour = phenotype. For real data, V1 checks that leading PCs recover known superpopulation structure.")
        learning("Population structure", "Ancestry/lineage differences in allele frequency.", "PCA of normalised genotypes; leading PCs as covariates.",
                 "population_structure.py", "V1: leading PCs vs known 1000 Genomes superpopulation labels (kNN balanced accuracy, eta², silhouette).",
                 "Structure unrelated to the phenotype can still confound if phenotype prevalence differs by ancestry.", "Why include PCs as covariates?")

# ------------------------------------------------------------------------------------ Association
assoc = None
with tabs[4]:
    if source == "Validated artifacts (real data)":
        st.subheader("GEUVADIS cis-eQTL scan (EUR / YRI analysed separately)")
        for pop in ("EUR", "YRI"):
            df, man, msg = artifacts.load_table_artifact(f"geuvadis_eqtl_{pop}.csv.gz", ART)
            if df is None:
                st.warning(f"{pop}: NOT RUN - {msg}")
                continue
            badge(capability.STATISTICALLY_DEFENSIBLE)
            st.write(f"**{pop}**: {len(df)} genes tested; eGenes (permutation BH q<0.05): {int((df['q_bh'] < 0.05).sum()) if 'q_bh' in df else 'n/a'}")
            st.dataframe(df.sort_values("p_nominal").head(50), hide_index=True, width="stretch")
            st.download_button(f"Download {pop} results", df.to_csv(index=False), f"geuvadis_eqtl_{pop}.csv", key=f"dl{pop}")
        v4 = capability.load_status(ART)["V4"]
        st.write(f"**V4 answer-key replication:** {status_chip(v4['status'])}" + ("" if v4["status"] != capability.NOT_RUN else " (not executed: no replication claim is made)"))
        st.caption("No 5e-8 GWAS threshold is used here; calibration is per-gene permutation + BH.")
    else:
        badge(data_tier)
        npc = st.slider("PC covariates", 0, 5, 2)
        if npc:
            if is_cohort:
                pcs_ = ps.genotype_pca(X.to_numpy(float), npc)[0]
            else:
                pcs_ = ps.feature_pca_v01(X.to_numpy(float), npc)[0]
            cov = pd.DataFrame(pcs_, columns=[f"PC{i + 1}" for i in range(pcs_.shape[1])])
        else:
            cov = None
        res = cached_scan(X, y, cov, tuple(features))
        assoc = association.annotate_with_metadata(res, meta)
        ok = assoc[assoc["status"] == "ok"]
        a, b = st.columns(2)
        a.plotly_chart(plots.manhattan_figure(ok), width="stretch")
        b.plotly_chart(plots.qq_figure(ok["p"]), width="stretch")
        lam = genomic_inflation_lambda(ok["p"])
        st.metric("Genomic inflation λ (diagnostic only, not a correction)", f"{lam:.3f}")
        bad = assoc[assoc["status"] != "ok"]
        if len(bad):
            st.warning(f"{len(bad)} features not testable (status shown in table): degenerate/separation/non-convergence.")
        viol = association.check_statistical_consistency(res)
        st.caption(f"Model: {res['model'].iloc[0]}; effect, SE, CI and p-value come from the SAME fitted model "
                   f"(internal consistency violations: {len(viol)}). Benchmark vs PLINK2 (V5): {status_chip(capability.load_status(ART)['V5']['status'])}")
        st.dataframe(assoc.head(30), hide_index=True, width="stretch")
        st.download_button("Download associations", assoc.to_csv(index=False), "variantbridge_associations.csv", "text/csv")
        st.warning("Exploratory screening on this dataset. Not a publication-grade GWAS: no relatedness model, no LD handling. Reproduce with PLINK2/REGENIE/SAIGE (human) or pyseer (microbial).")
        learning("Association", "Does genotype/feature differ with the phenotype?", "OLS (t-test) or logistic (Wald), covariate-adjusted; beta or log-odds; p-values multiple-testing corrected.",
                 "association.py: vectorised OLS (Frisch-Waugh-Lovell) and IRLS logistic; one fitted model yields effect, SE, CI and p.",
                 "V5: matches plink2 --glm on identical data/model; unit tests vs statsmodels.", "Unadjusted p with adjusted effect (the v0.1 flaw), population-structure confounding.",
                 "Why must beta, SE, CI and p come from the same model?")

# ------------------------------------------------------------------------------ Scores / Prediction
with tabs[5]:
    if is_cohort:
        st.subheader("Genetic Expression Score (predicted expression) - NOT a PRS")
        badge(capability.STATISTICALLY_DEFENSIBLE)
        v6 = capability.load_status(ART)["V6"]
        if v6["status"] == capability.NOT_RUN:
            st.warning("NOT RUN: the eQTLGen → GEUVADIS expression-score workflow has not been executed (scripts/08). No performance is claimed.")
        else:
            be = v6["biological_evidence"]
            st.write(f"Software implementation checks: **{v6['implementation']['status']}**.")
            st.write(f"Biological predictive evidence (reported separately): mean per-gene r = **{be['observed_mean_r']:.3f}** "
                     f"(95% bootstrap CI {be['bootstrap_ci95_mean_r'][0]:.3f} to {be['bootstrap_ci95_mean_r'][1]:.3f}); permutation-null 95% range "
                     f"{be['null_ci95'][0]:.3f} to {be['null_ci95'][1]:.3f}; one-sided empirical p = {be['empirical_p_one_sided']:.4f} "
                     f"({be['n_permutations']} permutations, {be['n_genes_scored']} genes, {be['n_samples']} samples).")
            st.info(be["interpretation"]); st.caption(be["caveat"])
        st.divider()
        st.subheader("Polygenic score method demonstration - SIMULATED")
        badge(capability.SIMULATED)
        st.caption(prs.DISEASE_PRS_LIMITATION)
        Gd, yd, Gt, yt, beta = demo.simulate_discovery_target(seed=3)
        prs.assert_disjoint_cohorts(Gd.index, Gt.index)
        thr_p = st.select_slider("Discovery p-value threshold (fixed in advance; never tune on the target)", [1e-2, 1e-3, 1e-4, 1e-5], 1e-4)
        disc = association.scan(Gd.to_numpy(), yd.to_numpy(), None, list(Gd.columns))
        w = prs.select_weights(disc, thr_p)
        if len(w):
            sc = prs.polygenic_score(Gt[w["variant"]].to_numpy(), w["effect"].to_numpy())
            ev = prs.evaluate_polygenic_score(sc, yt.to_numpy(), n_permutations=300, n_bootstrap=200)
            st.write(f"{len(w)} variants selected in the discovery cohort (n={len(Gd)}); evaluated in the disjoint target cohort (n={len(Gt)}). "
                     f"{ev['metric']} = **{ev['observed']:.3f}** (bootstrap 95% CI {ev['bootstrap_ci95'][0]:.3f}-{ev['bootstrap_ci95'][1]:.3f}); "
                     f"permutation empirical p = {ev['empirical_p_one_sided']:.4f}. SIMULATED: demonstrates the leakage-controlled workflow only.")
        else:
            st.write("No discovery variants pass this threshold.")
        learning("Expression score / PRS", "Genotype-based prediction of a molecular or complex trait.", "Weighted sum of allele dosages from independent discovery weights, evaluated in a disjoint target against a permutation null.",
                 "expression_score.py (expression) and prs.py (complex traits); clumping + thresholding.", "V6: implementation checks and biological performance reported separately.",
                 "Discovery/target overlap, tuning on the target, ancestry mismatch.", "What would make this score leak information?")
    else:
        st.subheader("Held-out interpretable AMR prediction (exploratory)")
        badge(data_tier)
        out = cached_rf(X, y)
        a, b, c = st.columns(3)
        a.metric("ROC-AUC", f"{out['roc_auc']:.3f}"); b.metric("Accuracy", f"{out['accuracy']:.3f}"); c.metric("Balanced accuracy", f"{out['balanced_accuracy']:.3f}")
        fi = out["importance"]
        st.plotly_chart(px.bar(fi.head(20).sort_values("importance"), x="importance", y="feature", orientation="h", error_x="sd"), width="stretch")
        st.dataframe(out["report"], width="stretch")
        st.warning("Predictive importance is not causal evidence; lineage, linkage, sampling and phenotype quality can generate apparent predictors.")

# ------------------------------------------------------------------------------- Evidence cards
with tabs[6]:
    if assoc is None:
        st.info("Evidence Cards are generated from an association table of the selected demo/uploaded data.")
    else:
        ok = assoc[assoc["status"] == "ok"].reset_index(drop=True)
        chosen = st.selectbox("Candidate", ok.head(25)["variant"])
        row = ok[ok["variant"] == chosen].iloc[0]
        ann = str(row["known_amr"]) if "known_amr" in row.index and pd.notna(row.get("known_amr")) else None
        card = evidence.make_card(row, phen_col, "PC covariates as selected (Association tab)", data_tier, known_annotation=ann)
        a, b = st.columns(2)
        a.markdown(f"### {card.candidate}")
        a.write(f"**Effect:** {card.effect:.3f} ({card.effect_scale}); **95% CI** [{card.ci_low:.3f}, {card.ci_high:.3f}]")
        a.write(f"**p:** {card.p:.3g}"); a.write(f"**Multiple testing:** {card.multiple_testing}"); a.write(f"**Known annotation:** {card.known_database_evidence}")
        b.write(f"**Claim currently justified:** {card.claim_justified}")
        b.write("**Claims NOT justified:** " + "; ".join(card.claims_not_justified))
        b.write("**Evidence missing:** " + "; ".join(card.evidence_missing))
        st.download_button("Download evidence card", evidence.render_text(card), f"{chosen.replace(':', '_')}_evidence.txt")

# ------------------------------------------------------------------------------ Methods & limits
with tabs[7]:
    st.download_button("Download methods & validation status report", reporting.methods_report(ART), "variantbridge_methods_status.md")
    st.markdown("""### Scope
**CohortOmics:** provenance → QC → relatedness → PCA → association (cis-eQTL on GEUVADIS; GWAS-style teaching module on simulated data) → Genetic Expression Score → evidence reporting.
**EvoResist-AI:** pathogen feature QC → lineage/population structure → association → held-out interpretable prediction → evidence cards (methodology unchanged from v0.1).

### Limitations
- Demo datasets are **simulated** and are never validation evidence.
- Real-data validations (V1–V4, V6) show **NOT RUN** until the offline scripts are executed on the downloaded datasets.
- No relatedness/LD-aware mixed models in the app; not a replacement for PLINK2, SAIGE, REGENIE, GEMMA, pyseer.
- ML feature importance is not mechanistic or causal evidence.
- A rigorous disease-trait PRS needs a controlled-access target cohort; none is claimed here.
- No clinical or therapeutic decisions should be based on this application.
""")
