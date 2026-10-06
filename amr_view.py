"""EvoResist-AI v0.2 artifact viewer (presentation only). Consumes compact validated artifacts from
artifacts/amr/ produced offline by scripts/amr_*.py. Nothing is downloaded or computed here.

Four evidence dimensions are always rendered in separate labelled sections and never combined:
association_evidence | predictive_contribution | biological_plausibility | causal_status."""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from variantbridge.amr import config as amr_config, reporting as amr_reporting, ui_data, validation as amr_validation
from variantbridge.amr.config import load_config, breakpoint_block

ART_AMR = Path(__file__).resolve().parent / "artifacts" / "amr"
ENV_OVERRIDE = "VARIANTBRIDGE_AMR_ARTIFACTS"       # test hook: point the viewer at another artifact directory


def _chip(s):
    return {"PASSED": "✅ PASSED", "FAILED": "❌ FAILED", "INCOMPLETE": "🟡 INCOMPLETE", "UNVERIFIED": "🟠 UNVERIFIED"}.get(s, f"⚪ {s}")


def render(tabs, root: Path | None = None):
    root = Path(os.environ[ENV_OVERRIDE]) if root is None and os.environ.get(ENV_OVERRIDE) else (root or ART_AMR)
    data = ui_data.load_all(root)
    cfg = load_config()
    sc = ui_data.scopes(data)
    synthetic = "synthetic" in sc
    # ------------------------------------------------------------------ overview
    with tabs[0]:
        st.subheader("EvoResist-AI v0.2 — E. coli × ciprofloxacin (offline results)")
        if not sc:
            st.warning("NOT RUN: no validated EvoResist artifacts found in artifacts/amr/. Run scripts/amr_01 → amr_06 offline "
                       "(see README). The app never downloads or processes genomic data.")
        if synthetic:
            st.error("SYNTHETIC DATA: these artifacts come from the software plumbing demonstration. They are NOT evidence about real resistance biology.")
        st.markdown("**Cohort architecture:** discovery = PMID 38052776 · replication = PMID 34485958 · exploratory/generalisation = remaining eligible BV-BRC data (descriptive only).")
        st.markdown("**Primary model:** interval-censored Gaussian regression of log2 MIC (censoring is explicit; raw MIC, boundary, direction and transforms are stored separately).")
        bp = breakpoint_block(cfg, "discovery")
        st.info(f"Discovery-cohort breakpoint status: **{bp['status']}**. No S/I/R analysis is shown while this is UNVERIFIED; no breakpoint is ever inferred.")
        ph = data["phenotypes"]
        if ph is not None:
            c = st.columns(3)
            for col, t in zip(c, amr_config.TIERS):
                col.metric(f"{t} isolates", int((ph["cohort_tier"] == t).sum()))
    # ------------------------------------------------------------------ provenance & status
    with tabs[1]:
        st.subheader("Validation status (real-data scope only)")
        rows = [{"id": k, "validation": amr_validation.AMR_VALIDATION_IDS[k], "status": _chip(v["status"]), "note": v.get("note", "")} for k, v in data["status"].items()]
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        st.caption("Absent ⇒ NOT RUN. Synthetic runs are never shown as validation. Passing software tests is not evidence that EvoResist-AI recovers real resistance biology.")
        a6 = data["status"]["A6"]
        if a6.get("status") == "UNVERIFIED":
            with st.expander("A6 — what was and was not found for the discovery breakpoint"):
                for k in ("source", "found", "not_found", "also_noted", "consequence"):
                    if a6.get(k):
                        st.markdown(f"**{k.replace('_', ' ').capitalize()}.** {a6[k]}")
        key, sha = amr_config.load_answer_key()
        st.markdown(f"**Pre-registered answer key** registered {key['registered_on']} · SHA-256 `{sha}`")
        for k in ("phenotypes", "association", "importance"):
            m = data.get(k + "_manifest")
            if m:
                with st.expander(f"Manifest: {k} (data_scope = {m.get('data_scope')})"):
                    st.json({kk: vv for kk, vv in m.items() if kk not in ("software_versions",)})
        for k, msg in data["messages"].items():
            st.caption(f"{k}: {msg}")
    # ------------------------------------------------------------------ QC
    with tabs[2]:
        st.subheader("Cohort assignment, duplicates and genome QC")
        ex = data["exclusions"]
        if ex is None:
            st.info("NOT RUN: exclusion table unavailable.")
        else:
            st.dataframe(ex["reason"].value_counts().rename("n").reset_index(), width="stretch", hide_index=True)
        m = data["phenotypes_manifest"]
        if m:
            st.json(m["model_specification"].get("qc", {}))
            st.markdown("**Threshold provenance (every threshold is configuration, none immutable):**")
            st.dataframe(pd.DataFrame(m["qc_thresholds"]), width="stretch", hide_index=True)
    # ------------------------------------------------------------------ population structure / censoring
    with tabs[3]:
        st.subheader("Lineages and MIC censoring")
        if ph is None:
            st.info("NOT RUN: phenotype artifact unavailable.")
        else:
            lin = ph.groupby(["cohort_tier", "mlst"], dropna=False).size().rename("n").reset_index()
            st.plotly_chart(px.bar(lin.sort_values("n", ascending=False).head(40), x="mlst", y="n", color="cohort_tier", title="Isolates per ST (top 40)"), width="stretch")
            dist = ph.groupby(["cohort_tier", "mic_operator", "mic_boundary_mg_l"]).size().rename("n").reset_index()
            dist["reported"] = dist["mic_operator"].astype(str) + dist["mic_boundary_mg_l"].astype(str)
            st.plotly_chart(px.bar(dist, x="reported", y="n", color="censoring" if "censoring" in dist else "mic_operator", facet_col="cohort_tier", title="Reported MIC (operators kept; censored values are never pooled with exact values)"), width="stretch")
    # ------------------------------------------------------------------ association
    with tabs[4]:
        st.subheader("Discovery association (genome-wide; no biological pre-filter)")
        a = data["association"]
        if a is None:
            st.info("NOT RUN: association artifact unavailable.")
        else:
            m = data["association_manifest"]
            st.caption(f"Model: {m['model_specification']['model']} · effect = log2-MIC difference per feature present · lineage fixed effects · BH-FDR over {len(a)} tested features · "
                       f"answer-key hash recorded before the run: `{m['answer_key_sha256'][:16]}…`")
            if m.get("consistency_violations"):
                st.error(f"Statistical-consistency violations: {m['consistency_violations'][:5]}")
            a = a.assign(neglog10p=-np.log10(a["p_value"].clip(lower=1e-300)))
            st.plotly_chart(px.scatter(a, x="effect", y="neglog10p", hover_name="feature_id", color=(a["q_value"] < 0.05).map({True: "q<0.05", False: "n.s."}), title="Effect vs significance"), width="stretch")
            st.dataframe(a[["feature_id", "feature_type", "n_present", "effect", "se", "ci_low", "ci_high", "p_value", "q_value", "unadjusted_p", "status"]].head(200), width="stretch", hide_index=True)
            rp = data["replication"]
            st.subheader("Replication of pre-specified discovery hits (PMID 34485958)")
            st.dataframe(rp, width="stretch", hide_index=True) if rp is not None else st.info("NOT RUN: replication artifact unavailable.")
            akv = data["answer_key_validation"]
            if akv:
                st.markdown(f"**Answer-key validation:** {_chip(akv['status'])} (data_scope = {akv['data_scope']})")
                st.json(akv["criteria"])
    # ------------------------------------------------------------------ prediction
    with tabs[5]:
        st.subheader("Predictive contribution (separate from association)")
        imp, pm = data["importance"], data["importance_manifest"]
        if pm is None:
            st.info("NOT RUN: prediction artifact unavailable.")
        else:
            res = pm["results"]
            c = st.columns(3)
            c[0].metric("Lineage-blocked extreme-contrast AUC", f"{res['lineage_blocked']['mean_enet_extreme_auc']:.3f}")
            c[1].metric("Interval-aware MAE (log2 MIC)", f"{res['lineage_blocked']['mean_interval_mae']:.3f}")
            d2r = res.get("discovery_to_replication")
            c[2].metric("Discovery → replication AUC", "NOT RUN" if not d2r else f"{d2r['enet_extreme_auc']:.3f}")
            st.caption("Extreme-contrast AUC is breakpoint-free and is NOT a clinical S/R classifier. Feature importance reflects predictive utility, not causal contribution.")
            with st.expander("Random-split result — LEAKAGE DEMONSTRATION ONLY, never performance"):
                st.write(res["random_split_leakage_demo"]); st.write(res["leakage_gap"])
            st.dataframe(imp.head(50), width="stretch", hide_index=True)
    # ------------------------------------------------------------------ evidence cards
    with tabs[6]:
        st.subheader("Evidence Cards (schema v1.1)")
        if data["cards_rejected"]:
            st.error(f"{len(data['cards_rejected'])} card(s) violated the claim restrictions and are not shown.")
        cards = data["cards"]
        if not cards:
            st.info("NOT RUN: no valid Evidence Cards (script amr_06).")
        else:
            pick = st.selectbox("Candidate", [f"{c['card_id']} · {c['candidate']['feature_id']}" for c in cards])
            c = cards[[f"{x['card_id']} · {x['candidate']['feature_id']}" for x in cards].index(pick)]
            if c["generated"]["data_scope"] == "synthetic":
                st.error("SYNTHETIC DATA: method demonstration only.")
            cols = st.columns(2)
            with cols[0].container(border=True):
                st.markdown("### 1 · Association evidence")
                st.json({**c["association_evidence"], **{"survives_lineage_adjustment": c["population_adjusted_evidence"]["survives_adjustment"]}})
            with cols[1].container(border=True):
                st.markdown("### 2 · Predictive contribution")
                st.json(c["predictive_contribution"] if c["predictive_contribution"] else "not used in the predictive baseline")
            cols = st.columns(2)
            with cols[0].container(border=True):
                st.markdown("### 3 · Biological plausibility (interpretation only)")
                st.json(c["biological_plausibility"])
            with cols[1].container(border=True):
                st.markdown("### 4 · Causal status")
                st.json(c["causal_status"])
            st.markdown("**Claims currently justified**"); [st.success(x) for x in c["claim_currently_justified"]]
            st.markdown("**Claims NOT justified**"); [st.error(x) for x in c["claim_not_justified"]]
            with st.expander("Uncertainty · evidence still missing · suggested validation"):
                st.write(c["uncertainty"]); st.write(c["evidence_still_missing"]); st.json(c["suggested_experimental_validation"])
    # ------------------------------------------------------------------ methods
    with tabs[7]:
        st.subheader("Methods and limits")
        st.code(amr_reporting.methods_report(cfg))
        st.markdown("""
- Censored MIC is modelled as censored; a naive linear model on `<=`/`>` boundary values is never fitted.
- Breakpoint-dependent analyses (S/I/R labels, binary logistic) are NOT RUN while the breakpoint is UNVERIFIED.
- The discovery sampling frame is enriched for 3GC-resistant isolates; associations are internal to the cohort.
- Exploratory-tier data are descriptive only and never enter headline association or predictive metrics.
- Linked features cannot be separated; association is not causation; prediction is not mechanism.
- A passing software test or a recovered *planted* determinant on synthetic data is not biological validation.""")
