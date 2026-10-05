import math, numpy as np, pandas as pd, streamlit as st, plotly.express as px, plotly.graph_objects as go
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, classification_report
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance

st.set_page_config(page_title="VariantBridge",page_icon="🧬",layout="wide")
st.title("🧬 VariantBridge")
st.caption("CohortOmics + EvoResist-AI • Interpretable genomics-to-evidence workbench")
st.info("Research/education prototype. Outputs are exploratory; they are not clinical, diagnostic, therapeutic, or causal claims.")

def human_demo(seed=7,n=320,m=180):
    r=np.random.default_rng(seed); ids=[f"H{i:04d}" for i in range(n)]
    ch=r.integers(1,23,m); pos=r.integers(1_000_000,180_000_000,m)
    v=[f"chr{a}:{b}:A:G" for a,b in zip(ch,pos)]; maf=r.uniform(.05,.48,m); anc=r.integers(0,2,n)
    G=np.zeros((n,m),float)
    for j,p in enumerate(maf):
        pp=np.clip(p+(anc-.5)*(.18 if j<25 else .02),.01,.95); G[:,j]=r.binomial(2,pp)
    G[r.random(G.shape)<.008]=np.nan; F=np.where(np.isnan(G),np.nanmean(G,0),G)
    age=r.normal(48,13,n).clip(18,85); sex=r.integers(0,2,n)
    logit=-2.1+.035*(age-45)+.25*sex+.65*F[:,5]-.55*F[:,17]+.45*F[:,49]
    y=r.binomial(1,1/(1+np.exp(-logit)))
    ge=pd.DataFrame(G,columns=v); ge.insert(0,"sample_id",ids)
    ph=pd.DataFrame({"sample_id":ids,"phenotype":y,"age":age.round(1),"sex":sex})
    me=pd.DataFrame({"variant":v,"chrom":ch,"pos":pos}); return ge,ph,me

def amr_demo(seed=11,n=300,m=140):
    r=np.random.default_rng(seed); ids=[f"ISO{i:04d}" for i in range(n)]; lin=r.integers(0,3,n)
    g=[f"gene_{i:03d}" for i in range(m)]; X=np.zeros((n,m))
    for j,p in enumerate(r.uniform(.05,.5,m)):
        X[:,j]=r.binomial(1,np.clip(p+(lin-1)*r.uniform(-.10,.10),.01,.95))
    causal=[2,12,37,75]; z=-1.7+1.15*X[:,2]+.9*X[:,12]-.75*X[:,37]+.7*X[:,75]+.3*(lin==2)
    y=r.binomial(1,1/(1+np.exp(-z))); ge=pd.DataFrame(X,columns=g); ge.insert(0,"sample_id",ids)
    ph=pd.DataFrame({"sample_id":ids,"phenotype":y,"lineage":lin})
    me=pd.DataFrame({"variant":g,"chrom":"chromosome","pos":np.arange(1,m+1)*10000,
                     "known_amr":["Known" if i in causal[:2] else "Uncertain" for i in range(m)]})
    return ge,ph,me

def prep(ge,ph,pcol):
    m=ge.merge(ph,on="sample_id"); f=[c for c in ge if c!="sample_id"]
    X=m[f].apply(pd.to_numeric,errors="coerce"); X=X.fillna(X.mean()).fillna(0)
    y=pd.to_numeric(m[pcol],errors="coerce"); k=y.notna()
    return m.loc[k].reset_index(drop=True),X.loc[k].reset_index(drop=True),y.loc[k].reset_index(drop=True),f

def pcs(X,ids,n=5):
    p=PCA(n_components=min(n,X.shape[1],X.shape[0]-1)); a=p.fit_transform(StandardScaler().fit_transform(X))
    d=pd.DataFrame(a,columns=[f"PC{i+1}" for i in range(a.shape[1])]); d.insert(0,"sample_id",list(ids))
    return d,p.explained_variance_ratio_

def scan(X,y,meta,cov=None):
    binary=set(y.unique()).issubset({0,1}); rows=[]
    C=pd.DataFrame(index=X.index) if cov is None else cov.reset_index(drop=True)
    for col in X:
        x=X[col].astype(float)
        if x.nunique()<2: continue
        try:
            Z=pd.concat([x.rename("variant"),C],axis=1)
            if binary:
                mdl=LogisticRegression(max_iter=1200).fit(Z,y); eff=float(mdl.coef_[0][0]); _,p=stats.pointbiserialr(y,x)
            else:
                mdl=LinearRegression().fit(Z,y); eff=float(mdl.coef_[0]); _,p=stats.pearsonr(x,y)
            rows.append((col,eff,float(p)))
        except: pass
    d=pd.DataFrame(rows,columns=["variant","effect","p"])
    if d.empty:return d
    d["neglog10p"]=-np.log10(d.p.clip(1e-300)); d["bonferroni"]=np.minimum(d.p*len(d),1)
    if meta is not None and "variant" in meta:d=d.merge(meta,on="variant",how="left")
    if "chrom" not in d:d["chrom"]="NA"
    if "pos" not in d:d["pos"]=np.arange(len(d))
    return d.sort_values("p").reset_index(drop=True)

def qq(d):
    p=np.sort(d.p.clip(1e-300,1)); e=-np.log10((np.arange(1,len(p)+1)-.5)/len(p)); o=-np.log10(p)
    f=go.Figure([go.Scatter(x=e,y=o,mode="markers",name="Observed")]); L=max(e.max(),o.max())
    f.add_trace(go.Scatter(x=[0,L],y=[0,L],mode="lines",name="Expected")); f.update_layout(title="QQ plot",xaxis_title="Expected -log10(p)",yaxis_title="Observed -log10(p)"); return f

with st.sidebar:
    branch=st.radio("Branch",["CohortOmics — Human cohorts","EvoResist-AI — Pathogen/AMR"])
    source=st.radio("Data",["Demo data","Upload CSV"])
    if source=="Upload CSV":
        gf=st.file_uploader("Genotype/features CSV",type="csv"); pf=st.file_uploader("Phenotype CSV",type="csv"); mf=st.file_uploader("Metadata CSV (optional)",type="csv")
    st.caption("Uploads: sample_id + numeric features; phenotype file: sample_id + phenotype.")

if source=="Demo data": ge,ph,meta=human_demo() if branch.startswith("Cohort") else amr_demo()
else:
    if not gf or not pf: st.warning("Upload genotype/features and phenotype CSVs."); st.stop()
    ge=pd.read_csv(gf); ph=pd.read_csv(pf); meta=pd.read_csv(mf) if mf else None
if "sample_id" not in ge or "sample_id" not in ph: st.error("Both files require sample_id."); st.stop()
pcol="phenotype" if "phenotype" in ph else st.selectbox("Phenotype column",[c for c in ph if c!="sample_id"])
M,X,y,features=prep(ge,ph,pcol)

tabs=st.tabs(["Overview","QC","Population Structure","Association","Risk / Prediction","Evidence Cards","Methods & Limits"])
with tabs[0]:
    a,b,c,d=st.columns(4); a.metric("Samples",len(ge)); b.metric("Features",len(features)); c.metric("Phenotypes",len(ph)); d.metric("Mode","Human" if branch.startswith("Cohort") else "AMR")
    st.write("VariantBridge separates association, prediction, biological interpretation and causal validation instead of treating a model output as biological proof.")
    with st.expander("Preview"): st.dataframe(ge.head(),use_container_width=True); st.dataframe(ph.head(),use_container_width=True)

with tabs[1]:
    raw=ge[features].apply(pd.to_numeric,errors="coerce"); sm=raw.isna().mean(1); vm=raw.isna().mean(0)
    freq=raw.mean(0)/(2 if np.nanmax(raw.values)>1 else 1); maf=np.minimum(freq,1-freq)
    st.dataframe(pd.DataFrame({"Metric":["Samples","Variants/features","Mean sample missingness","Mean feature missingness","Median minor frequency"],
                               "Value":[len(ge),len(features),sm.mean(),vm.mean(),np.nanmedian(maf)]}),hide_index=True,use_container_width=True)
    a,b=st.columns(2)
    with a: st.plotly_chart(px.histogram(sm,nbins=30,title="Sample missingness"),use_container_width=True)
    with b: st.plotly_chart(px.histogram(vm,nbins=30,title="Feature missingness"),use_container_width=True)
    st.plotly_chart(px.histogram(maf,nbins=30,title="Minor allele / feature frequency"),use_container_width=True)
    st.caption("Thresholds are not hard-coded because appropriate QC depends on organism, assay, cohort and downstream model.")

with tabs[2]:
    P,ev=pcs(X,M.sample_id,5); P["phenotype"]=y.astype(str)
    st.plotly_chart(px.scatter(P,x="PC1",y="PC2",color="phenotype",hover_name="sample_id",title=f"PCA • PC1 {ev[0]*100:.1f}% • PC2 {ev[1]*100:.1f}%"),use_container_width=True)
    st.dataframe(pd.DataFrame({"PC":[f"PC{i+1}" for i in range(len(ev))],"Explained variance":ev}),hide_index=True)

with tabs[3]:
    npc=st.slider("PC covariates",0,5,2); P,_=pcs(X,M.sample_id,max(npc,1)); cov=P[[f"PC{i+1}" for i in range(npc)]] if npc else None
    assoc=scan(X,y,meta,cov)
    if assoc.empty: st.warning("No testable features.")
    else:
        a,b=st.columns(2)
        with a:
            st.plotly_chart(px.scatter(assoc,x="pos",y="neglog10p",color=assoc.chrom.astype(str),hover_name="variant",title="Association overview"),use_container_width=True)
        with b: st.plotly_chart(qq(assoc),use_container_width=True)
        st.dataframe(assoc.head(30),hide_index=True,use_container_width=True)
        st.download_button("Download associations",assoc.to_csv(index=False),"variantbridge_associations.csv","text/csv")
        st.warning("v1 p-values are exploratory screening statistics. Publication-grade work should use validated cohort/organism-specific association software and models.")

with tabs[4]:
    if branch.startswith("Cohort"):
        st.subheader("Polygenic-style score demonstration")
        n=st.slider("Top variants",5,min(50,len(assoc)),min(20,len(assoc)))
        top=assoc.head(n); fs=[v for v in top.variant if v in X]; w=top.set_index("variant").loc[fs,"effect"]; score=X[fs].dot(w)
        out=pd.DataFrame({"score":score,"phenotype":y.astype(str)})
        st.plotly_chart(px.histogram(out,x="score",color="phenotype",barmode="overlay"),use_container_width=True)
        st.warning("This is not a clinically validated PRS. Real PRS requires independent effect estimates, ancestry-aware validation, calibration and an independent target cohort.")
    else:
        st.subheader("Held-out interpretable AMR prediction")
        Xtr,Xte,ytr,yte=train_test_split(X,y,test_size=.25,stratify=y,random_state=42)
        model=RandomForestClassifier(n_estimators=300,class_weight="balanced",random_state=42,n_jobs=-1).fit(Xtr,ytr)
        pr=model.predict(Xte); pb=model.predict_proba(Xte)[:,1]
        a,b,c=st.columns(3); a.metric("ROC-AUC",f"{roc_auc_score(yte,pb):.3f}"); b.metric("Accuracy",f"{accuracy_score(yte,pr):.3f}"); c.metric("Balanced accuracy",f"{balanced_accuracy_score(yte,pr):.3f}")
        pi=permutation_importance(model,Xte,yte,n_repeats=6,random_state=42,scoring="roc_auc")
        fi=pd.DataFrame({"feature":X.columns,"importance":pi.importances_mean,"sd":pi.importances_std}).sort_values("importance",ascending=False)
        st.plotly_chart(px.bar(fi.head(20).sort_values("importance"),x="importance",y="feature",orientation="h",error_x="sd"),use_container_width=True)
        st.dataframe(pd.DataFrame(classification_report(yte,pr,output_dict=True)).T,use_container_width=True)
        st.warning("Predictive importance is not causal evidence; lineage, linkage, sampling and phenotype quality can generate apparent predictors.")

with tabs[5]:
    st.subheader("Evidence Card")
    chosen=st.selectbox("Candidate",assoc.head(25).variant); r=assoc[assoc.variant==chosen].iloc[0]
    status="Stronger exploratory signal" if r.bonferroni<.05 else "Exploratory / requires validation"
    a,b=st.columns(2)
    with a:
        st.markdown(f"### {chosen}"); st.write(f"**p:** {r.p:.3g}"); st.write(f"**Effect:** {r.effect:.3f}"); st.write(f"**Adjusted p:** {r.bonferroni:.3g}")
        if "known_amr" in r: st.write(f"**Known AMR annotation:** {r.known_amr}")
    with b:
        st.write(f"**Evidence level:** {status}")
        st.write("**Supported:** prioritisation for follow-up in this dataset.")
        st.write("**Not supported:** causality, clinical utility, or a confirmed molecular mechanism.")
        st.write("**Next evidence:** independent replication, domain-appropriate association modelling, functional interpretation, and experimental validation.")
    card=f"Candidate: {chosen}\np: {r.p:.6g}\nEffect: {r.effect:.6g}\nAdjusted p: {r.bonferroni:.6g}\nEvidence: {status}\nNot a causal or clinical claim.\n"
    st.download_button("Download evidence card",card,f"{chosen.replace(':','_')}_evidence.txt")

with tabs[6]:
    st.markdown("""### Scope
**CohortOmics:** QC → PCA/population structure → exploratory association → polygenic-style score → evidence reporting.

**EvoResist-AI:** pathogen feature QC → lineage/population structure → association → held-out interpretable prediction → evidence cards.

### Limitations
- Demo datasets are simulated.
- v1 is not a replacement for PLINK2, SAIGE, REGENIE, GEMMA, pyseer or equivalent validated tools.
- Relatedness, LD, lineage, batch effects and phenotype quality need dataset-specific treatment.
- ML feature importance is not mechanistic or causal evidence.
- Experimental validation depends on the biological claim.
- No clinical or therapeutic decisions should be based on this application.
""")
