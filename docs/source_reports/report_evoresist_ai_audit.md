# EvoResist-AI Scientific Audit: Atria Alignment, Pathogen–Resistance System Selection, and Pilot Design

**Prepared for:** Yashwant Nama — VariantBridge EvoResist-AI branch (Atria University Life Sciences PhD application)
**Date:** 5 October 2026
**Scope:** Scientific methodology, literature and dataset audit. No production code. Every unverifiable claim is marked **UNVERIFIED**.
**Verification level:** Faculty claims from primary/authoritative web sources; dataset claims verified live via BV-BRC, ENA and NCBI APIs on 5 Oct 2026, plus a completed pilot download of 10 E. coli genomes with matched ciprofloxacin MIC phenotypes.

---

## 1. Atria University faculty / research-direction map

Legend: **[Led]** = personally led work (first/corresponding author or clear ownership); **[Co]** = co-authored; **[Inst]** = institutional theme or university description; **[Future]** = future-facing/aspirational description. Method relevance is stated conservatively: a past co-authored paper does not imply current use of that method.

### 1.1 Named faculty

**Dr Asha Velayudhan Nair** — Associate Professor, Life Sciences (at Atria since Aug 2024) [6]
- Training: PhD in Biochemistry, University of Bristol (Overseas Centenary Scholarship); postdocs at University of Cambridge and IIT Kharagpur [4, 12].
- **[Led/Co]** Molecular biology of multidrug transporters, all with H.W. van Veen and S. Murakami groups: NorM (V. cholerae MATE transporter) energetics, JBC 2014 [25]; LmrP (MFS) active-site plasticity, Sci Rep 2016 [25]; LmrA (ABC exporter) nucleotide energetics, Science Advances 2018 [25]; engineered MATE transporters, Communications Biology 2021 [25]; review chapter "Structural and functional landscape of MFS and MATE efflux pumps" [25]; drug-dependent inhibition of the S. pneumoniae ABC transporter PatAB [7].
- **[Inst]** Atria profile: antimicrobial resistance, drug delivery systems, biofilm-mediated resistance [1, 3]. Public talk "Decoding Resistance: Innovations Against AMR" (2026) covering biofilms, multi-omics and Big Data, Trojan-horse nanotechnology [12, 30].
- **[Led]** Actively recruiting a PhD student (advertised Apr 2025) on *mechanisms of multidrug transporters and their role in biofilm formation, and strategies to combat resistance*; desirable skills: bacterial culture, PCR, cloning [13]. Presented Atria work at NCCR-ABiS 2025 (DBT-sponsored national conference on combating AMR via biofilm inhibition) [29].
- **Method relevance (conservative):** wet-lab molecular microbiology and transporter biochemistry — efflux assays, biofilm assays, MIC testing, cloning. **No evidence she currently performs microbial GWAS or ML.** Her value to EvoResist-AI is biological interpretation of efflux/transporter hits and supervision of experimental validation.

**Dr Saurabh Mahajan** — Associate Professor and Head, Centre of Excellence in Life Sciences [2]
- Training: PhD in Biology, TIFR Mumbai [2].
- **[Led]** Bacterial genome evolution using large public genome sets and phylogenetic comparative methods: "Evolutionary jumps in bacterial GC content" (G3 2022, with D. Agashe, NCBS-TIFR) [32, 33]; codon-usage bias (GBE 2018) [31, 34]. **[Co]** E. coli/B. subtilis gene-regulatory-network modelling with metabolic feedback (FBA-based) [36]; SARS-CoV-2 variant impact on molecular diagnostics [31, 35].
- **[Inst/Future]** Heads the CoE; leads an open-source global biotechnology-tracking platform; industry project on systems-biology modelling of a bioleaching bacterium [2, 8].
- **Method relevance:** comparative bacterial genomics, population/phylogenetic thinking, bioinformatics and data science. The strongest computational-genomics alignment for EvoResist-AI. **No evidence he currently runs AMR GWAS specifically** — his genomics is evolutionary/comparative, which is adjacent, not identical.

**Dr Kristi Kabyashree** — Assistant Professor [1, 5]
- Training: PhD Tezpur University; postdoc NIPGR New Delhi [4].
- **[Led]** Virulence of the plant pathogen *Ralstonia solanacearum* and host interaction; interests: molecular biology, microbiology, bacterial genetics, plant–microbe interactions [4, 9].
- **Method relevance:** bacterial genetics and wet-lab microbiology. Relevant to experimental-validation culture, not to AMR genomics directly. Note: an Atria undergraduate intern annotated *R. solanacearum* genomes including "antimicrobial resistance traits" [14, 15] — this is **[Inst]** student project work, not evidence of faculty-led AMR genomics.

**Dr Abhishek Sinha** — Associate Professor [1, 16]
- Training: molecular parasitology (*Giardia*, proteasome biology); postdocs in Israel (VATAT) and the Netherlands (Oncode Institute, Leiden UMC) [16, 18].
- **[Inst]** Current stated interests: nanoscience, cancer biology, medicinal chemistry, computational biology [16]. **[Led/Inst]** Supervising an M.Sc. intern project "Assaying the Effectiveness of Small Molecules and Peptides Against Pathogenic Bacteria" (microbiological/biochemical antimicrobial assays, advertised Oct 2025) [44].
- **Method relevance:** antimicrobial compound screening and phenotypic assays — a plausible experimental-validation partner, but his publication record is not AMR.

**Dr Dhruti Patwardhan** — Assistant Professor, member of SynbioHub [11]
- **[Inst]** Leads development of a platform integrating curated AI models for synthetic biology; presented "AI for synthetic biology" at CDAC Accelerating Biology 2026 [11, 21]. Background: neuroscience (IISc), then product management/senior scientist at nference integrating genetics data into a clinical platform [19].
- **Method relevance:** AI/data-platform engineering and data-product discipline. Relevant to EvoResist-AI's tool engineering, not to AMR biology.

**Dr Sudha Warrier** — Professor and Program Director, Life Sciences [11, 17]
- Reproductive biology PhD (University of Madras); Duke postdoc; stem cells, cancer, regenerative medicine [11, 17]. **Institutional leadership role; not scientifically relevant to AMR genomics.**

### 1.2 Wider scan (separated from the named list)

- **Dr Soham Choudhuri** — Assistant Professor. PhD IIIT Hyderabad; postdoc TIFR Hyderabad. **[Led]** AI-driven drug discovery, protein/peptide generation models, single-cell multi-omics for target identification, AI virtual-cell models [39]. Method relevance: interpretable ML/AI methodology — relevant to the predictive-baseline module.
- **Dr Upasana Mohapatra** — Assistant Professor. Plant science/agriculture/biotechnology [40, 42]; supervising a bacterial-detection-kit internship project [44]. Marginal relevance.
- **[Inst] Institutional signals:** the CoE PhD programme states "AI x Biology is a persistent theme" across planetary/human/plant health, and welcomes self-defined topics [43]; CoE bioinformatics interns work on bacterial comparative genomics (DNA-repair enzymes, *Ralstonia* annotation) using HMMER/BLASTP, Python/R/Bash and GCP [14]. These are institutional themes and student projects, **not** faculty-led AMR-genomics programmes.

### 1.3 Implications for EvoResist-AI

1. The only faculty-led AMR research at Atria is **molecular efflux-transporter and biofilm biology (Dr Nair)** — a wet-lab, mechanism-first programme [13, 25, 29].
2. The computational strength is **bacterial comparative genomics and data science (Dr Mahajan)** plus **AI methods (Dr Choudhuri, Dr Patwardhan)** [8, 31–34, 39].
3. Therefore the defensible PhD narrative is: *computational population-genomic analysis of AMR that generates efflux-centred, experimentally testable hypotheses validated in a molecular microbiology lab* — exactly the handoff's "I do not want my research to end at a prediction" framing. A generic "AI predicts AMR" project would align with no one's actual research.

---

## 2. Comparison of candidate pathogen–resistance systems

All BV-BRC numbers below were queried live from the `genome_amr` API (laboratory-method records only, computational predictions excluded) on 5 Oct 2026. "Unique genomes" = distinct `genome_id`s.

| Criterion | **E. coli — ciprofloxacin** | **M. tuberculosis — rifampicin/isoniazid (CRyPTIC)** | **S. aureus — methicillin/fluoroquinolone** | **N. gonorrhoeae — azithromycin/cefixime** | **A. baumannii — colistin/carbapenem** |
|---|---|---|---|---|---|
| AMR mechanism | gyrA/parC QRDR mutations; PMQR (qnr, aac(6')-Ib-cr); AcrAB-TolC efflux upregulation via acrR/marR/soxR mutations [86, 87] | rpoB RRDR (codons 426–452), katG S315T, inhA promoter; exceptional mutation catalogue [50, 54] | mecA/SCCmec (near-deterministic for methicillin); grlA/gyrA + NorA efflux for fluoroquinolones | mtrCDE efflux + mtrR promoter/ coding mutations; penA mosaic; gyrA [UNVERIFIED — literature knowledge, not re-verified here] | AdeABC efflux; colistin heteroresistance (lipid A modification); OXA carbapenemases [UNVERIFIED here] |
| Atria faculty alignment | **Strong**: efflux biology (Nair), bacterial genomics (Mahajan), antimicrobial assays (Sinha), AI (Choudhuri/Patwardhan) | Weak: no TB faculty; BSL-3 blocks wet-lab work | Moderate: Gram-positive transporter link via Nair's PatAB work [7]; NorA efflux | Moderate (efflux-centric mtrCDE) but no culture capability | Moderate (AdeABC efflux) |
| Open isolate-level genomes | **Verified**: 10,904 unique genomes with lab AST in BV-BRC; assemblies + SRA reads open | **Verified**: 12,289 isolates, fully open (EBI FTP) [54] | Large public genome pools (e.g., ~99k NCBI assemblies [55]; Staphopia 43k [58]) but AST linkage is the bottleneck | **Verified**: 5,461 unique genomes with lab AST in BV-BRC (species taxid 485; an initial query under an incorrect taxid returned 0 and was re-checked); published collections also exist but **UNVERIFIED** here | BV-BRC: 384 colistin / 1,839 imipenem lab records (verified); modest |
| AST/MIC quality | 4,681 MIC + 6,538 disk-diffusion + 1,094 agar-dilution + 3,402 broth-dilution records; R/S calls present | **Best-in-class**: uniform MIC (UKMYC5/6 plates, 13 drugs), QC'd, binary calls [54] | BV-BRC: 0 MIC records for cipro/methicillin/cefoxitin; 254 oxacillin MIC | Deep but non-quantitative in BV-BRC: ceftriaxone 5,864; azithromycin 5,513; ciprofloxacin 2,844; cefixime 2,545 lab records — **0 MIC-method records** (categorical R/S or zone data) | BV-BRC: 0 MIC records for the drugs checked |
| Sample size (usable) | ~10,904 (AST-linked); ~4,700 MIC | 12,289 | Low hundreds to ~1,000 with AST (**UNVERIFIED** for specific studies, e.g., Earle 2016 collection) | ~5,461 AST-linked; 0 MIC | <400 (colistin) |
| Reference strain/build | E. coli K-12 MG1655 (ASM584v2 / GCF_000005845.2) or ST131 reference; standard, annotated | H37Rv NC_000962.3 (CRyPTIC uses it for catalogues) [54] | N315 / NCTC 8325 [59] | FA1090/NCCP11945 **UNVERIFIED** | ATCC 17978 **UNVERIFIED** |
| Genomic format | FASTA assemblies (draft, 74–162 contigs in pilot) + SRA reads; VCFs derivable | VCFs provided + ENA reads [54] | FASTA/reads | — | — |
| Phenotype format | MIC (mg/L) + CLSI R/S calls, per-isolate, with PMID provenance | MIC + binary R/S in CSV [54] | Mostly disk-diffusion zones (harder to model) | Categorical R/S (no MIC in BV-BRC) | — |
| Positive controls | gyrA S83L/D87N, parC S80I/E84V; acrR/marR truncations; qnr genes [86, 87] | rpoB 426–452, katG S315T, inhA promoter [50] | mecA presence; grlA S80F/Y, gyrA S84L | mtrR promoter deletions | adeRS mutations; lpx genes |
| Lineage/population structure | Diverse STs (pilot: ST10, ST131-like O25:H4, ST404, ST58...), strong structure — ideal for demonstrating structure control | Hyperclonal, low recombination — different (easier) structure problem | Clonal complexes (CC5/CC8/CC30...) | Strong structure + high recombination | IC1/IC2 clones |
| Microbial-GWAS feasibility | High (n large, MIC quantitative, structure demonstrable) | High (many published analyses) | Moderate | Moderate (large n, but categorical phenotypes only) | Low–moderate |
| Interpretable ML without leakage | Feasible with lineage-aware splits | Feasible | Feasible | Feasible (binary outcome) | Feasible but small n |
| Novelty potential | Efflux-regulator mutations (acrR/marR) are under-catalogued vs QRDR — realistic novelty space | Limited (catalogue near-saturated) | Limited for mecA; some for efflux | Some | Some |
| Testable follow-up hypotheses | Allele-swap, efflux assays, biofilm link — BSL-1/2 | BSL-3 — **not feasible at Atria** | BSL-2 feasible | Fastidious culture, specialized | BSL-2, heteroresistance assays tricky |
| Download/compute | ~5 MB/genome assembly; full 10.9k ≈ 55–60 GB; pilot subsets (500–2,000) run on 16 CPU/64 GB | VCFs + tables, tens of GB | — | — | — |
| Licensing/redistribution | BV-BRC public genomes (open); NCBI/ENA open data; no sequence redistribution restriction; cite source studies | Fully open ("fully open source" per consortium) [54] | Open | — | — |

### Verdict

- **E. coli–ciprofloxacin** is the only candidate that is simultaneously (a) verified at scale, (b) efflux-centric, (c) wet-lab-validatable at Atria, and (d) rich enough to demonstrate every methodology step honestly — in particular, it is the only large verified set with *quantitative* MIC data.
- **N. gonorrhoeae–azithromycin/cefixime** is a genuine runner-up: deep verified lab AST (5,461 genomes) and the most directly efflux-centric biology (mtrCDE/mtrR), but BV-BRC carries no MIC-method records for it (categorical phenotypes only), published-collection access is unverified, and experimental validation requires fastidious culture with no capability path at Atria.
- **M. tuberculosis–CRyPTIC** is the best pure dataset and an excellent *secondary benchmark* (recover known determinants in a hyperclonal organism), but fails the Atria experimental-alignment test (BSL-3).
- **S. aureus and A. baumannii** fail primarily on verified open AST/MIC depth (S. aureus: AST linkage bottleneck; A. baumannii: few hundred colistin records).

---

## 3. Recommended pilot system and justification

**Recommended: Escherichia coli — ciprofloxacin (fluoroquinolone) resistance.**

Justification:
1. **Verified data depth.** 10,904 unique genomes with laboratory AST for ciprofloxacin in BV-BRC (4,681 MIC records), traceable to major open surveillance studies: Lancet Microbe 2024 multi-country cohort (PMID 38219757; 4,286 records), Nat Commun 2023 national surveillance with standardized quantitative susceptibility testing (PMID 38052776; 2,429), Genome Res 2017 invasive E. coli England (PMID 28720578; 1,504), Lancet Microbe 2021 nosocomial transmission (PMID 34485958; 1,092). Pilot download of 10 isolates verified the genome↔phenotype join end-to-end (see §5).
2. **Efflux-centric answer key aligned with Atria.** Ciprofloxacin resistance is the textbook system where efflux-pump dysregulation (AcrAB-TolC via acrR/marR/soxR mutations) combines with target-site mutations (gyrA/parC QRDR) [86, 87]. Dr Nair's entire career is multidrug-transporter molecular biology [25]; the pilot's likely novel-candidate space (efflux-regulator variants, their epistasis with QRDR mutations) speaks directly to her PhD advertisement on multidrug transporters and biofilms [13].
3. **Experimental validation is genuinely feasible.** E. coli is BSL-1/2; MIC broth microdilution (CLSI), efflux assays (ethidium bromide accumulation with/without PAβN), allele replacement by recombineering or CRISPR, and biofilm assays (crystal violet) are standard molecular-microbiology experiments compatible with the Atria CoE lab as described in Dr Nair's recruitment advert (bacterial culture, PCR, cloning) [13].
4. **Methodological honesty is demonstrable.** The population is structured across many STs with strong lineage–phenotype confounding — exactly the condition under which naive GWAS and random train/test splits fail, so the pilot can *show* the failure and the correction rather than merely assert it.
5. **Compute is tractable.** ~5 MB per draft assembly; a 1,000-isolate pilot subset is ≈5 GB and runs pyseer comfortably on a 16-core/64 GB machine; the full 10.9k set (≈55–60 GB) is an HPC-scale but optional extension.

---

## 4. Exact open datasets and download/accession information

### 4.1 Primary: BV-BRC aggregated E. coli AST + genomes (verified 5 Oct 2026)

- **Source:** BV-BRC (PATRIC) AMR Phenotype data, https://www.bv-brc.org — API endpoint `https://www.bv-brc.org/api/genome_amr/` (RQL queries) and `https://www.bv-brc.org/api/genome_sequence/` (per-genome contig sequences as JSON; reconstruct FASTA from the `sequence` field).
- **Verified query:** `and(eq(taxon_id,562),eq(antibiotic,"ciprofloxacin"),eq(evidence,"Laboratory Method"))&limit(100000)` → 15,777 records / 10,904 unique genomes; 4,681 MIC records; phenotype calls: 2,066 Resistant / 6,559 Susceptible / 120 Intermediate / 7,032 uncalled.
- **Provenance (top PMIDs on records):** 38219757 (Lancet Microbe 2024), 38052776 (Nat Commun 2023), 28720578 (Genome Res 2017), 34485958 (Lancet Microbe 2021), 36165686 (AAC 2022).
- **Per-genome metadata available:** MLST (Achtman), cgMLST, CheckM completeness/contamination, serotype in genome_name (O:H), collection date, isolation country, SRA accession, genome status.
- **Pilot verification (completed):** 10 isolates (5 Resistant, MIC 16–256 mg/L; 5 Susceptible, ≤0.03 mg/L; all PMID 38219757, UK isolates with ERR SRA accessions) downloaded as FASTA (4.9–5.3 MB each, 74–162 contigs, 100% CheckM completeness, valid ACGTN content). Join key: BV-BRC `genome_id` (e.g., `562.22918`).
- **License:** BV-BRC serves public NCBI/ENA data; sequences are open. Redistribution of BV-BRC's *aggregated* database itself is not implied — cite BV-BRC and the source PMIDs; do not re-host the aggregate. (BV-BRC terms of use should be re-read at deployment; marked as a to-do, not a blocker.)
- **Caveat:** records mix studies, labs, methods (MIC vs disk diffusion) and years — this heterogeneity is a *feature for teaching batch effects* but must be handled (see §5, §8).

### 4.2 Secondary benchmark: CRyPTIC M. tuberculosis compendium (verified via publication)

- **Source:** CRyPTIC Consortium, PLOS Biology 2022, doi:10.1371/journal.pbio.3001721 [50, 54].
- **Data:** 12,289 isolates with matched WGS + MICs to 13 anti-tubercular drugs; 6,814 resistant to ≥1 drug; 2,129 RR/MDR/pre-XDR/XDR.
- **Access:** fully open at `ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/` — `CRyPTIC_reuse_table_20221019.csv` (MICs, binary R/S, phenotype QC metrics, ENA sample IDs, VCF paths), plus GENOTYPES/VARIANTS/MUTATIONS/SAMPLES/PHENOTYPES tables; code at github.com/kerrimalone/Brankin_Malone_2022 [54].
- **Use in EvoResist-AI:** known-answer benchmark in a hyperclonal organism (rpoB/katG recovery) — *Portfolio demonstration only* for the Atria narrative (no wet-lab follow-up possible).

### 4.3 Annotation databases (answer-key sources)

- **CARD (McMaster):** curated AMR models + ARO. **License nuance (verified):** data downloads are restricted to non-commercial academic use; commercial use requires a licence; the ARO ontology alone is CC-BY 4.0 [74, 75, 77]. Acceptable for an academic PhD pilot; do not bundle CARD data files in a redistributable app without checking terms.
- **NCBI AMRFinderPlus + RefGene database:** US-government public-domain software and data — the safer default for redistribution. **UNVERIFIED in this audit** (not queried live); verify current licence text at deployment.
- **ResFinder (CGE):** open-source; **UNVERIFIED** here.

---

## 5. Exact statistical/computational workflow

```
[1] Acquire: BV-BRC genome_amr (ciprofloxacin, laboratory evidence)
      + genome metadata (CheckM, MLST, country, date, SRA)
      + genome FASTA (assemblies)
[2] Phenotype harmonisation:
      - keep MIC records as primary (quantitative); convert to log2 MIC
      - R/S calls per CLSI breakpoints, recorded per record
      - study/batch covariate from PMID / source study
      - de-duplicate: one genome per isolate; flag mixed-methods records
[3] Genomic QC:
      - assembly QC: CheckM completeness ≥95%, contamination ≤2%,
        n50/contig-count sanity, genome length 4.5–5.8 Mb
      - remove lab strains/outliers; verify species identity (k-mer or ANI)
[4] Variation:
      - SNP/indel calls vs K-12 MG1655 (GCF_000005845.2) after read or
        assembly mapping; exclude recombinant regions (e.g., Gubbins) or
        use assembly-based variants (Ma et al. approach within pyseer)
      - gene presence/absence: panaroo or Roary pangenome
      - unitigs: unitig-caller → pyseer unitig mode
[5] Lineage/population structure:
      - mash/Jaccard distance tree + MLST/cgMLST lineage labels
      - PCA/MDS on gene presence-absence
      - lineage labels (ST or clonal-cluster) used as fixed effects;
        similarity matrix used for pyseer --lmm
[6] Known-determinant annotation (answer key, BEFORE association):
      - AMRFinderPlus/ResFinder/CARD-RGI screen for qnr, aac(6')-Ib-cr,
        efflux-pump alleles; explicit QRDR scan of gyrA (S83, D87) and
        parC (S80, E84)
[7] Association (pyseer; Lees et al. 2018 [71]):
      - unitig GWAS with --lmm (or lineage fixed effects + --pc-exclude)
      - gene PA GWAS (burden-style) with same structure control
      - report beta/SE/p per feature from the SAME model (fixing the
        v0.1 effect-vs-p mismatch flagged in the handoff)
[8] Multiple testing / structure control:
      - Bonferroni or BH-FDR over tested features (state which)
      - permutation-based significance where feasible
      - report genomic-inflation-style diagnostics adapted to mGWAS
[9] Interpretable predictive baseline:
      - elastic net on [gene PA + QRDR alleles + top unitigs], standardised
      - nested CV; feature selection inside folds
      - splits: (a) random split — run ONLY to demonstrate leakage;
        (b) lineage-blocked split (hold out entire STs) — the honest number
      - metrics: outer-fold AUC / balanced accuracy / MAE(log2 MIC)
[10] Comparison with known biology:
      - overlap of association hits with answer key (§6)
      - rank recovery of gyrA/parC among top hits
[11] Candidate prioritisation:
      - score = association evidence + structure-robustness + annotation
        novelty (not in CARD/AMRFinder) + efflux-pathway membership
[12] Evidence Card (schema: evidence_card_schema.json) →
[13] Experimentally testable hypothesis (template, §7)
```

**Software choice.** Use **pyseer** (unitig and gene-PA modes, LMM with population structure, lineage effects) [71] — the established microbial-GWAS tool — **not** human-GWAS packages: PLINK-style models assume unlinked loci and weak structure, which is exactly what bacteria violate. Complements: **Saber & Shapiro 2020** benchmark shows multi-locus elastic-net models recover causal variants better than single-locus tests in bacteria [66] — justifying the interpretable-ML baseline; **Earle et al. 2016** provides the lineage-effects framework [73]; **Collins & Didelot 2018** offers a phylogeny-aware alternative if recombination is severe [67]; **panfeed** can map k-mers to gene clusters for interpretability [72].

**Explicit treatment of the requested pitfalls:**

| Pitfall | Treatment |
|---|---|
| Bacterial clonality | Clonality makes features correlated genome-wide; use LMM/lineage effects, interpret hits as *lineage-associated* unless homoplasy is shown |
| Lineage confounding | Lineage fixed effects + LMM; report both unadjusted and adjusted statistics; flag hits that vanish after adjustment |
| Accessory vs SNP variation | Run both gene-PA and unitig/SNP GWAS separately; never merge them into one p-value space |
| Rare variants | Frequency filter (e.g., ≥1% or ≥5 isolates); report filtered-out count; rare-variant burden left as future work |
| Linkage | Do not claim independence of linked QRDR sites; cluster correlated features and report clusters |
| Phenotype imbalance | Report R/S/intermediate counts; use balanced-accuracy and class-weighted models; consider MIC as quantitative outcome (more power) |
| Train/test leakage | Random splitting in structured bacteria inflates performance because clones are shared across folds; the primary metric is lineage-blocked CV; random-split results shown only as the leakage demonstration |
| Multiple testing | Pre-registered correction (BH-FDR q<0.05 primary; Bonferroni sensitivity) |

**Classification of every proposed analysis:**

| Analysis | Classification |
|---|---|
| QC + structure description of the verified dataset | **Statistically defensible** |
| pyseer unitig/gene-PA GWAS with LMM + FDR, recovering gyrA/parC | **Statistically defensible** |
| Known-answer benchmark on CRyPTIC TB | **Statistically defensible** (as benchmark) |
| Elastic-net predictive baseline with lineage-blocked nested CV | **Statistically defensible** |
| Random-split performance reported as headline | **Not recommended** (leakage demonstration only) |
| "Novel resistance gene" claims from single-study associations | **Portfolio demonstration only** |
| Causal/mechanistic language for any association | **Not recommended** (violates Evidence Card rules) |
| PRS-style polygenic score for MIC | **Portfolio demonstration only** (no independent target cohort) |
| Any clinical claim (treatment guidance) | **Not recommended** |

---

## 6. Known-answer validation plan

**Answer key (must be recovered before any novel candidate is surfaced):**
1. **gyrA QRDR** mutations (S83L, S83A, D87N/Y/H) and **parC** (S80I, E84V/K) — the dominant ciprofloxacin determinants in E. coli; a systematic review quantifies gyrA/parC as producing the largest MIC increases, with efflux contributing smaller gains and combined mechanisms reaching ~250–4,000-fold [86, 87].
2. **PMQR genes** where present: qnrA/B/S, aac(6')-Ib-cr — lower-level resistance [87].
3. **Efflux dysregulation:** non-synonymous/truncating mutations in **acrR** and **marR** (local repressors of acrAB and the mar operon) — expected as secondary hits with smaller effect sizes [86, 87].
4. **Cross-check:** CARD/AMRFinderPlus annotation of the same genomes must reproduce the known-gene list; disagreement between GWAS hits and annotation is itself a reportable QC result.

**Pass criteria (pre-registered):**
- gyrA/parC QRDR alleles appear among top FDR-significant unitig/SNP hits in the structure-adjusted model.
- Their effect direction matches literature (QRDR alleles → higher MIC).
- In the CRyPTIC benchmark, rpoB RRDR and katG S315T are recovered as top hits.
- Predictive baseline achieves outer-fold AUC materially above chance **under lineage-blocked CV**; the gap between random-split and lineage-blocked AUC is reported as the leakage lesson.

**Failure handling:** if the answer key is not recovered, the pipeline is not ready to interpret anything novel — diagnose (phenotype noise, structure control, calling artefacts) before proceeding. This is the handoff's "recover established determinants first" rule made operational.

---

## 7. Experimental-validation roadmap appropriate to an Atria PhD

Framed to Dr Nair's actual capabilities (bacterial culture, PCR, cloning, transporter assays, biofilm assays) [13, 25, 29]:

**Phase 0 (in silico, months 0–6).** This pilot: verified data, GWAS, answer-key recovery, Evidence Cards. Deliverable: methods + honest negative results.

**Phase 1 (in vitro confirmation, months 6–18).**
- Assemble a small local panel: 30–50 E. coli isolates (clinical/environmental collaborators or culture-collection strains) — *requires a source agreement; do not claim feasibility without it*.
- Broth-microdilution MIC per CLSI (ciprofloxacin); WGS (outsourced); QRDR Sanger sequencing.
- Test: do local isolates reproduce the genotype–phenotype map learned from public data? (External validation of the *predictive* claim, not the causal one.)

**Phase 2 (mechanistic, months 18–36).**
- **Allele replacement** (recombineering or CRISPR) of candidate efflux-regulator alleles (e.g., acrR truncation) into a clean E. coli background; MIC shift measured. This is the only step that supports causal language, and only for the mutated allele in that background.
- **Efflux phenotyping:** ethidium-bromide accumulation/efflux ± PAβN (EPI) to show the phenotype is efflux-mediated.
- **Biofilm link (Nair alignment):** crystal-violet biofilm assays on isogenic pairs differing at candidate alleles — connects AMR genotype to the biofilm-mediated resistance theme of the host lab [13, 29].

**Phase 3 (generalisation, optional).** Repeat for a second antibiotic (e.g., chloramphenicol or tetracycline — classic AcrAB-TolC substrates) to test substrate-class hypotheses.

**Honesty constraints:** Phase 1 requires an isolate source and sequencing budget (UNVERIFIED — must be confirmed with the supervisor). Phase 2 requires molecular-microbiology bench capability that matches the skills listed in Dr Nair's PhD advert [13]. Nothing in the computational pilot should be described as "validated experimentally" until Phase 2 data exist.

---

## 8. Risks and limitations

1. **Phenotype heterogeneity.** BV-BRC aggregates studies with different labs, protocols, breakpoints and years. Mixing MIC and disk-diffusion is statistically hazardous; the pilot must either restrict to MIC records or model study as a batch covariate. Residual batch effects are the top threat to validity.
2. **Aggregation provenance.** BV-BRC records carry PMIDs but the underlying per-study QC is not re-verifiable from the aggregate alone; for publication-grade claims, re-derive phenotypes from the primary studies (e.g., the Nat Commun 2023 standardized-susceptibility subset).
3. **Structure confounding is severe.** In E. coli, STs differ in both resistance prevalence and genome content; unadjusted GWAS will return lineage markers. The design treats this as a teaching outcome, but it caps mechanistic interpretability of accessory-genome hits.
4. **Efflux-regulator effects are small and polygenic.** acrR/marR effects are modest vs QRDR [86]; the pilot may honestly fail to rank them significantly — this must be reported as a negative result, not hidden.
5. **Novelty risk.** With ~10k genomes already deeply mined, "novel" candidates are likely rare-effect or annotation-artefact findings; the Evidence Card's "claim NOT justified" field exists for this.
6. **Compute.** Full 10.9k-genome unitig GWAS needs >64 GB RAM for LMM kinship (pyseer --lmm memory scales with n²); plan HPC or a 1–2k subset for the sandbox pilot.
7. **Licensing.** CARD data files are non-commercial-only [74, 77]; use AMRFinderPlus (public domain, UNVERIFIED) for anything redistributable.
8. **NCBI Pathogen Detection API** was unreachable from this sandbox (DNS blocked) — its AST tables could not be verified here (**UNVERIFIED**); BV-BRC was used instead.
9. **No clinical claims.** Breakpoint calls in the data are surveillance labels, not treatment guidance.

---

## 9. EvoResist-AI input/output schemas

### Inputs

```yaml
# config.yaml
dataset:
  name: bvbrc_ecoli_cipro
  source_api: https://www.bv-brc.org/api/genome_amr/
  query: 'and(eq(taxon_id,562),eq(antibiotic,"ciprofloxacin"),eq(evidence,"Laboratory Method"))'
  citation: "BV-BRC (PATRIC/BV-BRC); source PMIDs: 38219757, 38052776, 28720578, 34485958"
  license: "Open public genomes; cite BV-BRC + source studies; no aggregate re-hosting"
  download_date: null          # filled by pipeline
reference:
  name: Escherichia coli K-12 MG1655
  accession: GCF_000005845.2   # ASM584v2
phenotype_table:               # one row per isolate
  columns:
    genome_id: {type: str, key: true}        # BV-BRC genome_id, e.g. 562.22918
    mic: {type: float|str, nullable: true}   # mg/L; "<=0.03" style strings parsed
    resistant_phenotype: {type: enum[Resistant, Susceptible, Intermediate, null]}
    laboratory_typing_method: {type: enum[MIC, Disk diffusion, Agar dilution, Broth dilution]}
    study_pmid: {type: str}
    batch: {type: str}                       # derived from PMID/study
manifest:
  qc:
    checkm_completeness_min: 95
    checkm_contamination_max: 2
    genome_length_range_mb: [4.5, 5.8]
  association:
    tool: pyseer
    mode: [unitig, gene_pa]
    structure: lmm            # or lineage_fixed_effects
    maf_min: 0.01
    correction: bh_fdr
    alpha: 0.05
  ml:
    model: elastic_net
    cv: nested_lineage_blocked
    leakage_demo_random_split: true   # reported separately, never as headline
```

### Outputs

```yaml
association_results:           # one row per feature
  columns:
    feature_id: str            # unitig / gene / SNP
    feature_type: enum[unitig, gene_pa, snp]
    beta: float
    standard_error: float
    ci_low: float
    ci_high: float
    p_value: float
    q_value: float             # BH-FDR
    lineage_adjusted: bool
    unadjusted_p: float        # for transparency
    annotation_source: str     # AMRFinderPlus/CARD/none
    in_answer_key: bool
predictive_baseline:
  outer_fold_metrics: list     # per-fold AUC, balanced accuracy
  split_strategy: enum[lineage_blocked, random_demo]
  top_features: list           # with coefficients, selection inside folds
evidence_cards:                # see evidence_card_schema.json
provenance_manifest:           # per handoff §12: dataset, URL, citation, license,
                               # download date, QC thresholds, covariates, model,
                               # seed, versions, analysis date, limitations
```

---

## 10. Implementation Specification for Claude

Claude implements engineering only; the scientific methodology above is fixed and must not be silently redesigned.

**Module layout** (extends the existing `src/variantbridge/` package; no changes to CohortOmics modules):

```
src/variantbridge/amr/
├── data_io.py        # BV-BRC fetchers: genome_amr query, genome metadata,
│                     # genome_sequence JSON→FASTA reconstruction; caching to data/
├── phenotypes.py     # MIC string parsing ("<=", ">="), log2 transform, R/S
│                     # harmonisation, batch covariate from PMID, de-duplication
├── qc.py             # CheckM/length/contig filters; species sanity check
├── structure.py      # mash/sketch distances, tree, PCA on gene PA, lineage labels
├── annotation.py     # AMRFinderPlus (default) / CARD-RGI (opt-in) wrappers;
│                     # QRDR scan of gyrA/parC against reference coordinates
├── association.py    # pyseer subprocess wrapper (unitig + gene_pa + --lmm);
│                     # parses pyseer TSV into association_results schema;
│                     # effect/SE/CI/p from the SAME model (fixes v0.1 mismatch)
├── prediction.py     # elastic-net baseline; nested lineage-blocked CV;
│                     # random-split demo flagged leakage_demo, reported separately
├── evidence.py       # builds Evidence Cards per evidence_card_schema.json;
│                     # enforces claim-justified / claim-NOT-justified logic
└── reporting.py      # provenance manifest (handoff §12 fields), methods report
```

**Interfaces:** `data_io.fetch_phenotypes(config) -> DataFrame`; `association.run_pyseer(features, phenotypes, structure, config) -> association_results`; `evidence.build_card(hit, annotation, prediction, literature) -> dict` validating against the JSON schema. All long-running steps write checkpoints to disk and are resumable.

**Tests (known-answer, small, deterministic):**
1. `test_phenotypes.py`: MIC string parsing (`"<=0.03"`, `"256"`, `">32"`), log2 transform, imbalance report.
2. `test_qc.py`: synthetic genomes failing each QC gate are excluded.
3. `test_association.py`: on a synthetic panmatrix with a planted determinant, pyseer wrapper + FDR recovers the planted feature as top hit; effect/SE/CI/p come from one model.
4. `test_prediction.py`: lineage-blocked CV on synthetic structured data yields lower AUC than random split — asserting the leakage demonstration works.
5. `test_evidence.py`: a card with an unadjusted-only association cannot claim "claim_justified: causal/mechanistic"; schema validation rejects cards missing uncertainty fields.
6. Pilot smoke test: the 10 verified pilot genomes + `pilot_phenotypes.csv` run through QC→structure→(mock) association end-to-end.

**Non-negotiables carried over from the handoff:** no fabricated results; no causal language; no random-split headline metrics; provenance manifest on every run; negative results reported; no production code before this audit's approval — this document is that approval basis.

---

### Summary of verified vs unverified items

| Item | Status |
|---|---|
| BV-BRC E. coli ciprofloxacin AST counts, unique genomes, provenance PMIDs | **Verified live (5 Oct 2026)** |
| Pilot 10-genome download + phenotype join + FASTA integrity + metadata | **Verified live** |
| CRyPTIC compendium contents and open FTP access | **Verified via publication** [54] |
| CARD licence split (data non-commercial; ARO CC-BY 4.0) | **Verified** [74–77] |
| Faculty roles, training, led vs co-authored work | **Verified from primary/authoritative pages** |
| NCBI Pathogen Detection AST availability | **UNVERIFIED** (API unreachable from sandbox) |
| AMRFinderPlus / ResFinder licence texts | **UNVERIFIED** |
| Earle 2016 S. aureus dataset access conditions | **UNVERIFIED** |
| Atria wet-lab facility scope and isolate-source agreements | **UNVERIFIED** (must be confirmed with supervisor) |
