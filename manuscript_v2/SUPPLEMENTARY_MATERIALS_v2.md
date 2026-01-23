# Supplementary Materials v2.0

**Title**: Ethylene-Induced Isoflavonoid Biosynthesis in Soybean: A Two-Track Multi-Omics Analysis  
**Version**: 2.0 (Two-Track Restructured)  
**Date**: 2026-01-21

---

## Supplementary Figures

### Figure S1. Database Coverage and Metabolite Categorization

**A.** KEGG database mapping coverage showing 36.7% of metabolites (29/79) successfully mapped.

**B.** Distribution of significant vs. non-significant metabolites stratified by KEGG mapping status.

**C.** Top 10 unmapped significant metabolites ranked by P-value (6''-O-acetyldaidzin, 6''-malonylgenistin, etc.).

**D.** Distribution of log₂ fold changes for KEGG-mapped vs. unmapped metabolites.

---

### Figure S2. Two-Track Analysis Framework Overview (NEW)

Schematic diagram illustrating the conceptual separation between Track A and Track B analyses:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         TWO-TRACK ANALYSIS FRAMEWORK                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  TRACK A: Biosynthetic Pathway Analysis                              │   │
│  │                                                                      │   │
│  │  Question: "Which enzymes synthesize these metabolites?"             │   │
│  │                                                                      │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────────┐  │   │
│  │  │ KEGG Reactions   │───→│      GNN         │───→│ Prioritized   │  │   │
│  │  │ STRING PPI       │    │     (HGT)        │    │   Enzymes     │  │   │
│  │  └──────────────────┘    └──────────────────┘    └───────┬───────┘  │   │
│  │                                                          │          │   │
│  │                                                          ▼          │   │
│  │                                               ┌───────────────────┐ │   │
│  │                                               │   Proteomics      │ │   │
│  │                                               │   Validation      │ │   │
│  │                                               │   (Expression)    │ │   │
│  │                                               └───────────────────┘ │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  TRACK B: Metabolite-Protein Interaction Screening                   │   │
│  │                                                                      │   │
│  │  Question: "Do metabolites bind to and regulate other proteins?"     │   │
│  │                                                                      │   │
│  │  ┌──────────────────┐    ┌──────────────────┐    ┌───────────────┐  │   │
│  │  │ Metabolite 3D    │───→│   Molecular      │───→│ Predicted     │  │   │
│  │  │ Protein Struct.  │    │   Docking        │    │   Binding     │  │   │
│  │  └──────────────────┘    └──────────────────┘    └───────┬───────┘  │   │
│  │                                                          │          │   │
│  │                                                          ▼          │   │
│  │                                               ┌───────────────────┐ │   │
│  │                                               │   Experimental    │ │   │
│  │                                               │   Validation      │ │   │
│  │                                               │   (SPR/ITC/MST)   │ │   │
│  │                                               │   ★ REQUIRED ★    │ │   │
│  │                                               └───────────────────┘ │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ⚠️ These tracks are INDEPENDENT and answer DIFFERENT questions            │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Figure S3. Track A: GNN Performance and Validation

**A.** Hits@20 comparison across methods:
- Random baseline: 5.8%
- Adamic-Adar heuristic: 14.9%
- HGT (Enhanced): 77.6%

**B.** Ablation study results showing contribution of different graph components.

**C.** Proteomics validation: Fold change of GNN-prioritized enzymes in ethylene vs. control samples.

**D.** Correlation between GNN prediction rank and proteomics fold change (r=0.78, P<0.01).

---

### Figure S4. Track A: Biosynthetic Pathway Activation

**A.** Complete isoflavonoid biosynthesis pathway with enzyme fold changes overlaid:
```
Phenylalanine
      │
      ▼ PAL (13×↑)
Cinnamic acid
      │
      ▼ C4H
p-Coumaric acid
      │
      ▼ 4CL (15×↑)
p-Coumaroyl-CoA
      │
      ▼ CHS (7×↑)
Naringenin chalcone
      │
      ▼ CHI (34×↑)
Naringenin
      │
      ▼ IFS (9×↑)
Genistein/Daidzein
      │
      ▼ IFR (84×↑)
Isoflavonoid products
```

**B.** Enzyme-metabolite coherence analysis (all 6 pairs show concordant upregulation).

**C.** Fisher combined P-value calculation demonstrating multi-omics convergence.

---

### Figure S5. Track B: Docking Predictions (Exploratory)

> [!CAUTION]
> Track B predictions are hypothetical and require experimental validation.
> These results are INDEPENDENT from Track A and do NOT validate GNN predictions.

**A.** Binding poses for Daidzein-FNR interaction (-7.80 kcal/mol).

**B.** Binding poses for Formononetin-Kinase interaction (-7.60 kcal/mol).

**C.** Distribution of docking scores across all tested metabolite-protein pairs.

**D.** Comparison with known ligand binding sites (when available).

---

### Figure S6. Methodological Clarification: Why Docking ≠ GNN Validation (NEW)

Conceptual diagram explaining why molecular docking cannot validate GNN predictions:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│        WHY DOCKING CANNOT VALIDATE GNN-PREDICTED ENZYME RELATIONSHIPS       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  GNN Prediction:                                                            │
│  "IFS is functionally related to Daidzein"                                  │
│       ↓                                                                     │
│  Correct interpretation: IFS catalyzes Daidzein synthesis                   │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                     ENZYMATIC CATALYSIS                              │   │
│  │                                                                      │   │
│  │     Liquiritigenin (기질)  ───IFS───→  Daidzein (생성물)              │   │
│  │           │                               │                          │   │
│  │           ▼                               ▼                          │   │
│  │      IFS가 결합함              IFS에서 방출됨 (결합 안 함)            │   │
│  │                                                                      │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│  ❌ Wrong validation approach:                                              │
│     Docking IFS with Daidzein                                               │
│     → Tests: "Can the PRODUCT bind to the ENZYME?"                          │
│     → This tests product inhibition, NOT catalytic activity!                │
│                                                                             │
│  ✓ Correct validation approaches:                                           │
│     1. Proteomics: Is IFS upregulated when Daidzein increases?              │
│     2. Enzyme assay: Does purified IFS convert substrate → Daidzein?        │
│     3. Genetics: Does IFS knockout eliminate Daidzein?                      │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Figure S7. Quality Control Summary

**A.** P-value distribution histogram.

**B.** Log₂ fold change distribution.

**C.** MA plot colored by significance.

**D.** Effect size distribution (Cohen's d).

---

## Supplementary Tables

### Table S1. Complete KEGG Pathway Enrichment Results

All 44 KEGG pathways tested with P-values, FDR, Bonferroni corrections.

**File**: `results/kegg_pathway_statistical_enhanced.csv`

---

### Table S2. Complete PlantCyc Pathway Enrichment Results

All 268 PlantCyc pathways tested.

**File**: `results/plantcyc_pathway_statistical_enhanced.csv`

---

### Table S3. Complete Differential Metabolite Analysis

All 79 metabolites with fold changes and P-values.

**File**: `results/Supplementary_Table_S3_All_Metabolites.csv`

---

### Table S4. Track A: GNN Model Details (NEW)

| Parameter | Value |
|-----------|-------|
| **Graph Structure** | |
| Enzyme nodes | 3,425 |
| Metabolite nodes | 200 |
| Tier-R edges (reaction-grounded) | 372 |
| Tier-P edges (pathway-supported) | 5,528 |
| PPI edges | ~50,000 |
| **Model Architecture** | |
| Embedding dimension | 64 |
| Number of layers | 2 |
| Attention heads | 4 |
| **Training** | |
| Split strategy | Node-disjoint (10% test enzymes) |
| Optimizer | Adam (lr=0.005) |
| Epochs | 100 |
| **Performance** | |
| Hits@20 | 77.6% |
| Improvement over random | 13.4× |

---

### Table S5. Track A: GNN-Proteomics Validation (NEW)

| GNN-Predicted Enzyme | Prediction Rank | Proteomics FC | P-value | Validated |
|---------------------|-----------------|---------------|---------|-----------|
| IFS (Isoflavone synthase) | Top-5 | 9× | 0.006 | ✓ |
| IFR (Isoflavone reductase) | Top-10 | 84× | 0.037 | ✓ |
| CHI (Chalcone isomerase) | Top-10 | 34× | 0.047 | ✓ |
| CHS (Chalcone synthase) | Top-15 | 7× | <0.05 | ✓ |
| 4CL (4-Coumarate:CoA ligase) | Top-15 | 15× | <0.05 | ✓ |
| PAL (Phenylalanine ammonia-lyase) | Top-20 | 13× | <0.05 | ✓ |

**Validation rate**: 6/6 (100%) of GNN-prioritized enzymes show significant upregulation

---

### Table S6. Track B: Docking Results (Exploratory) (NEW)

> [!WARNING]
> These predictions require experimental validation. Docking does not constitute proof of biological interaction.

| Metabolite | Protein Target | UniProt ID | Binding Energy (kcal/mol) | Proposed Mechanism |
|------------|----------------|------------|--------------------------|-------------------|
| Daidzein | FNR | A0A0R4J4B4 | -7.80 | Allosteric modulation of redox sensing |
| Formononetin | Ser/Thr Kinase | I1JK97 | -7.60 | Feedback regulation of signaling |
| Genistein | --- | --- | --- | Not screened |

**Required validation experiments**:
- Surface Plasmon Resonance (SPR)
- Isothermal Titration Calorimetry (ITC)
- Microscale Thermophoresis (MST)
- CRISPR knockout studies

---

### Table S7. Enzyme Abundance Changes in Isoflavonoid Biosynthesis

| Enzyme | Gene ID | Log₂ FC | Linear FC | P-value | Pathway Role |
|--------|---------|---------|-----------|---------|--------------|
| PAL | Glyma.02G042500 | 3.72 | 13× | <0.05 | Entry point |
| 4CL | Glyma.11G070500 | 3.89 | 15× | <0.05 | CoA activation |
| CHS | Glyma.01G228700 | 2.89 | 7× | <0.05 | Chalcone synthesis |
| CHI | Glyma.10G292200 | 5.08 | 34× | <0.05 | Ring closure |
| IFS | Glyma.13G173500 | 3.22 | 9× | 0.006 | Isoflavone core |
| IFR | Glyma.01G211800 | 6.39 | 84× | 0.037 | Reduction |

---

## Supplementary Methods

### Track A: GNN Analysis Details

**Graph Construction**:
1. STRING PPI edges filtered at confidence ≥700, textmining excluded
2. KEGG gene-EC mapping via REST API (`/link/ec/gmax`)
3. EC-reaction-metabolite mapping with currency metabolite filtering
4. Tier-R edges: Direct reaction-grounded (weight=1.0)
5. Tier-P edges: Same EC class first 3 digits (weight=0.5)

**Node Features**:
- All node types initialized with 64-dim random normal vectors
- Learnable embeddings updated during training

**Training Protocol**:
- Link prediction on Enzyme-Metabolite edges
- 10% validation, 10% test (node-disjoint)
- Negative sampling ratio: 1:1
- Binary cross-entropy loss

**Validation**:
- Computational: Hits@K on held-out edges
- Biological: Proteomics expression correlation

---

### Track B: Docking Analysis Details

**Protein Preparation**:
- AlphaFold2 structures (v4/v6) downloaded
- Hydrogen addition and charge assignment (OpenBabel)
- Active site prediction (when known) or blind docking

**Ligand Preparation**:
- 3D conformers from KEGG MOL files
- Gasteiger partial charges
- Energy minimization (MMFF94)

**Docking Protocol**:
- AutoDock Vina v1.2
- Exhaustiveness: 8
- Grid: Entire protein surface (blind docking)
- Output: Top 10 poses per ligand

**Interpretation Guidelines**:
- <-7.0 kcal/mol: Strong predicted binding
- -7.0 to -5.0 kcal/mol: Moderate binding
- >-5.0 kcal/mol: Weak/no binding

> [!NOTE]
> Docking scores correlate poorly with experimental Kd values.
> These predictions are hypothesis-generating only.

---

## Supplementary Note: Relationship Between Track A and Track B

### Why These Tracks Are Separate

Track A and Track B address **fundamentally different biological questions**:

| Aspect | Track A | Track B |
|--------|---------|---------|
| **Question** | Which enzymes synthesize metabolites? | Do metabolites regulate proteins? |
| **Input** | Metabolite of interest | Metabolite structure |
| **Method** | GNN on pathway network | Molecular docking |
| **Output** | Prioritized enzyme list | Binding predictions |
| **Validation** | Proteomics expression | Binding assays (SPR/ITC) |
| **Interpretation** | Biosynthetic pathway | Allosteric/signaling roles |

### Common Misconception

A frequent error in computational biology is attempting to "validate" GNN-predicted enzyme-metabolite relationships using docking:

```
❌ INCORRECT:
GNN predicts "IFS related to Daidzein"
→ Dock Daidzein into IFS active site
→ Good docking score = "validation"

This is wrong because:
- IFS binds SUBSTRATE (Liquiritigenin), not PRODUCT (Daidzein)
- Good docking to product suggests product INHIBITION, not synthesis
- The GNN prediction is about pathway membership, not physical binding
```

### Correct Validation Strategies

| Track | Correct Validation | Why |
|-------|-------------------|-----|
| Track A | Proteomics | Does predicted enzyme show expression change? |
| Track A | Enzyme assay | Does purified enzyme catalyze the reaction? |
| Track A | Genetics | Does knockout eliminate the metabolite? |
| Track B | SPR/ITC | Does metabolite physically bind protein? |
| Track B | Mutagenesis | Does binding site mutation abolish effect? |
| Track B | Cell assay | Does metabolite treatment affect protein function? |

---

## Data Availability

- **Metabolomics**: MetaboLights MTBLS531
- **Proteomics**: PRIDE PXD006989
- **STRING PPI**: STRING v12.0 (taxon 3847)
- **KEGG Pathways**: KEGG API (accessed 2026-01)
- **Analysis Code**: [Repository URL]
- **GNN Model**: Available upon request

---

*Supplementary materials v2.0 complete: 2026-01-21*
