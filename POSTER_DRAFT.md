# Conference Poster Draft
**Title**: Ethylene Primes Isoflavonoid Biosynthesis in Soybean: Integrating Multi-Omics with Graph Neural Networks to Uncover Metabolite-Mediated Signaling

---

## 📌 1. Background & Hypothesis
- **Ethylene's Role**: Key stress hormone, but its system-wide impact on secondary metabolism is complex.
- **Problem**: Traditional omics yield lists of genes/metabolites but miss *regulatory mechanisms*.
- **Hypothesis**: Ethylene signaling orchestrates a **coordinated system-wide reprogramming** of secondary metabolism, specifically priming the isoflavonoid pathway for stress defense via novel regulator interactions.

## 🛠️ 2. Workflow (Methodology)
*(Visual Suggestion: Flowchart connecting the three blocks)*
1.  **Multi-Omics**:
    - Metabolomics (LC-MS/MS, n=79)
    - Proteomics (Shotgun, n=6000+)
2.  **Graph Learning (GNN)**:
    - Heterogeneous Graph Transformer (HGT)
    - **Enhanced Graph**: STRING PPI + KEGG + Tier-P + Exp. Mets
3.  **Validation**:
    - **Docking**: AutoDock Vina + AlphaFold2 (v6)

## 📊 3. Key Results

### A. Massive & Specific Metabolic Priming
- **Isoflavonoid Conjugates Skyrocket**:
    - 6''-O-acetyldaidzin: **~5,000-fold** (P < 10⁻⁸)
    - 6''-malonylgenistin: **~4,300-fold**
    - Aglycones (active forms): Unchanged (~1.1-fold)
- **Conclusion**: Plants stockpile "ready-to-use" defense compounds (Priming).

### B. Coordinated Multi-Omics
- **Enzymes match Metabolites**: All pathway enzymes (PAL, CHS, IFS, IFR) significantly upregulated (7× ~ 84×).
- **Fisher's Combined P**: **1 × 10⁻¹²** (Extremely non-random activation).

### C. GNN & Novel Signaling Hypotheses
*(Visual: Network graph showing Daidzein connecting to FNR)*
- **Prioritized Screening**: Filtered candidate pool by **13-fold** (vs Random), prioritizing biologically relevant enzymes (Hits@20 = 77.6% vs 5.8%).
- **Novel Predictions (Docking Validated)**:
    1.  **Daidzein ↔ FNR** (-7.8 kcal/mol): Accumulating Daidzein may modulate chloroplast redox sensing (**Retrograde Signaling**).
    2.  **Formononetin ↔ Kinase** (-7.6 kcal/mol): Direct **Allosteric Feedback** to signaling pathways.

## 📝 4. Discussion & Impact
- **System-Level Insight**: Ethylene doesn't just "turn on" genes; it reconfigures the cellular signaling landscape via metabolite effectors.
- **Methodological Advance**: First successful application of **GNN + Docking** to discover non-genomic regulation in plant metabolomics.
- **Limitations**: In vitro validation (e.g., binding assays) is required to confirm predicted interactions. Phenotype data (mutants) is the next critical step.
- **Application**: Engineering stress-resilient crops with enhanced "priming" capabilities.

---
**Contact**: [User Name] | [Lab/Affiliation]
**Code**: https://github.com/[Lab]/ethylene-gnn
