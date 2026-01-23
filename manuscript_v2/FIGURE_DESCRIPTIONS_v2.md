# Figure Descriptions for Manuscript v2.0

**Version**: 2.0 (Two-Track Restructured)  
**Date**: 2026-01-21

---

## Main Text Figures

### Figure 1. Study Overview and Two-Track Analysis Framework

**Panel A**: Experimental design showing ethylene treatment and sample collection timeline.

**Panel B**: Two-track analysis framework:
- **Track A**: GNN-based biosynthetic enzyme prioritization → Proteomics validation
- **Track B**: Molecular docking → Experimental validation (required)

**Panel C**: Data integration overview showing metabolomics (MTBLS531), proteomics (PXD006989), STRING PPI, and KEGG pathways.

**Key Message**: This study employs two complementary but independent analytical approaches.

---

### Figure 2. Metabolomics Results: Isoflavonoid Accumulation

**Panel A**: Volcano plot showing differential metabolite abundance (x-axis: Log2FC, y-axis: -log10(P-value)).
- Highlight: Conjugated isoflavonoids (Log2FC ~12)
- Highlight: Basal aglycones (Log2FC ~0.14, highly significant)

**Panel B**: Bar chart comparing fold changes of conjugated vs. basal isoflavonoids.

**Panel C**: Heatmap of top 20 significantly changed metabolites.

---

### Figure 3. Pathway Enrichment Analysis

**Panel A**: KEGG pathway enrichment results (significant: map01110, P=0.030).

**Panel B**: PlantCyc pathway concordance showing biological consistency despite non-significance.

**Panel C**: Venn diagram of KEGG-mapped vs. unmapped metabolites.

---

### Figure 4. Track A: GNN Model Architecture and Performance

**Panel A**: Heterogeneous graph structure:
- Node types: Enzyme, Metabolite
- Edge types: PPI, Enzyme-Metabolite (Tier-R, Tier-P)

**Panel B**: HGT model architecture schematic.

**Panel C**: Performance comparison (Hits@20):
- Random: 5.8%
- Adamic-Adar: 14.9%
- HGT: 77.6%

**Panel D**: Ablation study results.

---

### Figure 5. Track A: Proteomics Validation of GNN Predictions

**Panel A**: Enzyme fold changes in isoflavonoid biosynthesis pathway:
```
PAL (13×) → 4CL (15×) → CHS (7×) → CHI (34×) → IFS (9×) → IFR (84×)
```

**Panel B**: Correlation between GNN prediction rank and proteomics fold change.

**Panel C**: Enzyme-metabolite coherence analysis (all pairs concordant).

**Panel D**: Fisher combined P-value calculation (P = 1×10⁻¹²).

**Key Message**: GNN predictions are validated by proteomics expression data, not docking.

---

### Figure 6. Track B: Metabolite-Protein Interaction Predictions (Exploratory)

> [!WARNING]
> Track B results are independent from Track A and require experimental validation.

**Panel A**: Daidzein-FNR docking pose with binding energy (-7.80 kcal/mol).

**Panel B**: Formononetin-Kinase docking pose (-7.60 kcal/mol).

**Panel C**: Proposed biological mechanisms (hypothetical):
- Daidzein → FNR → Chloroplast redox sensing
- Formononetin → Kinase → Signaling feedback

**Panel D**: Required experimental validation:
- SPR, ITC, MST for binding confirmation
- CRISPR knockouts for functional validation

---

### Figure 7. Summary: Relationship Between Tracks A and B

**Panel A**: Conceptual diagram showing the different questions addressed by each track:
- Track A: "Which enzymes synthesize metabolites?" ← Pathway analysis
- Track B: "Do metabolites regulate other proteins?" ← Signaling hypothesis

**Panel B**: Warning diagram explaining why docking cannot validate GNN:
```
IFS + Liquiritigenin (substrate) → Daidzein (product)
      ↑ binds here                    ↑ released here

Docking IFS + Daidzein tests "product inhibition"
NOT "biosynthetic capability"
```

---

## Supplementary Figures

### Figure S1. Database Coverage Analysis

Same as v1.0

---

### Figure S2. Two-Track Framework Detailed Schematic (NEW)

Full-page diagram with detailed workflow for both tracks.

---

### Figure S3. Track A: Complete GNN Analysis Details

**A.** Graph construction pipeline
**B.** Node feature initialization
**C.** Message passing illustration
**D.** Link prediction training

---

### Figure S4. Track A: Pathway Activation Network

Network visualization showing all pathway enzymes with expression changes overlaid.

---

### Figure S5. Track B: Docking Methodology (NEW)

**A.** Protein structure preparation workflow
**B.** Ligand conformer generation
**C.** Blind docking grid setup
**D.** Score distribution

---

### Figure S6. Methodological Clarification Diagram (NEW)

Detailed explanation of why GNN predictions and docking address different questions.

---

### Figure S7. Quality Control Metrics

Same as v1.0 (Figure S5)

---

## Figure Production Notes

### Color Scheme

| Element | Color | Hex |
|---------|-------|-----|
| Track A elements | Blue | #2196F3 |
| Track B elements | Orange | #FF9800 |
| Upregulated | Red | #F44336 |
| Downregulated | Blue | #2196F3 |
| Metabolites | Green | #4CAF50 |
| Enzymes | Purple | #9C27B0 |

### Software Requirements

- Network visualization: Cytoscape 3.9+
- Docking visualization: PyMOL 2.5+ or ChimeraX
- Statistical plots: Python (matplotlib/seaborn) or R (ggplot2)
- Pathway diagrams: BioRender or Illustrator

### Resolution Requirements

- Main figures: 300 dpi minimum
- Supplementary: 150 dpi minimum
- Vector format preferred for diagrams

---

## Figure Checklist

- [ ] Figure 1: Study overview and two-track framework
- [ ] Figure 2: Metabolomics volcano plot and heatmap
- [ ] Figure 3: Pathway enrichment results
- [ ] Figure 4: GNN architecture and performance
- [ ] Figure 5: Proteomics validation (Track A completion)
- [ ] Figure 6: Docking predictions (Track B, with warnings)
- [ ] Figure 7: Track relationship summary
- [ ] Figure S1: Database coverage (from v1)
- [ ] Figure S2: Two-track detailed schematic
- [ ] Figure S3: GNN analysis details
- [ ] Figure S4: Pathway network
- [ ] Figure S5: Docking methodology
- [ ] Figure S6: Methodological clarification
- [ ] Figure S7: Quality control (from v1)

---

*Figure descriptions v2.0 complete: 2026-01-21*
