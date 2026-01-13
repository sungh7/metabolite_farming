# Integrated Multi-Omics Analysis of Ethylene-Induced Isoflavonoid Biosynthesis in Soybean

Multi-omics integration analysis revealing coordinated activation of isoflavonoid biosynthesis pathway in ethylene-treated soybean leaves.

## Overview

This repository contains the complete analysis pipeline and publication materials for integrating metabolomics (MetaboLights MTBLS531) and proteomics (PRIDE PXD006989) data to characterize ethylene-induced metabolic reprogramming in *Glycine max*.

## Key Findings

- **Pathway Activation**: Ethylene treatment significantly induces secondary metabolite biosynthesis (KEGG map01110, P=0.030)
- **Conjugate-Selective Accumulation**: Malonylated/acetylated isoflavonoid conjugates increase >4,000-fold while basal aglycones show minimal change (~1.1-fold)
- **Multi-Omics Concordance**: 100% directional agreement between 6 biosynthetic enzymes (IFS, IFR, CHI, CHS, PAL, 4CL) and corresponding metabolites
- **Cross-Database Validation**: Both KEGG and PlantCyc independently identify isoflavonoid biosynthesis as top-ranked pathway

## Data Sources

| Dataset | Accession | Description |
|---------|-----------|-------------|
| Metabolomics | [MTBLS531](https://www.ebi.ac.uk/metabolights/MTBLS531) | 79 metabolites, soybean leaf ethylene response |
| Proteomics | [PXD006989](https://www.ebi.ac.uk/pride/archive/projects/PXD006989) | Leaf proteome, ethylene/ABA treatments |

## Repository Structure

```
ethylene/
├── data/
│   ├── experimental/          # Raw data from public databases
│   └── processed/             # Processed differential expression files
├── src/
│   ├── process_proteomics.py          # Proteomics processing pipeline
│   ├── generate_pathway_figures.py    # Main figure generation
│   ├── generate_enhanced_figures.py   # Multi-omics integration figures
│   ├── statistical_analysis.py        # Statistical analysis
│   └── plantcyc_pathway_enrichment.py # PlantCyc enrichment analysis
├── results/
│   ├── figures/pathway_analysis/      # Publication-ready figures (PNG/PDF)
│   └── *.csv                          # Statistical results tables
├── MANUSCRIPT_SECTIONS_COMPLETE.md    # Full manuscript text
├── SUPPLEMENTARY_MATERIALS_TEXT.md    # Supplementary methods and tables
├── COVER_LETTER.md                    # Journal submission cover letter
└── FINAL_VALIDATION_CHECKLIST.md      # Data verification checklist
```

## Requirements

```
python >= 3.10
pandas >= 2.0
numpy >= 1.24
scipy >= 1.10
matplotlib >= 3.7
seaborn >= 0.12
```

Install dependencies:
```bash
pip install -r requirements.txt
```

## Reproducing the Analysis

### 1. Process Proteomics Data
```bash
python src/process_proteomics.py
```

### 2. Generate Pathway Figures
```bash
python src/generate_pathway_figures.py
python src/generate_enhanced_figures.py
python src/generate_supplementary_materials.py
```

### 3. Run Statistical Analysis
```bash
python src/statistical_analysis.py
```

## Key Results

### Enzyme Fold Changes (Ethylene vs Control)

| Enzyme | Gene ID | Fold Change | P-Value |
|--------|---------|-------------|---------|
| IFR | Glyma.01G211800 | 6.39 | 0.037 |
| CHI | Glyma.10G292200 | 5.08 | 0.047 |
| 4CL | Glyma.11G070500 | 3.89 | 0.005 |
| PAL | Glyma.02G042500 | 3.72 | 0.007 |
| IFS1 | Glyma.13G173500 | 3.22 | 0.006 |
| CHS | Glyma.01G228700 | 2.89 | 0.007 |

### Top Metabolite Changes

| Metabolite | Log2FC | P-Value |
|------------|--------|---------|
| 6''-O-Acetyldaidzin | 12.30 | 1.7e-8 |
| 6''-O-Malonylgenistin | 12.13 | 8.3e-8 |
| 6''-O-Malonyldaidzin | 11.90 | 3.2e-8 |

## Figures

All publication figures available in `results/figures/pathway_analysis/`:

- **Figure 1**: KEGG pathway enrichment analysis
- **Figure 2**: Metabolite volcano plot
- **Figure 5**: Isoflavonoid biosynthesis pathway
- **Figure 6**: Enhanced pathway with proteomics integration
- **Figure 8**: Multi-omics integration summary
- **Figure 10**: Graphical abstract

## Citation

If you use this analysis or code, please cite:

> [Authors]. (2026). Integrated Multi-Omics Analysis Reveals Coordinated Activation of Isoflavonoid Biosynthesis in Ethylene-Treated Soybean Leaves. *[Journal]*.

## License

This project is available for academic and research purposes. Please cite the original data sources (MTBLS531, PXD006989) when using the data.

## Contact

For questions about this analysis, please open an issue or contact the corresponding author.

---

*Analysis performed: January 2026*
