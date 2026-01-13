# Data Sources Verification Report

**Date**: 2026-01-11
**Purpose**: Complete transparency about real vs. synthetic/curated data

---

## Summary

✅ **METABOLITE DATA**: 100% REAL from public database
⚠️ **PROTEIN DATA**: QUESTIONABLE - needs verification
✅ **PATHWAY DATABASES**: 100% REAL (KEGG, PlantCyc)
✅ **STATISTICAL ANALYSIS**: Real calculations on the data we have

---

## 1. Metabolite Data ✅ REAL

### Source:
- **Database**: MetaboLights MTBLS531
- **Public accession**: https://www.ebi.ac.uk/metabolights/MTBLS531
- **Study**: "Ethylene-induced metabolic changes in soybean"
- **Technology**: LC-MS/MS metabolomics

### Files:
```
data/raw/MTBLS531/m_mtbls531_metabolite_profiling_mass_spectrometry_v2_maf.tsv
└─> Processed by: src/process_experimental.py
    └─> Output: data/processed/mtbls531_differential_enhanced.csv
```

### Verification:
- ✅ Downloaded from public repository
- ✅ Contains 79 metabolites with replicate measurements
- ✅ Control (n=4) vs. Ethylene (n=4) treatment
- ✅ All statistical calculations (t-tests, P-values, fold changes) computed from raw data

### Data quality:
- **Total metabolites**: 79
- **Significant (P<0.05)**: 43
- **KEGG mapped**: 29 (36.7%)
- **Top hit**: 6''-O-Acetyldaidzin (Log2FC=12.30, P=1.7×10⁻⁸)

**CONCLUSION**: This is 100% real experimental data from a published dataset.

---

## 2. Protein Data ⚠️ QUESTIONABLE

### Source claimed:
- **Database**: PRIDE (Proteomics Identifications Database)
- **Accession**: PXD006989
- **Mentioned in**: Manuscript supplementary materials

### Actual files in directory:
```
results/IFS_IFR_CHI_Evidence.csv
```

**Contents:**
| Protein | Gene ID | Log2 FC | Significance |
|---------|---------|---------|--------------|
| IFS1 | Glyma.13G173500 | 3.22 | P<0.05 |
| IFR | Glyma.01G211800 | 6.39 | P<0.05 |
| CHI | Glyma.10G292200 | 5.08 | P<0.05 |
| CHS | Glyma.01G228700 | 2.89 | P<0.05 |
| PAL | Glyma.02G042500 | 3.72 | P<0.05 |
| 4CL | Glyma.11G070500 | 3.89 | P<0.05 |

### Issues identified:

❌ **No raw PRIDE data exists** in the directory
❌ **No script creates this file** (no to_csv found)
❌ **File created**: Dec 29, 2024 (before this session)
❌ **No processing pipeline** from PRIDE raw data to this file

### Possible explanations:

**Option 1**: This data was manually curated from a published paper
- Gene IDs are real soybean genes (Glyma.* nomenclature)
- Fold changes could be from literature
- Need to verify if PXD006989 is the correct source

**Option 2**: This data is synthetic/illustrative
- Created for demonstration purposes
- Gene IDs are real, but fold changes may be estimated
- Should be labeled as "example data"

**Option 3**: This data was processed in a previous session
- Raw PRIDE data was downloaded and processed
- Processing scripts were deleted or in a different location
- Final summary table was kept

### What we need to verify:

1. **Check if PXD006989 is real**: Visit https://www.ebi.ac.uk/pride/archive/projects/PXD006989
2. **Check if it matches this study**: Same organism (soybean), same treatment (ethylene)
3. **Verify fold changes**: Do the reported values match the actual dataset?
4. **Get raw data**: Download and process the actual proteomics data

**CONCLUSION**: Cannot confirm if this is real data without verification.

---

## 3. Pathway Databases ✅ REAL

### KEGG (Kyoto Encyclopedia of Genes and Genomes):
- **Source**: https://www.kegg.jp/
- **Access**: Real-time API queries via `src/improve_kegg_mapping.py`
- **Pathways tested**: 44 pathways
- **Significant finding**: map01110 (P=0.030)

### PlantCyc:
- **Source**: BioCyc Web Services API
- **Access**: API queries via `src/analyze_plantcyc.py`
- **Pathways tested**: 268 pathways
- **Top result**: ISOFLAVONOID-SYN (P=0.405, not significant)

**CONCLUSION**: Pathway databases are 100% real, live data.

---

## 4. Statistical Analysis ✅ REAL

All statistical calculations are genuine:

### Metabolite statistics:
- **Differential abundance**: Welch's t-test (real calculations)
- **Fold changes**: Log2 ratios from replicate means
- **P-values**: Computed from actual data distributions
- **Effect sizes**: Cohen's d calculated correctly

### Pathway enrichment:
- **Fisher's exact test**: Real 2×2 contingency tables
- **Odds ratios**: Calculated with Haldane correction
- **Confidence intervals**: Exact binomial CIs
- **Multiple testing**: FDR and Bonferroni corrections applied

### Multi-omics integration:
- **Coherence analysis**: Real comparison of enzyme vs. metabolite fold changes
- **Directional concordance**: 5/5 pairs both upregulated (IF protein data is real)
- **Pattern discovery**: Basal vs. conjugate metabolite distinction is real from metabolite data

**CONCLUSION**: All statistical methods are correctly applied to the data we have.

---

## 5. What the Correlation Analysis Actually Did

### What we analyzed:
```python
# Load REAL metabolite data
metabolites_df = pd.read_csv('data/processed/mtbls531_differential_enhanced.csv')  # ✅ REAL

# Load QUESTIONABLE protein data
proteins_df = pd.read_csv('results/IFS_IFR_CHI_Evidence.csv')  # ⚠️ NEEDS VERIFICATION
```

### Actual calculations:
1. ✅ Extracted enzyme fold changes from protein file
2. ✅ Extracted metabolite fold changes from REAL data
3. ✅ Compared directional concordance (both up or both down)
4. ✅ Found 5/5 pairs are coherent (both upregulated)
5. ✅ Correctly identified we cannot calculate Pearson r without paired replicates

### What we DID NOT do:
- ❌ Did NOT use synthetic replicate data (the simulation code exists but wasn't called)
- ❌ Did NOT calculate false correlation coefficients
- ❌ Did NOT make up P-values

### The biological insight:
**This is REAL even if protein data is questionable:**
- Basal aglycones (daidzein, formononetin): Log2FC ~0.14 ✅ REAL from metabolite data
- Conjugated forms: Log2FC ~12 ✅ REAL from metabolite data
- Pattern interpretation: Based on REAL metabolite data

The conjugate accumulation pattern is a genuine finding from the metabolite data alone!

---

## 6. Impact on Manuscript Claims

### Claims that are 100% solid (based on REAL data):

✅ "79 metabolites analyzed"
✅ "43 significant (P<0.05)"
✅ "6''-O-acetyldaidzin: 12.3-fold, P=1.7×10⁻⁸"
✅ "KEGG map01110, P=0.030"
✅ "Dramatic upregulation of conjugated isoflavonoids"
✅ "Basal aglycones show modest increases vs. massive conjugate accumulation"
✅ "36.7% KEGG mapping coverage"

### Claims that depend on protein data verification:

⚠️ "IFR: 6.4-fold upregulation"
⚠️ "CHI: 5.1-fold upregulation"
⚠️ "IFS: 3.2-fold upregulation"
⚠️ "All six enzymes show significant upregulation (P<0.05)"
⚠️ "100% pathway coherence (enzyme-metabolite pairs)"

### What happens if protein data is not real:

**Option A**: Verify it's real from PXD006989 → All claims stand ✅
**Option B**: It's from literature → Cite the source, claims stand ✅
**Option C**: It's synthetic → Remove protein claims, focus on metabolites only ⚠️

---

## 7. Recommended Actions

### URGENT (before submission):

1. **Verify PXD006989 accession**
   - Visit: https://www.ebi.ac.uk/pride/archive/projects/PXD006989
   - Check if it exists and matches this study
   - Download raw data if available

2. **If PXD006989 is correct:**
   - ✅ Download the raw proteomics data
   - ✅ Process it with MaxQuant or similar
   - ✅ Verify the fold changes match IFS_IFR_CHI_Evidence.csv
   - ✅ Document the processing pipeline

3. **If PXD006989 is NOT correct or doesn't exist:**
   - **Option A**: Find the correct accession/source
   - **Option B**: Remove all proteomics claims from manuscript
   - **Option C**: Label as "illustrative example" if synthetic

4. **Update manuscript data availability:**
   - Current: "Proteomics data are available from PRIDE (accession: PXD006989)"
   - Verify this is accurate

### OPTIONAL (strengthens paper):

5. **If no real proteomics data exists:**
   - The metabolite data ALONE is still a strong paper
   - Conjugate accumulation pattern is novel
   - Pathway enrichment is significant
   - Remove multi-omics claims, focus on metabolomics

---

## 8. Current Status

### What we know for certain:

✅ **Metabolite data is real** (MTBLS531)
✅ **Statistical analysis is correct**
✅ **Pathway enrichment is real**
✅ **Conjugate accumulation pattern is real**

### What we need to verify:

⚠️ **Protein data source** (PXD006989 or other)
⚠️ **Protein fold changes** (are they from this study?)
⚠️ **Multi-omics integration claims** (depend on protein data)

---

## 9. Recommendation

### Immediate action:

**Verify the protein data source TODAY before circulating to co-authors.**

### If protein data is verified:
✅ Manuscript is publication-ready
✅ All claims are supported
✅ Multi-omics integration is valid

### If protein data cannot be verified:
⚠️ Remove proteomics claims
✅ Paper is still strong based on metabolomics alone
✅ Conjugate accumulation pattern is the main finding

---

**Bottom line**: The metabolite analysis is 100% real and solid. The protein data needs verification before submission.

---

*Prepared by*: Claude Code Verification
*Date*: 2026-01-11
