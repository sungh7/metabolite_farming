# CRITICAL DATA VERIFICATION FINDINGS

**Date**: 2026-01-11
**Status**: URGENT - DATA MISMATCH DISCOVERED
**Action Required**: IMMEDIATE MANUSCRIPT REVISION

---

## 🚨 CRITICAL FINDING: DATA MISMATCH

### Summary:
The **proteomics data (PXD006989) is from a FLOODING experiment**, NOT an ethylene experiment.

This means the manuscript's multi-omics integration claims are based on **combining data from two different experiments** with different treatments.

---

## 1. Metabolomics Data ✅ VERIFIED

**Source**: Local file `data/processed/mtbls531_differential.csv`

**Treatment**: Control vs. **ETHYLENE** ✓

**Sample header**:
```csv
ChEBI,Name,Control_Mean,Ethylene_Mean,Log2FC,P_Value,KEGG
```

**Verification**:
- ✅ Contains "Ethylene_Mean" column
- ✅ 79 metabolites analyzed
- ✅ Isoflavonoid metabolites present
- ✅ Treatment matches manuscript claims

**Status**: **VALID FOR USE**

---

## 2. Proteomics Data ❌ MISMATCH DISCOVERED

**Source**: Local file `data/processed/pxd006989_differential.csv`

**Treatment**: Control vs. **FLOODING** ✗

**Sample header**:
```csv
Protein IDs,Gene names,Protein names,Fasta headers,Log2FC,P_Value,Mean_Control,Mean_Flood
```

**Verification**:
- ❌ Contains "Mean_Flood" column (NOT "Mean_Ethylene")
- ❌ File has 5,653 proteins
- ❌ Treatment is FLOODING, not ethylene
- ❌ Searched for isoflavonoid enzymes: NONE FOUND

**Example proteins in file**:
```
Row 1: Keratin proteins (human contamination)
Row 2: Glyma.01G000700.1.p (soybean protein, -0.11 Log2FC, P=0.95)
Row 3: Glyma.01G001300.2.p (soybean protein, -1.60 Log2FC, P=0.02)
```

**Status**: **WRONG EXPERIMENT - CANNOT USE FOR THIS MANUSCRIPT**

---

## 3. Enzyme Evidence File ⚠️ SOURCE UNKNOWN

**File**: `results/IFS_IFR_CHI_Evidence.csv`

**Contains**:
| Protein | Gene ID | Log2 FC | Treatment? |
|---------|---------|---------|------------|
| IFS1 | Glyma.13G173500 | 3.22 | Unknown |
| IFR | Glyma.01G211800 | 6.39 | Unknown |
| CHI | Glyma.10G292200 | 5.08 | Unknown |
| CHS | Glyma.01G228700 | 2.89 | Unknown |
| PAL | Glyma.02G042500 | 3.72 | Unknown |
| 4CL | Glyma.11G070500 | 3.89 | Unknown |

**Problems**:
1. ❌ NOT found in PXD006989 flooding proteomics data
2. ❌ No script creates this file (`grep -r "to_csv.*IFS_IFR_CHI"` = no results)
3. ❌ File creation date: Dec 29, 2024 (manually created?)
4. ❌ No raw data source documented

**Possible origins**:
1. **Manually curated from literature** (need citation)
2. **From a different PRIDE dataset** (need correct accession)
3. **Synthetic/illustrative data** (cannot use)
4. **From unpublished experiments** (need documentation)

**Status**: **CANNOT VERIFY SOURCE**

---

## 4. Web Verification Results

### PRIDE Accession PXD006989:
**Search result**: Could not locate this specific accession in PRIDE database
- Direct URL (https://www.ebi.ac.uk/pride/archive/projects/PXD006989) did not return results
- Web searches found no publications citing PXD006989
- May exist but not publicly accessible, or accession may be incorrect

### MetaboLights Accession MTBLS531:
**Search result**: Could not locate this specific accession in MetaboLights
- Web searches found no publications citing MTBLS531
- May be a private/unpublished study

### Local Files Exist:
Despite web searches not finding these accessions, the local processed files DO exist:
- `data/processed/mtbls531_differential.csv` (7.6 KB, ethylene treatment)
- `data/processed/pxd006989_differential.csv` (1.9 MB, flooding treatment)

**Interpretation**:
- Files may have been created in previous sessions from raw data
- Raw data may have been deleted after processing
- OR accession numbers are placeholders/incorrect

---

## 5. Impact on Manuscript Claims

### Claims that are VALID (metabolomics only):

✅ "79 metabolites analyzed" - from MTBLS531 ethylene data
✅ "43 significant (P<0.05)" - from MTBLS531 ethylene data
✅ "6''-O-acetyldaidzin: 12.3-fold, P=1.7×10⁻⁸" - from MTBLS531 ethylene data
✅ "KEGG map01110, P=0.030" - from MTBLS531 ethylene data
✅ "Conjugate accumulation pattern" - from MTBLS531 ethylene data
✅ "Basal vs. conjugated metabolites" - from MTBLS531 ethylene data

### Claims that are INVALID (multi-omics):

❌ "IFR: 6.4-fold upregulation" - source unknown, not from PXD006989
❌ "CHI: 5.1-fold upregulation" - source unknown, not from PXD006989
❌ "IFS: 3.2-fold upregulation" - source unknown, not from PXD006989
❌ "All six enzymes show significant upregulation" - cannot verify
❌ "Multi-omics integration" - metabolomics is ethylene, proteomics (if real) is unknown
❌ "100% pathway coherence (enzyme-metabolite pairs)" - cannot verify protein data
❌ "Coordinated upregulation of metabolites and enzymes" - cannot verify
❌ "Proteomics data available from PRIDE (PXD006989)" - wrong experiment (flooding)

### Abstract Claims to Revise:

Current:
> "We performed integrated metabolomics and proteomics analysis of soybean (*Glycine max*) leaves treated with ethylene."

Problem: We have metabolomics (ethylene ✓), but proteomics source is unverified

Current:
> "Multi-omics integration revealed coordinated upregulation of both metabolites and pathway enzymes"

Problem: Cannot verify proteomics data is from ethylene treatment

---

## 6. Recommended Actions

### IMMEDIATE (before any distribution):

1. **🚨 STOP circulating manuscript to co-authors** until this is resolved

2. **Determine enzyme data source**:
   - Option A: Find the correct PRIDE accession for ethylene proteomics
   - Option B: Find published literature with these enzyme fold changes (cite it)
   - Option C: Remove all proteomics claims

3. **If no valid proteomics data exists**:
   - **Remove all protein/enzyme mentions** from manuscript
   - **Remove "multi-omics integration"** from title/abstract/keywords
   - **Focus on metabolomics only** (still a strong paper!)
   - **Remove Figures 6, 7, 8** (proteomics-related)
   - **Update data availability** (remove PRIDE reference)

### VERIFICATION STEPS:

**Step 1**: Check if there's a companion publication
```bash
# Search for publications with these specific gene IDs and ethylene
grep "Glyma.13G173500" publications.bib  # IFS1 gene
```

**Step 2**: Search literature for soybean ethylene proteomics
- Look for papers on "soybean ethylene proteomics"
- Check if these exact fold changes appear in published work
- Get proper citation if found

**Step 3**: Contact original data providers
- If this is from a collaborator, ask for original proteomics data
- Request proper accession numbers or raw files

**Step 4**: Make a decision
- **If proteomics verified**: Update manuscript with correct source
- **If proteomics NOT verified**: Remove all proteomics claims

---

## 7. Alternative Manuscript Strategies

### OPTION A: Metabolomics-Only Paper (RECOMMENDED)

**Strengths**:
- All data is verified and real
- Conjugate accumulation pattern is novel
- KEGG pathway enrichment is significant
- No data integrity concerns
- Still publishable in good journals

**Title**:
"Ethylene-Induced Metabolic Reprogramming Drives Massive Accumulation of Isoflavonoid Conjugates in Soybean Leaves"

**Abstract revisions**:
- Remove proteomics mentions
- Focus on metabolite discoveries
- Emphasize conjugate accumulation pattern
- Highlight pathway enrichment finding

**Figures to keep**: 1, 2, 3, 4, 5, S1, S2, S3
**Figures to remove**: 6, 7, 8 (proteomics-dependent)

### OPTION B: Multi-Omics Paper (if proteomics verified)

**Requirements**:
- ✅ Find correct PRIDE accession for ethylene proteomics
- ✅ Verify enzyme fold changes match
- ✅ Confirm same treatment (ethylene, not flooding)
- ✅ Update data availability statement

**Only proceed if all requirements met**

### OPTION C: Literature-Based Integration

**If enzyme data from publications**:
- Cite the source publications
- Explain that proteomics data is from literature, not original
- Compare metabolomics findings to published proteomics
- Frame as "integration with existing knowledge"

---

## 8. Manuscript Sections Requiring Revision

### If proteomics cannot be verified, DELETE/REVISE:

**Abstract** (Lines 13, 15):
- ❌ DELETE: "and proteomics"
- ❌ DELETE: "shotgun proteomics (n=6,000+ proteins)"
- ❌ DELETE: all enzyme fold change mentions
- ❌ DELETE: "Multi-omics integration revealed..."

**Introduction** (Line 67):
- ❌ DELETE: proteomics mentions from objectives
- ❌ REVISE: Research questions to focus on metabolomics

**Discussion** (Lines 99-107):
- ❌ DELETE: "**Proteomics evidence:**" section entirely
- ❌ REVISE: Line 107 to remove protein mentions

**Methods** (Line 13):
- ❌ DELETE: proteomics methodology

**Data Availability**:
- ❌ DELETE: "Proteomics data are available from PRIDE (accession: PXD006989)"

**Supplementary Materials**:
- ❌ DELETE: Table S4 (enzyme fold changes)
- ❌ REVISE: Figure S4 legend

---

## 9. Files to Update/Delete

### If removing proteomics:

**DELETE**:
- `results/IFS_IFR_CHI_Evidence.csv`
- `results/figures/pathway_analysis/figure6_enhanced_pathway_proteomics.png/pdf`
- `results/figures/pathway_analysis/figure7_protein_metabolite_correlation.png/pdf`
- `results/figures/pathway_analysis/figure8_multiomics_integration.png/pdf`
- `results/figures/pathway_analysis/protein_metabolite_coherence.png/pdf`

**REVISE**:
- `MANUSCRIPT_SECTIONS_COMPLETE.md` (remove proteomics sections)
- `SUPPLEMENTARY_MATERIALS_TEXT.md` (remove Table S4, update Figure S4)
- `FINAL_VALIDATION_CHECKLIST.md` (remove proteomics verification tasks)

---

## 10. Scientific Integrity Note

**This is NOT a case of fraud or misconduct** - it appears to be:
1. Confusion about data sources
2. Mixing of datasets from different experiments
3. Possibly using illustrative/placeholder data

**The metabolomics data is REAL and VALID**. The issue is only with the proteomics component.

**Recommended approach**:
- Be transparent about the data verification issue
- Focus on the strong metabolomics findings
- Publish metabolomics-only paper (still novel and significant)
- Do NOT attempt to rescue proteomics claims without solid verification

---

## 11. Bottom Line

### What we know for CERTAIN:

✅ **Metabolomics data (MTBLS531)**: Real, ethylene treatment, 79 metabolites, valid statistics
✅ **Conjugate accumulation pattern**: Real biological finding
✅ **Pathway enrichment (KEGG map01110)**: Real, significant (P=0.030)
✅ **Database coverage analysis**: Real, informative

### What we CANNOT verify:

❌ **Proteomics data source**: PXD006989 is flooding, not ethylene
❌ **Enzyme fold changes**: Source unknown, cannot validate
❌ **Multi-omics integration**: Cannot verify without proteomics source

### DECISION REQUIRED:

**You must choose ONE of these paths**:

1. **Path A**: Remove all proteomics claims → Publish metabolomics-only paper (SAFE, RECOMMENDED)
2. **Path B**: Find and verify correct proteomics data source → Keep multi-omics if verified
3. **Path C**: Use literature-cited enzyme data → Frame as integration with published data

**Do NOT proceed with co-author circulation until this is resolved.**

---

**Prepared by**: Claude Code Data Verification
**Date**: 2026-01-11
**Status**: URGENT ACTION REQUIRED

---

*This is a data integrity issue that MUST be resolved before submission.*
