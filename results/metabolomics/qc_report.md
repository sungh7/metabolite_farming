# Metabolomics QC Report (MTBLS531)

**Dataset**: MTBLS531 (MetaboLights)
**Organism**: *Glycine max* (Soybean)
**Platform**: LC-MS Untargeted Metabolomics
**Conditions**: Ethylene (ET) vs Control (Ctrl)

---

## 1. Sample Overview

| Group | N (Samples) | Notes |
| :--- | :--- | :--- |
| Control | 3 | Untreated soybean leaves |
| Ethylene | 3 | Ethylene-treated soybean leaves |
| **Total** | **6** | Minimum n=3 per group (biological replicates) |

---

## 2. QC Metrics Summary

### 2.1 Detection Rate
*   **Feature Detection Threshold**: Features detected in >= 70% of samples per group.
*   **Result**: 150+ features passed detection threshold (exact count varies by preprocessing).

### 2.2 Blank Contamination
*   **Status**: ⚠️ **Not Explicitly Provided** in MTBLS531 metadata.
*   **Mitigation**: Features with very low intensities (< 1000 raw counts) were excluded as potential noise.
*   **Recommendation**: Future studies should include explicit blank injections for carryover assessment.

### 2.3 Pooled QC RSD (Coefficient of Variation)
*   **Status**: ⚠️ **Not Explicitly Provided** in MTBLS531 metadata.
*   **Assumption**: Standard LC-MS QC protocols assume RSD < 30% for technical replicates. MTBLS531 is a published study, so QC was likely performed by original authors.
*   **Recommendation**: When reusing public data, explicitly document assumed QC status.

### 2.4 Injection Order Drift
*   **Status**: ⚠️ **Not Assessed** (requires raw injection order metadata).
*   **Mitigation**: Log2 transformation and median normalization were applied to reduce intensity drift effects.

---

## 3. Normalization & Transformation

| Step | Method | Justification |
| :--- | :--- | :--- |
| **Log Transformation** | Log2 | Stabilize variance, approximate normality. |
| **Normalization** | Median Centering | Correct for sample-to-sample intensity variation. |
| **Missing Value Handling** | Min/2 Imputation | Replace missing values with half of the minimum detected value per feature. |

---

## 4. Differential Analysis Summary

| Metric | Value | Notes |
| :--- | :--- | :--- |
| **Fold Change Threshold** | Log2FC >= 1.0 | 2-fold change minimum. |
| **Significance Threshold** | FDR <= 0.05 | Benjamini-Hochberg correction. |
| **Upregulated Features (ET vs Ctrl)** | 7 | See `top_features_up_with_msi.tsv`. |

---

## 5. MSI Level Distribution

| MSI Level | Count | Description |
| :--- | :--- | :--- |
| **Level 2** (Putatively Annotated) | 2 | ChEBI + Name + KEGG ID. |
| **Level 3** (Putatively Characterized) | 5 | ChEBI + Name only. |
| **Level 4** (Unknown) | 0 | Not applicable (all top features have ChEBI). |

---

## 6. Limitations & Recommendations

1.  **Blank/QC Data Unavailable**: MTBLS531 does not include explicit blank or pooled QC injections in the public metadata. Future analyses should prioritize datasets with complete QC documentation.
2.  **MSI Level 1 Not Achievable**: No authentic standards were available for RT + MS/MS matching. All annotations are database-based (Level 2-3).
3.  **Time-Course Not Available**: MTBLS531 is a single-timepoint comparison (ET vs Ctrl). Time-course dynamics cannot be assessed.

---

**Report Generated**: 2026-01-02
**Pipeline**: NoTaMe-inspired (Feature Detection -> QC Filtering -> Normalization -> Differential)
