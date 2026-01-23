# Public Soybean Ethylene Time-Series Datasets

## Overview

This document lists publicly available soybean (Glycine max) transcriptome datasets
that could be used to supplement the current static (72h single time point) analysis
with temporal information.

---

## Most Relevant: Ethylene Treatment Time-Course

### Soybean Leaf Abscission Study (2020)

**Dataset**: RNA-seq of leaf abscission zones (LAZ) and petioles after ethylene treatment

**Time points**: 0, 12, 24, 48, and 72 hours

**Key findings**:
- 5,206 soybean TF genes analyzed
- 1,088 TFs differentially regulated (>8-fold) in LAZ over time
- Direct ethylene treatment time-course data

**Relevance**: **HIGH** - Direct ethylene response time-course, could validate temporal ordering of TF activity

**Reference**: Check GEO/SRA for "soybean ethylene abscission RNA-seq"

---

## Recent 2025 Dataset

### Temporal Shifts in Hormone Signaling Networks (2025)

**Source**: [MDPI International Journal of Molecular Sciences](https://www.mdpi.com/1422-0067/26/13/6455)

**Description**: Strand-specific, high-depth temporal transcriptome atlas of soybean inflorescences

**Time points**: Stadium 0 through Stadium 3 (developmental stages)

**Hormone pathways covered**:
- Auxin
- Cytokinin
- Gibberellin
- Abscisic acid
- **Ethylene**
- Jasmonate
- Salicylate

**Reference genome**: G. max GCF_000004515.6 (NCBI, November 2024)

**Relevance**: **MEDIUM** - Contains ethylene signaling data but developmental context differs

---

## Comprehensive Resources

### Soybean Expression Atlas v2

**Source**: [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.04.28.538661v1.full)

**Content**: 5,481 publicly available RNA-seq samples

**Features**:
- Transcript-level and gene-level abundance matrices
- Multiple conditions and treatments
- Searchable by gene/condition

**Relevance**: **HIGH** - May contain ethylene treatment samples for meta-analysis

### SoyBase Collections

**Source**: [SoyBase Data Collections](https://www.soybase.org/collections/)

**Content**: Curated soybean genomic and transcriptomic datasets

**Relevance**: **MEDIUM** - Comprehensive resource for finding additional datasets

---

## Recommended Integration Strategy

### For Current Project

1. **Primary target**: Soybean ethylene abscission time-course data (0-72h)
   - Direct temporal validation of TF activity ordering
   - Same ethylene treatment context

2. **Secondary**: Soybean Expression Atlas v2
   - Meta-analysis across conditions
   - Validation of TF-enzyme co-expression

3. **Tertiary**: 2025 hormone signaling dataset
   - Ethylene pathway temporal dynamics
   - Cross-reference with predicted TF targets

### Integration Approach

```python
# Pseudo-temporal ordering using external time-course data
# 1. Map genes from current data to time-course expression
# 2. Assign pseudo-time based on peak expression
# 3. Validate HGT predictions against temporal ordering
```

---

## Search Queries for Additional Data

### GEO Database
```
"Glycine max" AND "ethylene" AND "RNA-seq" AND "time course"
"soybean" AND "ethylene treatment" AND "transcriptome"
```

### SRA Database
```
(Glycine max[Organism]) AND ethylene AND RNA-Seq[Strategy]
```

### Literature Search
```
soybean isoflavonoid ethylene time-course transcriptome
Glycine max phenylpropanoid ethylene temporal gene expression
```

---

## Notes

- Current project data: Static comparison at 72h
- Time-course data would enable:
  - Temporal ordering validation
  - Granger causality analysis
  - Dynamic regulatory network construction
  - Stronger mechanism claims

---

*Document created for: GNN vs Transformer evaluation project*
*Purpose: Identify time-series data to supplement static analysis*
