# 🧬 이소플라보노이드 생합성 경로 예측 종합 보고서

**생성일**: 2026년 1월 22일  
**분석**: GNN 기반 Enzyme-Metabolite Link Prediction  

---

## 📊 예측 완료된 대사체

| Target Metabolite | 경로 | 출력 디렉토리 |
|-------------------|------|---------------|
| **Daidzein** | Isoflavonoid Biosynthesis | `results/isoflavonoid_prediction/` |
| **Genistein** | Isoflavonoid Biosynthesis | `results/genistein_prediction/` |
| **Formononetin** | Isoflavonoid Biosynthesis | `results/formononetin_prediction/` |
| **Glycitein** | Isoflavonoid Biosynthesis | `results/glycitein_prediction/` |

---

## 🔬 생합성 경로 개요

```
                    ISOFLAVONOID BIOSYNTHESIS PATHWAY
═══════════════════════════════════════════════════════════════════════

    L-Phenylalanine
           │
           ▼  [PAL: Phenylalanine ammonia-lyase]
    trans-Cinnamic acid
           │
           ▼  [C4H: Cinnamate 4-hydroxylase]
    p-Coumaric acid
           │
           ▼  [4CL: 4-Coumarate:CoA ligase]
    p-Coumaroyl-CoA  +  3 × Malonyl-CoA
           │
           ▼  [CHS: Chalcone synthase]
    Naringenin chalcone
           │
           ▼  [CHI: Chalcone isomerase]  ★★ (5.1× 증가, p<0.05)
    Naringenin / Liquiritigenin
           │
           ▼  [IFS: Isoflavone synthase]  ★★ (3.2× 증가, p<0.05)
    2-Hydroxyisoflavanone
           │
           ▼  [HID: 2-Hydroxyisoflavanone dehydratase]
           │
     ┌─────┴─────┬─────────────┐
     ▼           ▼             ▼
 DAIDZEIN    GENISTEIN    GLYCITEIN
     │           │             │
     ▼           ▼             ▼
 Formononetin  Biochanin A  (derivatives)

═══════════════════════════════════════════════════════════════════════
```

---

## 📋 단계별 Top-1 예측 효소

| Step | Enzyme Class | Top-1 Candidate | UniProt | GNN Score | Combined Score | Log2FC | P-value | Significant |
|------|--------------|-----------------|---------|-----------|----------------|--------|---------|-------------|
| 1 | PAL | I1NHH9_SOYBN | I1NHH9 | 0.7905 | 0.7905 | - | - | - |
| 2 | C4H | (No candidates) | - | - | - | - | - | - |
| 3 | 4CL | 4CL | H2BER4 | 0.7093 | 0.7093 | -1.41 | 0.109 | No |
| 4 | CHS | PROPEP914 | K7LSB9 | 0.7089 | 0.7089 | - | - | - |
| 5 | **CHI** | **CHI1B2-2** | **I1LFF1** | 0.7224 | **1.7377** | **+5.08** | **0.047** | **★ Yes** |
| 6 | **IFS** | **ifs1** | **Q9M6D6** | 0.7186 | **1.3632** | **+3.22** | **0.006** | **★ Yes** |
| 7 | HID | HIDH | Q5NUF3 | 0.7848 | 0.7848 | - | - | - |

---

## 🎯 주요 발견

### 1. 에틸렌 반응 효소 (Ethylene-responsive Enzymes)
| 효소 | Fold Change | 해석 |
|------|-------------|------|
| **CHI (CHI1B2-2)** | **34× 증가** | Chalcone → Flavanone 전환 단계의 핵심 조절점 |
| **IFS (ifs1)** | **9× 증가** | Flavanone → Isoflavone 전환의 핵심 효소 |

### 2. 경로 조절 모델
```
에틸렌 신호
    │
    ▼
ERF/EIN 전사인자
    │
    ├──────────────────────┐
    ▼                      ▼
  CHI 증가 (34×)       IFS 증가 (9×)
    │                      │
    ▼                      ▼
Naringenin 축적    Isoflavone 생합성 증가
                           │
                           ▼
              말로닐화 이소플라본 4,000× 축적
```

### 3. 관련 전사인자 (Top 5)
| TF | Enzyme Connections | Description |
|----|-------------------|-------------|
| I1KWF7 | 3 | 다중 효소 조절 |
| A0A0R0GPT0 | 3 | 다중 효소 조절 |
| A0A0R0EQI8 | 1 | 단일 효소 조절 |
| C6TDI4 | 1 | 단일 효소 조절 |
| C6TEX3 | 1 | 단일 효소 조절 |

---

## 📁 생성된 파일 목록

### Daidzein 예측
- [isoflavonoid_pathway_enzymes.csv](./isoflavonoid_prediction/isoflavonoid_pathway_enzymes.csv)
- [pathway_associated_tfs.csv](./isoflavonoid_prediction/pathway_associated_tfs.csv)
- [pathway_proteomics_validation.csv](./isoflavonoid_prediction/pathway_proteomics_validation.csv)
- [pathway_summary.txt](./isoflavonoid_prediction/pathway_summary.txt)
- [pathway_diagram.png](./isoflavonoid_prediction/pathway_diagram.png)

### Genistein 예측
- [isoflavonoid_pathway_enzymes.csv](./genistein_prediction/isoflavonoid_pathway_enzymes.csv)
- [pathway_associated_tfs.csv](./genistein_prediction/pathway_associated_tfs.csv)
- [pathway_diagram.png](./genistein_prediction/pathway_diagram.png)

### Formononetin 예측
- [isoflavonoid_pathway_enzymes.csv](./formononetin_prediction/isoflavonoid_pathway_enzymes.csv)
- [pathway_associated_tfs.csv](./formononetin_prediction/pathway_associated_tfs.csv)
- [pathway_diagram.png](./formononetin_prediction/pathway_diagram.png)

### Glycitein 예측
- [isoflavonoid_pathway_enzymes.csv](./glycitein_prediction/isoflavonoid_pathway_enzymes.csv)
- [pathway_associated_tfs.csv](./glycitein_prediction/pathway_associated_tfs.csv)
- [pathway_diagram.png](./glycitein_prediction/pathway_diagram.png)

---

## 🔧 재현 명령어

```bash
# Daidzein
python src/generate_pathway_report.py --target Daidzein --pathway isoflavonoid --top-k 5 --output-dir results/isoflavonoid_prediction

# Genistein
python src/generate_pathway_report.py --target Genistein --pathway isoflavonoid --top-k 5 --output-dir results/genistein_prediction

# Formononetin
python src/generate_pathway_report.py --target Formononetin --pathway isoflavonoid --top-k 5 --output-dir results/formononetin_prediction

# Glycitein
python src/generate_pathway_report.py --target Glycitein --pathway isoflavonoid --top-k 5 --output-dir results/glycitein_prediction
```

---

## 📈 통계 요약

| 항목 | 값 |
|------|-----|
| 예측된 대사체 | 4개 (Daidzein, Genistein, Formononetin, Glycitein) |
| 경로 단계 | 7단계 |
| 유의미한 효소 (p<0.05) | 2개 (CHI, IFS) |
| 최대 Fold Change | 34× (CHI) |
| 연결된 TF | 15개 |

---

*보고서 생성: GNN Link Prediction + Proteomics Integration*
