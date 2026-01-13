# 업무 지시서 (v1.0): Multi-Omics Ethylene Network

프로젝트: Untargeted 메타볼로믹스를 활용한 에틸렌 유도 콩 대사체 네트워크 규명과 기능성 대사산물 발굴
목표: 폭증 대사체 발굴, Multi-omics 기능적 연관 보강, HGT 링크 예측 신규성 검증.

---

## 0. 공통 원칙
1.  **메인 결론**: ET(ethylene/ACC) 축 실데이터에서만 도출. 외부 데이터는 "보강 근거(orthogonal support)"로만 사용.
2.  **표현 표준화**:
    *   금지: Validated the interaction / proves binding
    *   허용: supports functional association / co-regulation / prioritizes candidate regulator
3.  **재현성**: seed=42 고정, split 스키마 문서화, 로그 커밋.
4.  **Novel 정의**: (문헌/규제 DB 근거 없음) + (Text-mining 배제) + (외부 오믹스 지지 1개 이상).

---

## 1. 역할 및 산출물 (Roles & Deliverables)

### A. 메타볼로믹스 담당 (LC-MS)
*   **QC Rule**: 결측률 < 30% (Sample >= 70% detection), Log2 transform, Median normalization.
*   **Selection**: Log2FC >= 1.0 AND FDR <= 0.05.
*   **Outputs**:
    *   `data/metabolomics/processed/peak_table.tsv`
    *   `results/metabolomics/top_features_up.tsv`
    *   `results/metabolomics/annot_levels.tsv`

### B. 프로테오믹스 담당 (PXD)
*   **DiffExp**: |Log2FC| >= 0.58, Q <= 0.1 (or 0.05).
*   **Integration**: Sign-consistency Check (Protein ↑ & Metabolite ↑).
*   **Outputs**:
    *   `results/proteomics/protein_diffexp.tsv`
    *   `results/integration/sign_consistency_links.tsv`

### C. GNN/Network 담당 (HGT Baseline)
*   **Graphs**: Strict (No Text-Mining) vs Full.
*   **Baselines**: MLP, GCN/SAGE(Homo), SAGE(Hetero-Simple), Pathway Heuristic.
*   **Ablation**: Relation-Type Attention Removal etc.
*   **Case Study**: Top-1 Novel Pair using Strict Graph.
*   **Outputs**:
    *   `results/gnn/performance_table.tsv`
    *   `results/case_study/top1_novel_pair.md`

### D. 조절근거/Motif 담당 (Supplement)
*   **Protocol**: TSS upstream 2kb Scan for TF Motifs.
*   **Output**: `results/motif/motif_scan_summary.tsv`

### E. MD/Structure 담당 (Supplement)
*   **Protocol**: Docking -> MD (50ns).
*   **Metrics**: Ligand retention, Contact occupancy, RMSD.
*   **Status**: Optional (Post-MVP).

### F. Reporting/Editing (Defensive Writing)
*   **Key Phrases**:
    *   "Functional association" instead of "Interaction".
    *   "Proxy-labeled biochemical adjacency".
*   **Tables**: Baseline Comparison, Ablation Summary.

---

## 2. 통합 파이프라인 (Execution Order)
1.  **A/B**: Data Processing (Features & Proteins).
2.  **C**: Graph Construction (Strict/Full) -> GNN Training -> Baseline Comparison.
3.  **C/D**: Case Study Selection (Strict Top-1) -> Motif Support.
4.  **F**: Final Manuscript Assembly.

---

## 3. Directory Structure
*   `data/`: raw vs processed separation.
*   `results/`: TSV/CSV format.
*   `docs/`: Method descriptions.
*   `b_src/`: Source code (modularized by role A/B/C/D).

---

## 4. MVP Goals (Immediate Actions)
1.  [ ] **A**: Finalize "Top Up/Down Feature List" (Metabolomics).
2.  [ ] **C**: Build "Strict Graph" (No Text-Mining) & Run Baselines.
3.  [ ] **C**: Identify "Top-1 Novel Case Study".
4.  [ ] **F**: Update Methods with "Leakage/Split" defense text.
