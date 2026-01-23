# Manuscript v2 변경 요약

## 버전 정보
- **Version**: 2.0 (Two-Track Restructured)
- **Date**: 2026-01-21
- **Location**: `/data/ethylene/manuscript_v2/`

---

## 핵심 변경 사항

### 1. Two-Track 구조 도입

| Track | 목적 | 방법 | 검증 |
|-------|------|------|------|
| **Track A** | 생합성 경로 분석 | GNN (HGT) | Proteomics |
| **Track B** | 대사체-단백질 상호작용 | Molecular Docking | 실험 필요 (가설 단계) |

### 2. 주요 수정 내용

#### Abstract
- 두 트랙을 명확히 구분하여 기술
- 각 트랙의 목적과 결과를 별도로 제시

#### Methods
- **Section 2.5**: Track A - GNN-Based Biosynthetic Enzyme Prioritization
- **Section 2.6**: Track B - Metabolite-Protein Interaction Screening
- 각 트랙에 명확한 disclaimer 추가

#### Results
- **Section 3.4**: Track A Results (GNN + Proteomics 검증)
- **Section 3.5**: Track B Results (Docking, 탐색적 가설)
- 두 트랙이 서로 다른 질문에 답한다는 점 명시

#### Discussion
- **Section 4.1**: Track A 토론 (생합성 경로 조정)
- **Section 4.2**: Track B 토론 (신호전달 가설)
- **Section 4.3**: 방법론적 명확화 (왜 도킹이 GNN을 검증할 수 없는지)

### 3. 추가된 핵심 메시지

```
⚠️ "GNN이 예측한 효소-대사체 관계"와 "도킹이 테스트하는 물리적 결합"은
   서로 다른 생물학적 질문에 답한다.

✓ GNN 검증 방법: Proteomics, 효소 활성 측정, 유전학적 연구
✗ 도킹은 GNN 검증에 적합하지 않음
```

---

## 파일 구조

```
/data/ethylene/manuscript_v2/
└── MANUSCRIPT_v2.md     # 재구조화된 논문 (Two-Track)
```

---

## 향후 작업

- [x] Supplementary Materials 업데이트 ✓ `SUPPLEMENTARY_MATERIALS_v2.md`
- [x] Figures 재구성 (Track A/B 분리) ✓ `FIGURE_DESCRIPTIONS_v2.md`
- [x] 공저자 검토 요청 ✓ `COAUTHOR_REVIEW_REQUEST.md`

## 추가 작업 (선택)

- [ ] 실제 Figure 이미지 생성 (Cytoscape, PyMOL 등)
- [ ] 영문 교정 (Grammarly 또는 전문 서비스)
- [ ] 저널 투고 준비
