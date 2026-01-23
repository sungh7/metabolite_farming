# 공저자 검토 요청서

**To**: 공저자 전원  
**From**: [제1저자]  
**Date**: 2026-01-21  
**Subject**: 논문 v2.0 구조 변경 검토 요청

---

## 요청 사항

논문 원고에 중요한 구조적 변경을 가했습니다. 귀하의 검토와 의견을 요청드립니다.

---

## 변경 배경

### 발견된 문제

기존 논문(v1)에서 GNN 분석과 분자 도킹을 "통합 파이프라인"으로 제시했으나, 다음과 같은 방법론적 문제가 확인되었습니다:

1. **GNN 분석**은 대사 경로에서의 기능적 관계(효소-대사체)를 예측합니다.
2. **분자 도킹**은 최종 대사체(생성물)와 단백질 간의 물리적 결합을 테스트합니다.

> **핵심 문제**: 효소는 일반적으로 기질(substrate)에 결합하지, 생성물(product)에 결합하지 않습니다. 따라서 GNN이 예측한 효소에 생성물을 도킹하는 것은 "검증"이 아닙니다.

예시:
- GNN 예측: "IFS가 Daidzein과 관련있다" (= IFS가 Daidzein을 합성한다)
- 도킹 테스트: "Daidzein이 IFS에 결합하는가?" (= 생성물 억제?)
- **두 분석이 다른 질문에 답함**

---

## 제안된 해결책: Two-Track 구조

논문을 두 개의 독립적인 분석 트랙으로 재구성했습니다:

### Track A: 생합성 경로 분석 (Biosynthetic Pathway Analysis)

| 항목 | 내용 |
|------|------|
| **질문** | "에틸렌이 어떤 효소를 활성화해서 이소플라보노이드를 합성하는가?" |
| **방법** | GNN (Heterogeneous Graph Transformer) |
| **검증** | Proteomics 발현 데이터 |
| **결과** | PAL, 4CL, CHS, CHI, IFS, IFR 모두 상향조절 확인 (P<0.05) |

### Track B: 대사체-단백질 상호작용 스크리닝 (Metabolite-Protein Interaction)

| 항목 | 내용 |
|------|------|
| **질문** | "이소플라보노이드가 세포 내 다른 단백질에 결합해서 조절 효과를 일으키는가?" |
| **방법** | 분자 도킹 (AutoDock Vina) |
| **검증** | 실험적 검증 필요 (SPR, ITC, MST) |
| **결과** | Daidzein-FNR, Formononetin-Kinase 결합 예측 (가설 단계) |

**핵심**: 두 트랙은 서로 다른 질문에 답하며, 독립적으로 분석됩니다.

---

## 주요 변경 내용

### Abstract
- 두 트랙을 명확히 구분하여 기술
- 각 트랙의 목적과 결과를 별도로 제시

### Methods
- **Section 2.5**: Track A 방법론 (GNN)
- **Section 2.6**: Track B 방법론 (Docking) - GNN 검증이 아님을 명시

### Results
- **Section 3.4**: Track A 결과 + Proteomics 검증
- **Section 3.5**: Track B 결과 (탐색적 가설, 실험 검증 필요)

### Discussion
- **Section 4.3**: 방법론적 명확화 추가

### Supplementary Materials
- 새로운 Figure S2: Two-Track Framework
- 새로운 Figure S6: 왜 도킹이 GNN을 검증할 수 없는지
- 새로운 Table S4-S6: Track별 상세 데이터

---

## 검토 파일

| 파일 | 설명 |
|------|------|
| `manuscript_v2/MANUSCRIPT_v2.md` | 재구조화된 논문 본문 |
| `manuscript_v2/SUPPLEMENTARY_MATERIALS_v2.md` | 새로운 보충 자료 |
| `manuscript_v2/FIGURE_DESCRIPTIONS_v2.md` | 그림 설명 |
| `manuscript_v2/VERSION_CHANGES.md` | 변경 사항 요약 |

---

## 검토 요청 항목

다음 항목에 대해 의견을 부탁드립니다:

### 1. 구조 변경 승인
- [ ] Two-Track 구조에 동의하십니까?
- [ ] 트랙 명칭 (Track A/B)이 적절합니까?

### 2. Technical 검토
- [ ] GNN 분석 방법론 설명이 정확합니까?
- [ ] 도킹 분석의 한계 설명이 적절합니까?
- [ ] 통계 분석 및 검증 방법이 올바릅니까?

### 3. 해석 검토
- [ ] Track B 결과를 "가설"로 제시하는 것이 적절합니까?
- [ ] FNR-Daidzein 상호작용 해석이 생물학적으로 타당합니까?

### 4. 추가 제안
- [ ] 누락된 내용이 있습니까?
- [ ] 삭제해야 할 내용이 있습니까?
- [ ] 표현 수정이 필요한 부분이 있습니까?

---

## 응답 기한

**2026-01-28 (7일 후)**까지 의견을 회신해 주시면 감사하겠습니다.

---

## 연락처

질문이 있으시면 언제든 연락 주십시오.

- Email: [email@example.com]
- Phone: [phone number]

---

감사합니다.

[제1저자 서명]

---

*이 문서는 `/data/ethylene/manuscript_v2/` 디렉토리에서 확인하실 수 있습니다.*
