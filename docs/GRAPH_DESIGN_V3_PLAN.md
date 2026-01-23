# 그래프 설계 v3 수정 계획

## 핵심 변경 사항 요약

| 항목 | v2 (현재) | v3 (수정) |
|------|----------|----------|
| **Tier 반영** | edge_weight (미사용) | **edge_type 분리** |
| **Tier-P 기준** | EC 3자리 (24,416) | **reaction distance ≤2** (~2,500) |
| **Tier-P 역할** | supervision + message | **message only** |
| **TF 엣지 의미** | interacts (모호) | **associates** (명시) |
| **Filtered set** | 미정의 | **Tier-R only** |
| **Node features** | random embedding | **omics features** |

---

## 1. Edge Type 분리 (Q1 답변 반영)

### 변경 전
```python
data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = [1.0, 0.5, ...]
# HGTConv가 edge_weight 무시
```

### 변경 후
```python
# Tier-R: supervision용 (라벨로 사용)
data['Enzyme', 'catalyzes_R', 'Metabolite'].edge_index = tier_r_edges

# Tier-P: message passing only (라벨에서 제외)
data['Enzyme', 'catalyzes_P', 'Metabolite'].edge_index = tier_p_edges

# HGT가 자동으로 타입별 파라미터 학습
# W_catalyzes_R ≠ W_catalyzes_P
```

### 장점
1. **구현 단순화**: HGTConv 수정 불필요
2. **해석 용이**: 타입별 attention 분석 가능
3. **Semi-positive 자연 구현**: Tier-P를 train_pos_edges에서 제외하면 끝

---

## 2. Tier-P 기준 강화: Reaction Distance (Q3 답변 반영)

### 현재 문제
- EC 3자리 공유 = 24,416 엣지
- **Noise가 signal의 7.6배**
- 기질 특이성 무시로 false positive 다수

### 수정: Compound-only Reaction Graph

```python
# 1. Reaction graph 구축
adj = defaultdict(set)
for rxn, mets in reaction_to_mets.items():
    for sub in mets['substrates']:
        for prod in mets['products']:
            adj[sub].add(prod)
            adj[prod].add(sub)

# 2. Tier-P: distance ≤ 2
def is_tier_p(met1, met2, adj, max_dist=2):
    distances = bfs_distances(met1, adj, max_dist)
    return 0 < distances.get(met2, float('inf')) <= max_dist

# 3. 결과: ~2,500 pairs (10배 감소)
```

### 근거 명시 (문서용)
```
Tier-R: KEGG reaction에서 EC-metabolite가 직접 연결된 경우
        (반응 기반 직접 증거, supervision용)

Tier-P: KEGG reaction graph에서 거리 ≤ 2인 metabolite 쌍
        (반응 경로 인접성 기반 약한 prior, message passing용)
        Label noise 가능성이 있어 supervision에서 제외
```

---

## 3. Filtered Ranking (Q2 답변 반영)

### Protocol
```python
def evaluate_filtered(model, test_edges, known_positives):
    """
    known_positives: Tier-R edges only
    (KEGG reaction-grounded evidence)
    """
    for src, dst in test_edges:
        # 모든 candidate에 대해 score 계산
        scores = model.predict(src, all_metabolites)

        # Known positives 제거 (filtered)
        for known_dst in known_positives[src]:
            if known_dst != dst:  # test edge 자신은 유지
                scores[known_dst] = -inf

        # Ranking 계산
        rank = (scores > scores[dst]).sum() + 1
        ...
```

### 문서 명시
```
평가에서 Tier-R만 filtered set으로 사용.
Tier-P는 noise 가능성이 있어 filtered에서 제외.
이로 인해 일부 true positive가 후보로 남을 수 있으나,
이는 KEGG annotation incompleteness의 한계로 명시.
```

---

## 4. TF 엣지 의미 명확화

### 변경
```python
# 엣지 타입 이름
('TF', 'interacts', 'Enzyme')   # v2
('TF', 'associates', 'Enzyme')  # v3

# 또는 메타데이터
edge_metadata['evidence_type'] = 'functional_association'
```

### 문서 명시
```
STRING-DB 기반 TF-Enzyme 엣지는 functional association을 나타냄.
이는 "TF가 효소 유전자를 직접 조절한다"를 의미하지 않음.

TF 후보는 다음 검증 단계를 거쳐야 함:
1. Promoter motif 존재 확인 (JASPAR/PlantTFDB)
2. 발현 상관관계 (r > 0.5)
3. 문헌 증거

검증 전까지는 "가설 수준"으로만 해석.
```

---

## 5. 구현 순서

### Phase 1: 그래프 재구축 (Critical)

```bash
# 1. data_pipeline.py 수정
#    - Tier-R/P edge type 분리
#    - Tier-P 기준을 reaction distance로 변경
#    - TF 엣지 이름 변경

# 2. 그래프 재생성
python src/data_pipeline.py --output data/processed/graph_v3.pt
```

### Phase 2: 학습 코드 수정

```python
# train.py 수정
# - train_pos_edges = catalyzes_R만 사용
# - message passing은 catalyzes_R + catalyzes_P 모두
# - filtered evaluation 추가
```

### Phase 3: 노드 특성 추가

```python
# Metabolite features
met_features = torch.stack([
    log2fc,           # from MTBLS531
    -torch.log10(pvalue),
    pathway_multihot, # KEGG pathways
], dim=1)

# Enzyme features
enz_features = torch.stack([
    proteomics_log2fc,
    -torch.log10(proteomics_pvalue),
], dim=1)
```

---

## 6. 검증 계획

### A. Tier-P 변경 효과
```bash
# EC 3자리 vs Reaction distance 비교
python train.py --tier-p ec3     # 현재
python train.py --tier-p rxn2    # 수정
```

### B. Semi-positive 효과
```bash
# Tier-P를 라벨에 포함 vs 제외
python train.py --tier-p-supervision true   # 현재
python train.py --tier-p-supervision false  # 수정
```

### C. Filtered vs Raw ranking
```bash
# 평가 프로토콜 비교
python eval.py --ranking raw       # 현재
python eval.py --ranking filtered  # 수정
```

---

## 7. 예상 효과

| 지표 | v2 (현재) | v3 (예상) |
|------|----------|----------|
| Tier-P noise | 높음 | **~10배 감소** |
| 평가 신뢰도 | 낮음 (leakage 가능) | **높음 (filtered)** |
| 해석 가능성 | 낮음 | **타입별 분석 가능** |
| 재현성 | 중간 | **높음 (프로토콜 명시)** |

---

## 버전 정보

- 계획 버전: 3.0
- 작성일: 2026-01-23
- 상태: 승인 대기
