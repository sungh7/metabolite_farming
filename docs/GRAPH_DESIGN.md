# Ethylene-Isoflavonoid GNN 그래프 설계 문서

## 목차
1. [개요](#1-개요)
2. [그래프 구조](#2-그래프-구조)
3. [데이터 소스](#3-데이터-소스)
4. [ID 매핑 파이프라인](#4-id-매핑-파이프라인)
5. [엣지 티어 시스템](#5-엣지-티어-시스템)
6. [EC 매핑 및 네거티브 샘플링](#6-ec-매핑-및-네거티브-샘플링)
7. [그래프 빌드 파이프라인](#7-그래프-빌드-파이프라인)
8. [설정 및 하이퍼파라미터](#8-설정-및-하이퍼파라미터)

---

## 1. 개요

### 1.1 프로젝트 목표
에틸렌 처리에 의한 대두(Glycine max)의 이소플라보노이드 생합성 경로를 예측하기 위한 이종 그래프 신경망(Heterogeneous GNN) 구축.

### 1.2 핵심 과제
- **Link Prediction**: 효소(Enzyme)와 대사물질(Metabolite) 간의 촉매 관계 예측
- **Multi-omics 통합**: STRING-DB(단백질 상호작용), KEGG(대사 경로), 실험 데이터(MTBLS531) 통합

### 1.3 그래프 개요

```
┌─────────────────────────────────────────────────────────────────┐
│                     Heterogeneous Graph                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   [Enzyme] ←──interacts──→ [Enzyme]                            │
│      │                        │                                 │
│      │ interacts              │ interacts                       │
│      ↓                        ↓                                 │
│    [TF]                   [Protein]                            │
│                                                                 │
│   [Enzyme] ───catalyzes───→ [Metabolite]                       │
│             (Tier-R: 1.0)                                       │
│             (Tier-P: 0.5)                                       │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 2. 그래프 구조

### 2.1 노드 타입

| 노드 타입 | 설명 | 데이터 소스 | 노드 수 |
|----------|------|------------|--------|
| **Enzyme** | 대사 반응을 촉매하는 효소 | STRING-DB, KEGG | ~3,425 |
| **Metabolite** | 대사물질 (기질/생성물) | KEGG, MTBLS531 | ~306 |
| **TF** | 전사인자 (Transcription Factor) | STRING-DB | ~100+ |
| **Protein** | 일반 단백질 | STRING-DB | Variable |

### 2.2 노드 특성 (Features)

```python
# 모든 노드는 64차원 학습 가능한 임베딩으로 초기화
data['Enzyme'].x = torch.randn(num_enzymes, 64)
data['Metabolite'].x = torch.randn(num_metabolites, 64)

# 메타데이터
data['Metabolite'].compound_ids = ['C02495', 'C00858', ...]  # KEGG compound IDs
```

### 2.3 엣지 타입

| 엣지 타입 | 방향 | 가중치 | 설명 |
|----------|------|--------|------|
| `interacts` | Enzyme ↔ Enzyme | 700-1000 | STRING-DB PPI (양방향) |
| `interacts` | Enzyme ↔ TF | 700-1000 | 효소-전사인자 상호작용 |
| `interacts` | Enzyme ↔ Protein | 700-1000 | 효소-단백질 상호작용 |
| **`catalyzes`** | Enzyme → Metabolite | 0.5-1.0 | **핵심: 촉매 관계** |
| `rev_catalyzes` | Metabolite → Enzyme | 0.5-1.0 | 역방향 촉매 엣지 |

### 2.4 PyTorch Geometric 데이터 구조

```python
from torch_geometric.data import HeteroData

data = HeteroData()

# 노드
data['Enzyme'].x = torch.randn(3425, 64)
data['Enzyme'].num_nodes = 3425
data['Metabolite'].x = torch.randn(306, 64)
data['Metabolite'].num_nodes = 306
data['Metabolite'].compound_ids = [...]  # KEGG IDs

# 엣지
data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = torch.tensor([
    [enzyme_indices],   # source
    [metabolite_indices]  # target
], dtype=torch.long)
data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = torch.tensor([...])

# 역방향 엣지 (메시지 패싱용)
data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = ...
```

---

## 3. 데이터 소스

### 3.1 STRING-DB v12.0 (단백질 상호작용 네트워크)

**소스 파일:**
```
data/raw/
├── 3847.protein.info.v12.0.txt.gz      # 단백질 메타데이터
├── 3847.protein.links.full.v12.0.txt.gz # PPI 엣지 (신뢰도 점수)
└── 3847.protein.aliases.v12.0.txt.gz    # ID 별칭 (NCBI ↔ STRING)
```

**PPI 신뢰도 점수 계산 (Strict Mode):**
```python
# Text-mining 제외, 재계산
evidence_channels = [
    'neighborhood',      # 유전체 인접성
    'fusion',           # 유전자 융합
    'cooccurence',      # 계통 발생적 공존
    'coexpression',     # 공동 발현
    'experimental',     # 실험적 증거
    'database',         # 데이터베이스 증거
]
# combined_score >= 700 필터링 (고신뢰도)
```

**노드 분류 기준 (graph_builder.py):**
```python
node_classification = {
    'Enzyme': [
        'ETR1', 'ETR2', 'ERS1', 'EIN2', 'CTR1',  # 에틸렌 신호전달
        'PAL', 'C4H', '4CL', 'CHS', 'CHI', 'IFS', 'HID',  # 이소플라보노이드
        '*synthase', '*kinase', '*transferase', '*reductase'  # 일반 효소
    ],
    'TF': ['WRKY*', 'MYB*', 'bHLH*', 'ERF*', '*zinc finger*'],
    'Protein': 'default'
}
```

### 3.2 KEGG (대사 경로 데이터베이스)

**소스 파일:**

| 파일 | 레코드 수 | 내용 |
|-----|----------|------|
| `gene_ec_mapping.tsv` | 7,938 | Gene ID ↔ EC 번호 |
| `full_enzyme_metabolite_edges.tsv` | 5,302 | EC ↔ Metabolite 반응 |
| `kegg_uniprot_mapping.csv` | 34,857 | KEGG Gene ↔ UniProt |
| `metabolites.csv` | 26 | 실험 대사물질 (MTBLS531) |

**gene_ec_mapping.tsv 스키마:**
```
gene_id      ec
100037445    1.7.3.3
100037447    1.11.1.6
100801944    2.4.1.13
```

**full_enzyme_metabolite_edges.tsv 스키마:**
```
enzyme_ec    metabolite_id    reaction_id    is_substrate    is_product
1.1.1.1      C00772           R02246         False           True
1.1.1.1      C00418           R02246         True            False
2.4.1.170    C02495           R08023         False           True
```

### 3.3 MTBLS531 (실험 데이터)

**metabolites.csv:**
```csv
compound_id,name,chebi,log2fc,pvalue,pathways
C02495,Daidzein,CHEBI:28197,0.1417,7.39e-07,map00943;map01110
C00858,Formononetin,CHEBI:18088,0.1348,3.80e-08,map00943;map01110
C02659,Genistein,CHEBI:28088,0.2105,1.23e-06,map00943;map01110
C12625,Glycitein,CHEBI:5373,0.1892,2.45e-05,map00943
C00062,L-Arginine,CHEBI:16467,0.3289,7.29e-05,map00220
```

**핵심 이소플라보노이드:**
- **Daidzein (C02495)**: 대두 이소플라본의 주요 성분
- **Formononetin (C00858)**: 메틸화된 이소플라본
- **Genistein (C02659)**: 항산화 활성 이소플라본
- **Glycitein (C12625)**: 대두 특이적 이소플라본

---

## 4. ID 매핑 파이프라인

### 4.1 문제점
서로 다른 데이터베이스는 다른 ID 체계를 사용:
- **STRING-DB**: `3847.A0A075W8S1` (종.UniProt)
- **KEGG**: NCBI Gene ID (예: `100037445`)
- **UniProt**: `A0A075W8S1`

### 4.2 3단계 매핑 전략

```
┌─────────────────────────────────────────────────────────────────┐
│                    ID Mapping Pipeline                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  KEGG Gene ID                                                   │
│       │                                                         │
│       ├──[1차]──→ kegg_uniprot_mapping.csv ──→ UniProt ID      │
│       │                                            │            │
│       │                                            ↓            │
│       │                              enzyme_string_mapping.csv  │
│       │                                            │            │
│       │                                            ↓            │
│       │                                      Enzyme Index       │
│       │                                                         │
│       └──[2차]──→ protein.aliases.gz ──→ STRING ID             │
│                   (NCBI GeneID source)        │                 │
│                                               ↓                 │
│                                    string_to_enzyme_idx         │
│                                               │                 │
│                                               ↓                 │
│                                         Enzyme Index            │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 4.3 매핑 코드 (data_pipeline.py)

```python
# 1차: UniProt 매핑 (선호)
uniprot = kegg_to_uniprot.get(gene_id)
if uniprot and uniprot in uniprot_to_enzyme_idx:
    enzyme_idx = uniprot_to_enzyme_idx[uniprot]
    matched_uniprot += 1

# 2차: NCBI Direct 매핑 (폴백)
elif gene_id in ncbi_to_enzyme_idx:
    enzyme_idx = ncbi_to_enzyme_idx[gene_id]
    matched_ncbi += 1

else:
    unmatched += 1
```

### 4.4 매핑 통계 (일반적)

```
Matched via UniProt: 1,139  (85%)
Matched via NCBI:      131  (10%)
Unmatched:           6,667  (5%)
────────────────────────────
ECs with enzyme mapping: 290
Unique enzymes with EC: 1,147
```

---

## 5. 엣지 티어 시스템

### 5.1 2티어 증거 기반 아키텍처

효소-대사물질 엣지는 **생화학적 증거 수준**에 따라 분류:

| 티어 | 이름 | 증거 | 가중치 | 설명 |
|-----|------|------|--------|------|
| **R** | Reaction-grounded | KEGG 직접 반응 | **1.0** | EC가 해당 대사물질과 직접 반응 |
| **P** | Pathway-supported | EC class 유사성 | **0.5** | 동일 EC class (3자리)의 다른 반응 |

### 5.2 EC 번호 체계

```
EC X.X.X.X
   │ │ │ └── 특정 기질 (Serial number)
   │ │ └──── 아그룹 (Sub-subclass)
   │ └────── 서브클래스 (Subclass)
   └──────── 클래스 (Main class)

예: EC 1.1.1.1 (Alcohol dehydrogenase)
    │ │ │ └── 특정: ethanol → acetaldehyde
    │ │ └──── CH-OH 그룹에 작용
    │ └────── donor가 CH-OH
    └──────── 산화환원효소 (Oxidoreductases)
```

### 5.3 티어 구축 알고리즘

```python
# Tier-R: 직접 반응 링크
for ec, metabolite in kegg_reactions:
    if ec in ec_to_enzyme_indices:
        for enzyme_idx in ec_to_enzyme_indices[ec]:
            tier_r_edges.add((enzyme_idx, metabolite_idx))
            # weight = 1.0

# Tier-P: EC class 유사성 기반
for ec in ec_to_enzyme_indices:
    ec_class = '.'.join(ec.split('.')[:3])  # "1.1.1"

    for other_ec, other_metabolites in ec_to_metabolites.items():
        if other_ec.startswith(ec_class) and other_ec != ec:
            for met_idx in other_metabolites:
                for enz_idx in ec_to_enzyme_indices[ec]:
                    tier_p_edges.add((enz_idx, met_idx))
                    # weight = 0.5

# 중복 제거 (Tier-R 우선)
tier_p_edges = tier_p_edges - tier_r_edges
```

### 5.4 예시: Daidzein (C02495) 엣지

```
대사물질: Daidzein (C02495)
관련 EC: 2.4.1.170 (Isoflavone 4'-O-methyltransferase)

Tier-R 엣지 (직접 반응):
  2.4.1.170 → C02495  (weight=1.0)
  2.4.1.170 → C00858  (weight=1.0)  # Formononetin

Tier-P 엣지 (EC class 2.4.1.* 유사):
  2.4.1.123 → C02495  (weight=0.5)
  2.4.1.234 → C02495  (weight=0.5)
```

### 5.5 티어 통계

```
Tier-R edges: 3,195  (직접 반응, high confidence)
Tier-P edges: 24,416 (EC class 유사, medium confidence)
─────────────────────
Total edges: 27,611
Mean edge weight: 0.558
```

---

## 6. EC 매핑 및 네거티브 샘플링

### 6.1 문제점: False Negative

**기존 Hard Negative Sampling의 문제:**
```python
# 기존 방식: 인접 인덱스 샘플링
offset = torch.randint(1, 6, (num_neg,))  # ±1~5
neg_dst = (pos_dst + offset) % num_metabolites
```

**문제점:**
- 대사물질 인덱스는 **KEGG compound ID 알파벳순** 정렬
- 인접 인덱스가 **동일 pathway**일 확률 높음
- 측정된 False Negative Rate: **~35%**

### 6.2 EC-class 기반 샘플링 해결책

**원리:** 동일 EC class의 대사물질은 유사한 효소가 처리 → 이들을 **피해서** 샘플링

**그래프에 저장되는 매핑:**
```python
# EC → 대사물질 인덱스
data['Metabolite'].ec_to_indices = {
    '1.1.1.1': [44, 55, 89],      # EC 1.1.1.1에 관련된 대사물질
    '1.1.1.2': [23, 67],
    '2.4.1.170': [0, 5, 12],      # Daidzein 관련
    ...
}

# 대사물질 → EC classes (역매핑)
data['Metabolite'].met_to_ecs = {
    0: {'2.4.1.170', '1.14.14.87'},  # Daidzein의 EC classes
    1: {'2.4.1.170'},                 # Formononetin
    ...
}
```

### 6.3 EC-class 샘플링 알고리즘

```python
def _sample_ec_class(self, pos_edge_index, num_target_nodes, num_neg, device):
    neg_dst_list = []
    all_mets = set(range(num_target_nodes))

    for pos_dst_idx in pos_edge_index[1]:
        # 현재 대사물질의 EC classes
        current_ecs = self.met_to_ecs.get(pos_dst_idx, set())

        # 동일 EC class의 대사물질 제외
        excluded = set()
        for ec in current_ecs:
            excluded.update(self.ec_to_indices.get(ec, []))
        excluded.add(pos_dst_idx)  # 자기 자신 제외

        # 후보에서 랜덤 샘플링
        candidates = list(all_mets - excluded)
        if candidates:
            neg_dst_list.append(random.choice(candidates))
        else:
            # 폴백: 랜덤
            neg_dst_list.append(random.randint(0, num_target_nodes - 1))

    return neg_src, torch.tensor(neg_dst_list, device=device)
```

### 6.4 성능 비교

| 지표 | Hard Sampling | EC-class Sampling |
|-----|---------------|-------------------|
| False Negative Rate | ~35% | **< 1%** |
| 구현 복잡도 | 낮음 | 중간 |
| 추가 메모리 | 없음 | EC 매핑 저장 |

### 6.5 설정 (config.py)

```python
NEGATIVE_SAMPLING_CONFIG = {
    'strategy': 'ec_class',        # 권장: EC-class 기반
    'ec_class_fallback': 'random', # EC 정보 없을 때 폴백
    'hard_offset_min': 1,          # (레거시) hard sampling용
    'hard_offset_max': 5,
    'hard_ratio': 0.5,             # mixed 전략용
    'neg_ratio': 1.0,              # negative:positive 비율
}
```

---

## 7. 그래프 빌드 파이프라인

### 7.1 메인 파이프라인 (data_pipeline.py)

**실행:**
```bash
python src/data_pipeline.py --output data/processed/enhanced_bipartite_graph_v2.pt
```

**단계별 처리:**

```
Step 1: Base PPI 그래프 로드
        └── data/processed/strict_graph.pt

Step 2: STRING-NCBI 매핑 로드
        └── 3847.protein.aliases.v12.0.txt.gz
        └── 71,808 NCBI→STRING 매핑

Step 3: EC → Enzyme 매핑 생성
        ├── gene_ec_mapping.tsv (7,938 entries)
        ├── kegg_uniprot_mapping.csv
        └── 290 EC classes → 1,147 unique enzymes

Step 4: KEGG 엣지 로드
        └── full_enzyme_metabolite_edges.tsv (5,302 edges)

Step 5: 대사물질 선택
        ├── 실험 대사물질 (MTBLS531): 10개 (KEGG 엣지 있는 것)
        └── Top-degree 컨텍스트: 296개
        └── 총 306개 대사물질

Step 6: Tier-R 엣지 구축
        └── 3,195 직접 반응 엣지 (weight=1.0)

Step 7: Tier-P 엣지 구축
        └── 24,416 EC class 유사 엣지 (weight=0.5)

Step 8: EC 매핑 저장
        ├── ec_to_indices: 914 EC classes
        └── met_to_ecs: 306 metabolites

Step 9: HeteroData 저장
        └── enhanced_bipartite_graph_v2.pt
```

### 7.2 출력 그래프 구조

```python
HeteroData(
    Enzyme={
        x=[3425, 64],           # 노드 특성
        num_nodes=3425
    },
    Metabolite={
        x=[306, 64],
        num_nodes=306,
        compound_ids=[...],     # KEGG IDs
        ec_to_indices={...},    # EC → metabolite indices
        met_to_ecs={...}        # metabolite → EC classes
    },
    (Enzyme, interacts, Enzyme)={
        edge_index=[2, num_ppi]
    },
    (Enzyme, catalyzes, Metabolite)={
        edge_index=[2, 27611],
        edge_weight=[27611]     # 0.5 또는 1.0
    },
    (Metabolite, rev_catalyzes, Enzyme)={
        edge_index=[2, 27611],
        edge_weight=[27611]
    }
)
```

### 7.3 레거시 빌더 (참고용)

| 파일 | 상태 | 설명 |
|-----|------|------|
| `bipartite_builder.py` | Deprecated | 시뮬레이션 엣지 사용 |
| `tiered_bipartite_builder.py` | Deprecated | 시뮬레이션 매핑 |
| `enhanced_bipartite_builder.py` | Deprecated | data_pipeline.py의 전신 |

---

## 8. 설정 및 하이퍼파라미터

### 8.1 그래프 설정 (config.py)

```python
GRAPH_CONFIG = {
    # STRING-DB
    'ppi_threshold': 700,           # PPI 신뢰도 임계값 (0-1000)
    'strict_mode': True,            # Text-mining 증거 제외

    # 노드 특성
    'feature_dim': 64,              # 임베딩 차원

    # 엣지 가중치
    'tier_r_weight': 1.0,           # Reaction-grounded
    'tier_p_weight': 0.5,           # Pathway-supported

    # 대사물질 선택
    'max_context_metabolites': 300,

    # 제외할 Currency 대사물질
    'currency_metabolites': {
        'C00001',  # H2O
        'C00002',  # ATP
        'C00003',  # NAD+
        'C00004',  # NADH
        'C00005',  # NADPH
        'C00006',  # NADP+
        'C00008',  # ADP
        'C00009',  # Orthophosphate
        'C00010',  # CoA
        'C00011',  # CO2
        'C00013',  # Diphosphate
        'C00014',  # Ammonia
        'C00020',  # AMP
        'C00027',  # H2O2
        'C00044',  # GTP
        'C00080',  # H+
    }
}
```

### 8.2 학습 설정

```python
TRAINING_CONFIG = {
    'hidden_channels': 64,
    'out_channels': 64,
    'num_heads': 4,
    'num_layers': 2,
    'learning_rate': 0.01,
    'weight_decay': 1e-5,
    'epochs': 50,
    'patience': 10,
    'max_grad_norm': 1.0,
}
```

### 8.3 데이터 분할 설정

```python
SPLIT_CONFIG = {
    'strategy': 'node_split',
    'train_ratio': 0.70,
    'val_ratio': 0.15,
    'test_ratio': 0.15,
    'disjoint': True,  # Inductive 평가용
}
```

---

## 부록 A: 파일 구조

```
data/
├── raw/                              # STRING-DB 원본
│   ├── 3847.protein.info.v12.0.txt.gz
│   ├── 3847.protein.links.full.v12.0.txt.gz
│   └── 3847.protein.aliases.v12.0.txt.gz
│
├── kegg/                             # KEGG 데이터
│   ├── gene_ec_mapping.tsv           # 7,938 records
│   ├── full_enzyme_metabolite_edges.tsv  # 5,302 records
│   ├── kegg_uniprot_mapping.csv      # 34,857 records
│   └── metabolites.csv               # 26 experimental
│
└── processed/                        # 처리된 그래프
    ├── strict_graph.pt               # Base PPI
    ├── enhanced_bipartite_graph_v2.pt # 최종 그래프 (권장)
    ├── enzyme_string_mapping.csv     # 3,425 enzymes
    └── tf_string_mapping.csv         # TF 매핑

src/
├── config.py                 # 중앙 설정
├── data_pipeline.py          # 메인 그래프 빌더 (권장)
├── graph_builder.py          # STRING-DB PPI 처리
├── utils/
│   ├── negative_sampling.py  # 네거티브 샘플러
│   ├── data_split.py         # 데이터 분할
│   └── seed.py               # 시드 관리
└── train.py                  # 학습 파이프라인
```

---

## 부록 B: 그래프 통계 요약

| 항목 | 값 |
|-----|-----|
| **노드** | |
| Enzymes | 3,425 |
| Metabolites | 306 |
| - Experimental | 10 |
| - Context | 296 |
| **엣지** | |
| Tier-R (direct) | 3,195 |
| Tier-P (pathway) | 24,416 |
| Total | 27,611 |
| Mean weight | 0.558 |
| **EC 매핑** | |
| EC classes | 914 |
| Metabolites with EC | 306 (100%) |
| Avg ECs per metabolite | 7.9 |

---

## 부록 C: 사용법 예시

### 그래프 빌드
```bash
# EC 매핑 포함 그래프 생성
PYTHONPATH=/data/ethylene python src/data_pipeline.py \
    --output data/processed/enhanced_bipartite_graph_v2.pt
```

### 학습 실행
```bash
# EC-class 네거티브 샘플링으로 학습
PYTHONPATH=/data/ethylene python src/train.py \
    --graph data/processed/enhanced_bipartite_graph_v2.pt \
    --neg-strategy ec_class \
    --seeds 42,123,456
```

### 그래프 검증
```python
import torch

data = torch.load('data/processed/enhanced_bipartite_graph_v2.pt')

print(f"Enzymes: {data['Enzyme'].num_nodes}")
print(f"Metabolites: {data['Metabolite'].num_nodes}")
print(f"Edges: {data['Enzyme', 'catalyzes', 'Metabolite'].edge_index.shape[1]}")
print(f"EC classes: {len(data['Metabolite'].ec_to_indices)}")
```

---

*문서 버전: 2.0*
*최종 수정: 2025-01*
