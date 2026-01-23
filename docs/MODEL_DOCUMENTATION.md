# 🧠 GNN 모델 아키텍처 상세 문서

**프로젝트**: 에틸렌 유도 이소플라보노이드 생합성 예측 모델  
**작성일**: 2026년 1월 22일  
**목적**: 모델 아키텍처 및 코드 상세 설명

---

## 📋 목차

1. [모델 개요](#1-모델-개요)
2. [HGT (Heterogeneous Graph Transformer)](#2-hgt-heterogeneous-graph-transformer)
3. [Link Predictor](#3-link-predictor)
4. [보조 모델들](#4-보조-모델들)
5. [그래프 구축 파이프라인](#5-그래프-구축-파이프라인)
6. [학습 파이프라인](#6-학습-파이프라인)
7. [추론 파이프라인](#7-추론-파이프라인)
8. [전체 파이프라인 흐름](#8-전체-파이프라인-흐름)

---

## 1. 모델 개요

### 1.1 연구 목표

본 모델은 **대사체-효소 상호작용 예측**을 위한 Graph Neural Network 기반 Link Prediction 시스템입니다.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         MODEL PIPELINE OVERVIEW                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Input Data                                                                │
│   ┌────────────┐   ┌───────────┐   ┌───────────┐                           │
│   │ STRING DB  │   │   KEGG    │   │ Proteomics│                           │
│   │    PPI     │   │ Pathways  │   │ PXD006989 │                           │
│   └──────┬─────┘   └─────┬─────┘   └─────┬─────┘                           │
│          │               │               │                                  │
│          └───────────────┼───────────────┘                                  │
│                          ▼                                                  │
│   ┌──────────────────────────────────────────────────────────────────────┐  │
│   │              Heterogeneous Knowledge Graph                           │  │
│   │   Nodes: Enzyme (500+), Metabolite (200+), TF (100+), Protein        │  │
│   │   Edges: catalyzes, interacts, regulates (2-Tier: R/P)               │  │
│   └─────────────────────────┬────────────────────────────────────────────┘  │
│                             ▼                                               │
│   ┌──────────────────────────────────────────────────────────────────────┐  │
│   │                    HGT Encoder (2 Layers)                            │  │
│   │   • Type-specific Linear Projection                                  │  │
│   │   • Multi-head Attention (4 heads)                                   │  │
│   │   • Output: 64-dim node embeddings                                   │  │
│   └─────────────────────────┬────────────────────────────────────────────┘  │
│                             ▼                                               │
│   ┌──────────────────────────────────────────────────────────────────────┐  │
│   │                    Link Predictor                                    │  │
│   │   • Dot product scoring: score = embed(u) · embed(v)                 │  │
│   │   • Sigmoid → probability                                            │  │
│   └─────────────────────────┬────────────────────────────────────────────┘  │
│                             ▼                                               │
│   Output: Top-K predicted Enzyme-Metabolite interactions                    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 1.2 구현된 모델 목록

| 모델 | 클래스명 | 설명 | 역할 |
|------|----------|------|------|
| **HGT** ⭐ | `HGT` | Heterogeneous Graph Transformer | 메인 인코더 |
| HAN | `HAN` | Heterogeneous Attention Network | 비교 실험 |
| HeteroSAGE | `HeteroSAGE` | Heterogeneous GraphSAGE | Baseline |
| SimpleMLP | `SimpleMLP` | Multi-Layer Perceptron | Ablation |
| LinkPredictor | `LinkPredictor` | Dot product 기반 링크 예측 | 디코더 |

---

## 2. HGT (Heterogeneous Graph Transformer)

### 2.1 개요

HGT는 **이종 그래프**에서 노드/엣지 타입에 따라 다른 attention 파라미터를 학습하는 Graph Transformer입니다.

> **논문**: "Heterogeneous Graph Transformer" (WWW 2020)

### 2.2 아키텍처 다이어그램

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         HGT ARCHITECTURE                                     │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│ INPUT: x_dict                                                               │
│ ┌───────────────────────────────────────────────────────────────────────┐   │
│ │  {                                                                    │   │
│ │      'Enzyme':     [N_enz × 64]   ← 랜덤 초기화 임베딩                 │   │
│ │      'Metabolite': [N_met × 64]   ← 랜덤 초기화 임베딩                 │   │
│ │      'TF':         [N_tf × 64]    ← 랜덤 초기화 임베딩                 │   │
│ │      'Protein':    [N_prot × 64]  ← 랜덤 초기화 임베딩                 │   │
│ │  }                                                                    │   │
│ └───────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│ ┌───────────────────────────────────────────────────────────────────────┐   │
│ │  STAGE 1: Type-Specific Linear Projection (lin_dict)                  │   │
│ │  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐        │   │
│ │  │ W_enzyme: 64→64 │  │ W_met: 64→64    │  │ W_tf: 64→64     │        │   │
│ │  └─────────────────┘  └─────────────────┘  └─────────────────┘        │   │
│ │  각 노드 타입별로 독립적인 Linear layer 적용                           │   │
│ │  Output: x = ReLU(W_type · x)                                         │   │
│ └───────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│ ┌───────────────────────────────────────────────────────────────────────┐   │
│ │  STAGE 2: HGTConv Layer 1                                             │   │
│ │  ┌─────────────────────────────────────────────────────────────────┐  │   │
│ │  │  Multi-Head Type-Aware Attention (4 heads)                      │  │   │
│ │  │                                                                 │  │   │
│ │  │  For each edge_type (src_type, relation, dst_type):             │  │   │
│ │  │    • Q = W_Q[edge_type] × h[dst]  (Query from destination)      │  │   │
│ │  │    • K = W_K[edge_type] × h[src]  (Key from source)             │  │   │
│ │  │    • V = W_V[edge_type] × h[src]  (Value from source)           │  │   │
│ │  │                                                                 │  │   │
│ │  │    Attention = softmax(Q · K^T / √d)                            │  │   │
│ │  │    Message = Attention × V                                      │  │   │
│ │  │    Aggregate = Σ Messages (per destination node)                │  │   │
│ │  └─────────────────────────────────────────────────────────────────┘  │   │
│ └───────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│ ┌───────────────────────────────────────────────────────────────────────┐   │
│ │  STAGE 3: HGTConv Layer 2                                             │   │
│ │  [동일한 구조 - 2-hop 이웃 정보 집약]                                   │   │
│ └───────────────────────────────────────────────────────────────────────┘   │
│                              │                                              │
│                              ▼                                              │
│ OUTPUT: x_dict (Updated Embeddings)                                         │
│ ┌───────────────────────────────────────────────────────────────────────┐   │
│ │  {                                                                    │   │
│ │      'Enzyme':     [N_enz × 64]   ← 구조 정보 반영된 임베딩            │   │
│ │      'Metabolite': [N_met × 64]   ← 구조 정보 반영된 임베딩            │   │
│ │      'TF':         [N_tf × 64]    ← 구조 정보 반영된 임베딩            │   │
│ │  }                                                                    │   │
│ └───────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 2.3 코드 구현

**파일 위치**: [model.py](file:///data/ethylene/src/model.py)

```python
import torch
import torch.nn as nn
from torch_geometric.nn import HGTConv, Linear

class HGT(nn.Module):
    """
    Heterogeneous Graph Transformer
    
    핵심 특징:
    - 노드/엣지 타입별 서로 다른 attention 파라미터
    - Multi-head attention mechanism
    - Type-aware message passing
    
    Args:
        metadata: 그래프 메타데이터 (node_types, edge_types)
        in_channels: 입력 특징 차원 (64)
        hidden_channels: 은닉층 차원 (64)
        out_channels: 출력 차원 (64)
        num_heads: Attention head 수 (4)
        num_layers: Convolution 레이어 수 (2)
    """
    
    def __init__(self, metadata, in_channels, hidden_channels, 
                 out_channels, num_heads=4, num_layers=2):
        super().__init__()
        
        # 1. 노드 타입별 입력 투영 레이어
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:  # metadata[0] = node_types
            self.lin_dict[node_type] = Linear(in_channels, hidden_channels)
        
        # 2. HGT Convolution 레이어 스택
        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HGTConv(
                hidden_channels, 
                hidden_channels, 
                metadata, 
                heads=num_heads
            )
            self.convs.append(conv)
        
        # 3. 출력 레이어 (선택적)
        self.out_lin = Linear(hidden_channels, out_channels)
    
    def forward(self, x_dict, edge_index_dict):
        """
        Forward pass
        
        Args:
            x_dict: {node_type: [N × in_channels]} 노드 특징 딕셔너리
            edge_index_dict: {edge_type: [2, E]} 엣지 인덱스 딕셔너리
            
        Returns:
            x_dict: {node_type: [N × out_channels]} 업데이트된 노드 임베딩
        """
        # Step 1: 타입별 특징 투영 + ReLU 활성화
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
        }
        
        # Step 2: HGT Convolution 통과
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
        
        return x_dict
```

### 2.4 하이퍼파라미터

| 파라미터 | 값 | 설명 |
|----------|-----|------|
| `in_channels` | 64 | 입력 노드 특징 차원 |
| `hidden_channels` | 64 | 은닉층 차원 |
| `out_channels` | 64 | 출력 임베딩 차원 |
| `num_heads` | 4 | Multi-head attention 헤드 수 |
| `num_layers` | 2 | Convolution 레이어 수 (2-hop 이웃) |

---

## 3. Link Predictor

### 3.1 개요

Link Predictor는 HGT가 생성한 노드 임베딩으로부터 **두 노드 사이의 엣지 존재 확률**을 예측합니다.

### 3.2 아키텍처

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                      LINK PREDICTOR                                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Input: 노드 쌍 (source, destination)                                      │
│                                                                             │
│   ┌──────────────────┐  ┌──────────────────┐                               │
│   │  h_source [64]   │  │  h_dest [64]     │                               │
│   │  (예: Enzyme)    │  │  (예: Metabolite)│                               │
│   └────────┬─────────┘  └────────┬─────────┘                               │
│            │                     │                                          │
│            └──────────┬──────────┘                                          │
│                       │                                                     │
│                       ▼                                                     │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  Element-wise Product: h_source ⊙ h_dest                             │ │
│   │  [64] × [64] → [64]                                                   │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                       │                                                     │
│                       ▼                                                     │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  Sum Reduction: Σ(h_source ⊙ h_dest) → scalar                        │ │
│   │  This equals: dot_product(h_source, h_dest)                           │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                       │                                                     │
│                       ▼                                                     │
│   Output: Link Score (before sigmoid) → Probability (after sigmoid)        │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 3.3 코드 구현

```python
class LinkPredictor(nn.Module):
    """
    Dot Product 기반 Link Predictor
    
    Method: score(u, v) = embedding(u) · embedding(v)
    
    장점:
    - 파라미터 없음 (효율적)
    - 임베딩 유사도를 직접 활용
    - 대규모 그래프에 확장 가능
    """
    
    def __init__(self, in_channels):
        super().__init__()
        # 파라미터 없음 - dot product만 사용
        
    def forward(self, x_src, x_dst, edge_label_index):
        """
        Args:
            x_src: Source 노드 임베딩 [N_src × 64]
            x_dst: Destination 노드 임베딩 [N_dst × 64]
            edge_label_index: 예측할 엣지 [2, E]
            
        Returns:
            scores: 각 엣지의 점수 [E]
        """
        row, col = edge_label_index
        src_feats = x_src[row]    # [E × 64]
        dst_feats = x_dst[col]    # [E × 64]
        
        # Element-wise product + sum = dot product
        return (src_feats * dst_feats).sum(dim=-1)  # [E]
```

---

## 4. 보조 모델들

### 4.1 HAN (Heterogeneous Attention Network)

```python
class HAN(nn.Module):
    """
    Hierarchical Attention:
    1. Node-level attention
    2. Semantic-level attention (edge type)
    """
    def __init__(self, metadata, in_channels, hidden_channels, 
                 out_channels, num_heads=2, num_layers=2):
        super().__init__()
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = Linear(in_channels, hidden_channels)

        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HANConv(hidden_channels, hidden_channels, 
                          metadata, heads=num_heads)
            self.convs.append(conv)
        self.out_lin = Linear(hidden_channels, out_channels)

    def forward(self, x_dict, edge_index_dict):
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
        }
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
        return x_dict
```

### 4.2 HeteroSAGE

```python
class HeteroSAGE(nn.Module):
    """
    GraphSAGE의 이종 그래프 버전
    각 엣지 타입별로 SAGEConv 적용
    """
    def __init__(self, metadata, in_channels, hidden_channels, 
                 out_channels, num_layers=2):
        super().__init__()
        # 1. Projections
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = Linear(in_channels, hidden_channels)

        # 2. HeteroConvs
        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv_dict = {}
            for edge_type in metadata[1]:
                conv_dict[edge_type] = SAGEConv(hidden_channels, hidden_channels)
            self.convs.append(HeteroConv(conv_dict, aggr='sum'))
        self.out_lin = Linear(hidden_channels, out_channels)

    def forward(self, x_dict, edge_index_dict):
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
        }
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
            x_dict = {k: v.relu_() for k, v in x_dict.items()}
        return x_dict
```

### 4.3 SimpleMLP (Ablation 용)

```python
class SimpleMLP(nn.Module):
    """
    그래프 구조를 무시하는 MLP 베이스라인
    Ablation study에서 "그래프가 필요한가?" 검증용
    """
    def __init__(self, metadata, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = nn.Sequential(
                Linear(in_channels, hidden_channels),
                nn.ReLU(),
                Linear(hidden_channels, hidden_channels),
                nn.ReLU(),
                Linear(hidden_channels, out_channels)
            )

    def forward(self, x_dict, edge_index_dict):
        # edge_index_dict 무시 - 그래프 구조 사용 안함
        out_dict = {}
        for node_type, x in x_dict.items():
            if node_type in self.lin_dict:
                out_dict[node_type] = self.lin_dict[node_type](x)
        return out_dict
```

---

## 5. 그래프 구축 파이프라인

### 5.1 데이터 로더

**파일 위치**: [dataloader.py](file:///data/ethylene/src/dataloader.py)

```python
class StringDBLoader:
    """
    STRING Database에서 단백질-단백질 상호작용(PPI) 로드
    
    Data Sources:
    - 3847.protein.info.v12.0.txt.gz (단백질 정보)
    - 3847.protein.links.full.v12.0.txt.gz (상호작용)
    
    Species: Glycine max (Soybean) - ID: 3847
    """
    
    def __init__(self, raw_dir='data/raw'):
        self.raw_dir = raw_dir
        self.protein_info_path = os.path.join(raw_dir, 
            '3847.protein.info.v12.0.txt.gz')
        self.links_path = os.path.join(raw_dir, 
            '3847.protein.links.full.v12.0.txt.gz')
        self.protein_map = {}   # string_id -> preferred_name
        self.node_to_idx = {}   # string_id -> int index
        self.idx_to_node = {}   # int index -> string_id

    def load_protein_info(self):
        """단백질 정보 파싱: string_id -> name 매핑"""
        pass
    
    def load_interactions(self, threshold=700, strict=False):
        """
        상호작용 로드 및 필터링
        
        Args:
            threshold: 최소 combined_score (0-1000). 700 = high confidence
            strict: True면 text-mining 채널 제외
            
        Returns:
            edges: [(src_idx, dst_idx), ...]
        """
        pass
```

### 5.2 그래프 빌더

**파일 위치**: [graph_builder.py](file:///data/ethylene/src/graph_builder.py)

```python
def identify_node_type(annotation, preferred_name):
    """
    휴리스틱 기반 노드 타입 분류
    
    Priority:
    1. 에틸렌 신호전달 경로 (Signaling)
    2. 이소플라보노이드 효소 (Enzyme)
    3. 전사인자 (TF)
    4. 일반 효소 (Enzyme)
    5. 기타 (Protein)
    """
    ann = str(annotation).lower()
    name = str(preferred_name).lower()
    
    # Priority 1: Ethylene Signaling
    if any(k in name for k in ['etr1', 'ein2', 'ein3', 'ctr1', 'ebf1']):
        return 'Signaling'
        
    # Priority 2: Key Isoflavonoid Enzymes
    if any(k in name for k in ['pal', 'c4h', '4cl', 'chs', 'chi', 'ifs']):
        return 'Enzyme'

    # Priority 3: Transcription Factors
    if 'transcription factor' in ann or 'myb' in ann or 'wrky' in ann:
        return 'TF'
        
    # Priority 4: General Enzymes
    if 'synthase' in ann or 'kinase' in ann or 'transferase' in ann:
        return 'Enzyme'
        
    return 'Protein'


def build_graph(threshold=700, strict=False):
    """
    이종 그래프 구축
    
    Output Schema:
    ┌────────────────────────────────────────────────────────────┐
    │  HeteroData                                                │
    │  ├── Enzyme.x: [N_enz × 64]                               │
    │  ├── Metabolite.x: [N_met × 64]                           │
    │  ├── TF.x: [N_tf × 64]                                    │
    │  ├── Protein.x: [N_prot × 64]                             │
    │  ├── (Enzyme, interacts, Enzyme).edge_index               │
    │  ├── (TF, interacts, Enzyme).edge_index                   │
    │  └── ...                                                  │
    └────────────────────────────────────────────────────────────┘
    """
    pass
```

### 5.3 Tiered Bipartite Builder

**파일 위치**: [tiered_bipartite_builder.py](file:///data/ethylene/src/tiered_bipartite_builder.py)

```python
def build_tiered_bipartite_graph(graph_path, output_path, kegg_dir):
    """
    2-Tier 증거 구조의 이종 그래프 구축
    
    Tier Structure:
    ┌─────────────────────────────────────────────────────────────┐
    │  Tier-R (weight=1.0): Reaction-grounded edges (KEGG)       │
    │  - 직접적인 효소-대사체 반응 관계                            │
    │  - 높은 신뢰도                                              │
    ├─────────────────────────────────────────────────────────────┤
    │  Tier-P (weight=0.5): Pathway-supported edges              │
    │  - 동일 경로 멤버십 기반                                     │
    │  - 중간 신뢰도                                              │
    └─────────────────────────────────────────────────────────────┘
    
    Edge Weighting:
    - edge_weight: 학습 시 loss 가중치로 사용 가능
    - edge_tier: 분석용 레이블 ('R' or 'P')
    """
    pass
```

### 5.4 그래프 스키마

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        HETEROGENEOUS GRAPH SCHEMA                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   NODE TYPES:                                                               │
│   ┌─────────────┬──────────────────────┬─────────┬──────────┐              │
│   │ Type        │ Description          │ Dim     │ Count    │              │
│   ├─────────────┼──────────────────────┼─────────┼──────────│              │
│   │ Enzyme      │ 생합성 효소          │ 64      │ ~500     │              │
│   │ Metabolite  │ 대사체               │ 64      │ ~200     │              │
│   │ TF          │ 전사인자             │ 64      │ ~100     │              │
│   │ Signaling   │ 신호전달 단백질       │ 64      │ ~50      │              │
│   │ Protein     │ 기타 단백질          │ 64      │ ~1000    │              │
│   └─────────────┴──────────────────────┴─────────┴──────────┘              │
│                                                                             │
│   EDGE TYPES:                                                               │
│   ┌─────────────────────────────────────┬────────────────────────────────┐  │
│   │ Edge Type                           │ Description                    │  │
│   ├─────────────────────────────────────┼────────────────────────────────│  │
│   │ (Enzyme, catalyzes, Metabolite)     │ 효소-대사체 촉매 관계 [예측타겟]│  │
│   │ (TF, regulates, Enzyme)             │ 전사인자-효소 조절              │  │
│   │ (Enzyme, interacts, Enzyme)         │ 효소-효소 상호작용 (PPI)       │  │
│   │ (TF, interacts, TF)                 │ 전사인자 상호작용              │  │
│   │ (Metabolite, rev_catalyzes, Enzyme) │ 역방향 (message passing용)     │  │
│   └─────────────────────────────────────┴────────────────────────────────┘  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 6. 학습 파이프라인

### 6.1 기본 Trainer

**파일 위치**: [trainer.py](file:///data/ethylene/src/trainer.py)

```python
def train():
    """
    기본 학습 루프
    
    Pipeline:
    1. 그래프 로드 및 전처리
    2. 데이터 분할 (Train/Val/Test)
    3. 모델 초기화
    4. 학습 루프 (Negative Sampling + BCE Loss)
    5. 검증 (AUC-ROC)
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # 1. Load Graph
    data = torch.load('data/processed/graph.pt')
    data = T.ToUndirected()(data)
    data = T.AddSelfLoops()(data)
    data = data.to(device)
    
    # 2. Data Split
    transform = T.RandomLinkSplit(
        num_val=0.1,
        num_test=0.1,
        is_undirected=True,
        edge_types=[target_edge_type],
        add_negative_train_samples=False
    )
    train_data, val_data, test_data = transform(data)
    
    # 3. Model Init
    model = HGT(
        metadata=data.metadata(),
        in_channels=64,
        hidden_channels=64,
        out_channels=64,
        num_heads=4,
        num_layers=2
    ).to(device)
    predictor = LinkPredictor(64).to(device)
    
    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()), 
        lr=0.01
    )
    
    # 4. Training Loop
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()
        
        # Forward pass
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        # Positive Scores
        pos_out = predictor(x_dict[src_type], x_dict[dst_type], 
                           edge_label_index)
        
        # Negative Sampling
        neg_dst_idx = torch.randint(0, x_dict[dst_type].size(0), 
                                    (num_positives,), device=device)
        neg_out = predictor(x_dict[src_type], x_dict[dst_type], 
                           torch.stack([src_idx, neg_dst_idx]))
        
        # Binary Cross Entropy Loss
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() \
               -torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        
        loss.backward()
        optimizer.step()
```

### 6.2 Refined Trainer (Node-Disjoint Split)

**파일 위치**: [refined_trainer.py](file:///data/ethylene/src/refined_trainer.py)

```python
def train_refined(graph_path):
    """
    개선된 학습 - Node-Disjoint Split
    
    핵심 차이점:
    - 10% 효소를 테스트용으로 완전히 분리
    - 새로운 효소 발견 시나리오 시뮬레이션
    - Hard Negative Sampling ("같은 경로, 다른 대사체")
    
    Evaluation:
    - Hits@20: Top-20 예측 중 정답 포함 비율
    - MAP: Mean Average Precision
    """
    
    # Node-Disjoint Split
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool)
    train_enz_mask[indices[:split]] = True
    test_enz_mask = ~train_enz_mask
    
    # Hard Negative Sampling
    # 같은 경로의 인접 대사체를 negative로 사용
    neg_src = pos_edge_index[0]  # Same enzyme
    offset = torch.randint(1, 6, (num_pos,)) * (2 * torch.randint(0, 2, (num_pos,)) - 1)
    neg_dst = (pos_edge_index[1] + offset) % num_metabolites
```

### 6.3 Loss Function

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                      BINARY CROSS ENTROPY LOSS                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Positive Loss:                                                            │
│   L_pos = -log(σ(score_positive))                                          │
│   → 양성 엣지의 점수를 높이도록 학습                                         │
│                                                                             │
│   Negative Loss:                                                            │
│   L_neg = -log(1 - σ(score_negative))                                      │
│   → 음성 엣지의 점수를 낮추도록 학습                                         │
│                                                                             │
│   Total Loss:                                                               │
│   L = mean(L_pos) + mean(L_neg)                                            │
│                                                                             │
│   where σ = sigmoid function                                                │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 7. 추론 파이프라인

### 7.1 추론 스크립트

**파일 위치**: [inference.py](file:///data/ethylene/src/inference.py)

```python
def run_inference(graph_path, model_path, output_dir):
    """
    학습된 모델로 신규 상호작용 예측
    
    Pipeline:
    1. 그래프 및 모델 로드
    2. Forward pass → 임베딩 계산
    3. Verified TF 필터링
    4. 모든 TF-Enzyme 쌍 점수 계산
    5. 기존 엣지 마스킹 (Novel만 출력)
    6. Top-K 예측 저장
    """
    
    # 1. Load Model
    model = HGT(metadata, 64, 64, 64, num_heads=4, num_layers=2)
    model.load_state_dict(torch.load(model_path))
    model.eval()
    
    # 2. Forward Pass
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        tf_emb = x_dict['TF']       # [N_TF × 64]
        enz_emb = x_dict['Enzyme']  # [N_Enz × 64]
    
    # 3. Score Matrix
    scores = (tf_emb @ enz_emb.t()).sigmoid()  # [N_TF × N_Enz]
    
    # 4. Mask Existing Edges
    for (tf_idx, enz_idx) in existing_edges:
        scores[tf_idx, enz_idx] = 0
    
    # 5. Top-K Predictions
    topk_vals, topk_indices = torch.topk(scores.flatten(), 50)
    
    # 6. Save Results
    # → top_novel_pairs.tsv
    # → top1_novel_pair.md (Case Study)
```

### 7.2 출력 형식

```
# top_novel_pairs.tsv
Rank  TF_ID           TF_Name    Enzyme_ID       Enzyme_Name    Score
1     3847.WRKY_001   WRKY15     3847.IFS_042    IFS1           0.9234
2     3847.MYB_023    MYB12      3847.CHI_018    CHI            0.8912
...

# top1_novel_pair.md
## Case Study: Top-1 Novel Verified TF-Enzyme Pair

**Pair**: WRKY15 (TF) -- IFS1 (Enzyme)
**Score**: 0.9234 (Rank #1)

## Biological Hypothesis
This high-confidence link proposes a direct regulatory axis...
```

---

## 8. 전체 파이프라인 흐름

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     COMPLETE PIPELINE WORKFLOW                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   STEP 1: Data Preparation                                                   │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  python src/dataloader.py                                             │ │
│   │  → STRING DB PPI 로드 (3847.protein.*.gz)                             │ │
│   │  → data/processed/string_interactions.pkl                             │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                              ↓                                              │
│   STEP 2: Graph Construction                                                 │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  python src/graph_builder.py --threshold 700                          │ │
│   │  → 노드 타입 분류 (Enzyme, TF, Signaling, Protein)                    │ │
│   │  → data/processed/graph.pt                                            │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                              ↓                                              │
│   STEP 3: Bipartite Graph Extension                                          │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  python src/tiered_bipartite_builder.py                               │ │
│   │  → KEGG 대사체 노드 추가                                              │ │
│   │  → Tier-R/P 엣지 생성                                                 │ │
│   │  → data/processed/tiered_bipartite_graph.pt                           │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                              ↓                                              │
│   STEP 4: Model Training                                                     │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  python src/refined_trainer.py --graph data/processed/bipartite.pt   │ │
│   │  → HGT 모델 학습 (20 epochs)                                          │ │
│   │  → Hits@20, MAP 평가                                                  │ │
│   │  → data/models/refined_hgt.pth                                        │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                              ↓                                              │
│   STEP 5: Inference & Prediction                                             │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  python src/inference.py                                              │ │
│   │  → Top-50 신규 상호작용 예측                                          │ │
│   │  → results/case_study/top_novel_pairs.tsv                             │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                              ↓                                              │
│   STEP 6: Explainability (Optional)                                          │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  python src/gnnshap_explainability.py                                 │ │
│   │  → Shapley Value 기반 엣지 중요도 분석                                 │ │
│   │  → results/explainability/                                            │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                              ↓                                              │
│   STEP 7: Structural Validation                                              │
│   ┌───────────────────────────────────────────────────────────────────────┐ │
│   │  python src/run_docking.py                                            │ │
│   │  → AutoDock Vina 분자 도킹                                            │ │
│   │  → results/docking/                                                   │ │
│   └───────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 📊 평가 지표 요약

| 지표 | 설명 | 목표 |
|------|------|------|
| **AUC-ROC** | ROC 곡선 아래 면적 | ≥ 0.85 |
| **Hits@20** | Top-20 예측 중 정답 비율 | ≥ 0.50 |
| **MAP** | Mean Average Precision | ≥ 0.30 |
| **Fidelity** | 설명 신뢰도 (GNNShap) | > 0.30 |

---

## 📁 핵심 파일 목록

| 파일 | 경로 | 역할 |
|------|------|------|
| **model.py** | [src/model.py](file:///data/ethylene/src/model.py) | 모델 정의 (HGT, HAN, LinkPredictor) |
| **trainer.py** | [src/trainer.py](file:///data/ethylene/src/trainer.py) | 기본 학습 루프 |
| **refined_trainer.py** | [src/refined_trainer.py](file:///data/ethylene/src/refined_trainer.py) | 개선된 학습 (Node-Disjoint) |
| **dataloader.py** | [src/dataloader.py](file:///data/ethylene/src/dataloader.py) | STRING DB 데이터 로더 |
| **graph_builder.py** | [src/graph_builder.py](file:///data/ethylene/src/graph_builder.py) | 이종 그래프 구축 |
| **tiered_bipartite_builder.py** | [src/tiered_bipartite_builder.py](file:///data/ethylene/src/tiered_bipartite_builder.py) | 2-Tier 그래프 확장 |
| **inference.py** | [src/inference.py](file:///data/ethylene/src/inference.py) | 추론 및 예측 |

---

*문서 작성: 2026년 1월 22일*
