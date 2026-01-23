# Graph Design v3: Academically Defensible Architecture

This document describes the v3 graph design for the Ethylene-Isoflavonoid GNN project, which addresses key academic defensibility concerns.

## Overview

The v3 design separates **supervision edges** (Tier-R) from **context edges** (rxn_neighbor), enabling cleaner evaluation and preventing information leakage.

### Key Changes from v2

| Aspect | v2 | v3 |
|--------|-----|-----|
| **Tier-P Edge Type** | Enzyme→Metabolite | **Metabolite↔Metabolite** |
| **Tier-P Semantics** | EC-class similarity | **Reaction distance ≤2** |
| **Supervision** | Tier-R + Tier-P (weight 0.5) | **Tier-R only** |
| **Node Features** | `torch.randn` (fixed) | **nn.Embedding + omics** |
| **TF Edge Semantics** | `interacts` (ambiguous) | **`associates`** (functional) |
| **Evaluation** | Raw ranking | **Filtered ranking** |

---

## Graph Structure

### Node Types

1. **Enzyme** (~3,500 nodes)
   - Features: Learnable embedding (61d) + proteomics (3d: log2FC, -log10(p), abundance)
   - Source: STRING database (soybean PPI)

2. **Metabolite** (~306 nodes)
   - Features: Learnable embedding (61d) + metabolomics (3d: log2FC, -log10(p), n_pathways)
   - Source: KEGG compound database + MTBLS531 experimental data

3. **TF** (Transcription Factors, ~500 nodes)
   - Features: Learnable embedding (58d) + domain one-hot (6d: ERF, WRKY, MYB, bHLH, zinc_finger, other)
   - Source: STRING database annotations

### Edge Types

| Edge Type | Source | Target | Semantics | Use |
|-----------|--------|--------|-----------|-----|
| `ppi` | Enzyme | Enzyme | Protein-protein interaction | Message passing (all layers) |
| `catalyzes_R` | Enzyme | Metabolite | Direct enzymatic reaction (Tier-R) | **Supervision + Message passing** |
| `rev_catalyzes_R` | Metabolite | Enzyme | Reverse of catalyzes_R | Message passing |
| `rxn_neighbor` | Metabolite | Metabolite | Reaction distance ≤2 | **Message passing (Layer 1 only)** |
| `associates` | TF | Enzyme | Functional association | Message passing (Layer 1 only) |
| `rev_associates` | Enzyme | TF | Reverse of associates | Message passing (Layer 1 only) |

---

## Tier-R (Supervision Edges)

### Definition
Tier-R edges represent **direct enzymatic reactions** from KEGG:
- An enzyme (via EC number) catalyzes the production or consumption of a metabolite
- Source: KEGG `full_enzyme_metabolite_edges.tsv`

### Properties
- **Edge weight**: 1.0 (full confidence)
- **Usage**: Supervision (link prediction target) + message passing
- **Statistics**: ~1,500-2,000 edges

### Construction
```python
for each KEGG reaction (EC, metabolite):
    for enzyme with EC mapping:
        tier_r_edges.add((enzyme_idx, metabolite_idx))
```

---

## rxn_neighbor (Metabolite Context Edges)

### Definition
`rxn_neighbor` represents **reaction neighborhood proximity**:
- Distance 1: Metabolites participate in the same reaction
- Distance 2: Metabolites share a common neighbor (via BFS)

### Design Rationale

> "rxn_neighbor is NOT a directional flux model. It represents reaction neighborhood proximity as a structural prior. The directionality of metabolic flux is not relevant to this task, hence edges are undirected."

### Properties
- **Undirected**: Both (A, B) and (B, A) edges exist
- **No self-loops**: Validated at build time
- **Currency exclusion**: ATP, H2O, etc. are excluded from both keys and neighbors
- **Statistics**: ~2,500 edges (10x reduction from v2 Tier-P)

### Construction
```python
# Step 1: Build reaction graph
for each reaction:
    mets = metabolites_in_reaction (excluding currency)
    for all pairs (met_i, met_j) in mets:
        reaction_adj[met_i].add(met_j)

# Step 2: Extend to distance 2
for met in reaction_adj:
    for neighbor in reaction_adj[met]:
        for hop2 in reaction_adj[neighbor]:
            if hop2 != met:
                edges.add((met, hop2))
```

### Information Leakage Control

**Problem**: rxn_neighbor edges could allow the model to learn metabolite community structure rather than enzyme function.

**Solution**:
1. **Ablation requirement**: Results table must include both `rxn_neighbor ON` and `rxn_neighbor OFF`
2. **Edge dropout**: 30% of rxn_neighbor edges randomly dropped during training
3. **Layer restriction**: rxn_neighbor only used in Layer 1 (not multi-hop)

```python
V3_CONFIG = {
    'rxn_neighbor_dropout': 0.3,
    'rxn_neighbor_dropout_mode': 'per_run_fixed',  # reproducibility
    'layer1_edge_types': [..., 'rxn_neighbor', ...],
    'layer2_edge_types': [... (no rxn_neighbor) ...],
}
```

---

## Currency Metabolite Handling

### Excluded Metabolites
```python
CURRENCY_METABOLITES = {
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
```

### Validation
Currency metabolites are excluded from:
1. `reaction_adj` keys (source metabolites)
2. `reaction_adj` values (neighbor metabolites)
3. `rxn_neighbor` edge index (both src and dst)

```python
def validate_no_currency(reaction_adj, met_to_idx, rxn_neighbor_edge_index):
    for met_id, neighbors in reaction_adj.items():
        assert met_id not in CURRENCY_METABOLITES
        for neighbor in neighbors:
            assert neighbor not in CURRENCY_METABOLITES

    for idx in rxn_neighbor_edge_index.flatten():
        assert idx_to_met[idx] not in CURRENCY_METABOLITES
```

---

## Node Features

### Design Philosophy
- **Learnable embeddings** are model parameters (nn.Embedding)
- **Omics features** are fixed input data (stored in HeteroData)
- **Linear concat** (no MLP projection): omics features are few, MLP would over-parameterize

### Metabolite Features
```python
# Learnable: [num_mets, 61]
met_embedding = nn.Embedding(num_mets, 61)

# Omics: [num_mets, 3] (z-score normalized)
# - log2FC: Differential expression
# - -log10(pvalue): Statistical significance
# - n_pathways: Pathway annotation count

met_x = torch.cat([met_embedding(node_ids), omics_x], dim=-1)  # [num_mets, 64]
```

### Enzyme Features
```python
# Learnable: [num_enz, 61]
enz_embedding = nn.Embedding(num_enz, 61)

# Omics: [num_enz, 3] (z-score normalized)
# - log2FC: Proteomics differential
# - -log10(pvalue): Statistical significance
# - log1p(abundance): Expression level

enz_x = torch.cat([enz_embedding(node_ids), omics_x], dim=-1)  # [num_enz, 64]
```

### TF Features
```python
# Learnable: [num_tfs, 58]
tf_embedding = nn.Embedding(num_tfs, 58)

# Domain one-hot: [num_tfs, 6]
# Categories: ERF, WRKY, MYB, bHLH, zinc_finger, other

tf_x = torch.cat([tf_embedding(node_ids), domain_x], dim=-1)  # [num_tfs, 64]
```

---

## Layer-wise Edge Type Filtering

### Rationale
- **Layer 1**: All edge types provide local context
- **Layer 2+**: Only strong relationships (PPI, catalyzes_R) to prevent noise amplification

> "rxn_neighbor and TF associates are weak associations. Multi-hop propagation through weak edges causes noise diffusion and representation pollution."

### Configuration
```python
V3_CONFIG = {
    'layer1_edge_types': [
        ('Enzyme', 'ppi', 'Enzyme'),
        ('Enzyme', 'catalyzes_R', 'Metabolite'),
        ('Metabolite', 'rev_catalyzes_R', 'Enzyme'),
        ('Metabolite', 'rxn_neighbor', 'Metabolite'),
        ('TF', 'associates', 'Enzyme'),
        ('Enzyme', 'rev_associates', 'TF'),
    ],
    'layer2_edge_types': [
        ('Enzyme', 'ppi', 'Enzyme'),
        ('Enzyme', 'catalyzes_R', 'Metabolite'),
        ('Metabolite', 'rev_catalyzes_R', 'Enzyme'),
    ],
}
```

---

## Filtered Evaluation Protocol

### Problem
Raw ranking can be pessimistic: known positives may rank higher than the test positive, inflating rank estimates.

### Solution: Filtered Ranking
For each test edge (enzyme, metabolite):
1. Score all candidate metabolites
2. **Filter** (set to -inf) known Tier-R positives for that enzyme
3. Compute rank of the true metabolite

### Implementation
```python
def evaluate_filtered(model, test_edges, tier_r_lookup, candidate_mets):
    # Global→Local mapping for candidate subset
    global_to_local = {m: i for i, m in enumerate(candidate_mets)}

    for enz_idx, met_idx in test_edges:
        scores = model.predict(enz_idx, candidate_mets)  # [num_cand]

        # Filter known positives
        for known_met in tier_r_lookup.get(enz_idx, set()):
            if known_met != met_idx and known_met in global_to_local:
                scores[global_to_local[known_met]] = -inf

        # Compute rank
        rank = (scores > scores[global_to_local[met_idx]]).sum() + 1
```

### Two-Tier Evaluation
1. **Full (306 metabolites)**: Standard link prediction (Hits@K, MRR)
2. **Experimental (10 metabolites)**: Biological relevance (Precision, Recall)

---

## Task Definition

### Primary Task
**Enzyme → Metabolite link prediction**

Given an enzyme, rank all candidate metabolites by likelihood of catalytic relationship.

### Candidate Sets
- **Full**: All 306 metabolites in the graph
- **Experimental**: 10 metabolites from MTBLS531 with experimental validation

### Limitation
> KEGG annotation incompleteness may cause some false negatives. A predicted positive that is not in Tier-R is not necessarily wrong—it may be an unannotated true positive.

---

## Ablation Requirements

The following ablations must be included in the main results table:

| Setting | Description |
|---------|-------------|
| Full | All edge types enabled |
| No rxn_neighbor | `rxn_neighbor` edges removed |
| No TF | TF `associates` edges removed |
| No rxn_neighbor + No TF | Both removed |

### Expected Outcome
If rxn_neighbor contributes meaningfully:
- `Full` > `No rxn_neighbor` (modest improvement)
- If `Full` >> `No rxn_neighbor`, investigate potential leakage

---

## File Structure

```
src/
├── config.py              # V3_CONFIG added
├── data_pipeline_v3.py    # Graph builder (NEW)
├── model_v3.py            # HGTv3 with layer-wise filtering (NEW)
└── train_v3.py            # Filtered evaluation (NEW)

data/processed/
└── graph_v3.pt            # v3 graph output

docs/
└── GRAPH_DESIGN.md        # This document
```

---

## Usage

```bash
# 1. Build v3 graph
python src/data_pipeline_v3.py --output data/processed/graph_v3.pt

# 2. Train with full configuration
python src/train_v3.py --graph data/processed/graph_v3.pt --seeds 42,123,456

# 3. Ablation: no rxn_neighbor
python src/train_v3.py --graph data/processed/graph_v3.pt --no-rxn-neighbor

# 4. Verify graph
python -c "
import torch
d = torch.load('data/processed/graph_v3.pt')
print('Tier-R:', d['Enzyme', 'catalyzes_R', 'Metabolite'].edge_index.shape[1])
print('rxn_neighbor:', d['Metabolite', 'rxn_neighbor', 'Metabolite'].edge_index.shape[1])
print('Experimental mets:', d['Metabolite'].is_experimental.sum().item())
"
```

---

## Changelog

- **v3.0**: Initial academically defensible design
  - Separated Tier-P into Metabolite↔Metabolite (rxn_neighbor)
  - Added filtered evaluation protocol
  - Implemented layer-wise edge type filtering
  - Added omics features with learnable embeddings
  - Currency metabolite validation (key + neighbor)
