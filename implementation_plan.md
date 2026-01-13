# Implementation Plan: Ethylene-Metabolite Interaction Prediction

This plan targets the development of a Graph Neural Network (HAN) to predict key regulators in the Ethylene-Isoflavonoid pathway in Soybean.

## Goal
To build a computational pipeline that integrates heterogeneous biological data (PPI, GRN, Metabolic pathways) and learns to predict protein-metabolite functional associations using a Heterogeneous Attention Network (HAN).

## User Review Required
> [!IMPORTANT]
> **Data Availability Claim**: This plan assumes simultaneous access to public databases (SoyNet, KEGG) and future experimental data. Phase 1 focuses on **Public Data** to ensure the code is functional before wet-lab data arrives.

## Proposed Changes

We will create a structured Python project.

### Directory Structure
```
project_root/
├── data/
│   ├── raw/                # SoyNet, STRING, KEGG files
│   ├── processed/          # PyG graph objects
│   └── experimental/       # Future LC-MS/RNA-seq data
├── src/
│   ├── __init__.py
│   ├── dataloader.py       # Parsers for SoyNet/KEGG/PlantTFDB
│   ├── graph_builder.py    # Construct HeteroData for PyG
│   ├── model.py            # HAN (Heterogeneous Attention Network) implementation
│   ├── trainer.py          # Training loop (Link Prediction)
│   └── utils.py            # Metrics and visualization
├── notebooks/              # Exploratory Data Analysis (EDA)
└── requirements.txt
```

### Component Breakdown

#### [NEW] Data Ingestion Layer (`src/dataloader.py`)
- **Functionality**:
    - **STRING DB**: Parse `3847.protein.info.v12.0.txt.gz` (Metadata) and `3847.protein.links.full.v12.0.txt.gz` (PPI Network).
    - **SoyNet**: Download/Parse Gene-Gene associations.
    - **KEGG**: Parse pathway maps (Enzyme-Metabolite).
    - **PlantTFDB**: Parse TF-Target interactions.
- **Key Task**: Map all STRING protein IDs (e.g., `3847.GLYMA_...`) to a common space (Glyma IDs).

#### [NEW] Graph Construction (`src/graph_builder.py`)
- **Functionality**:
    - Instantiate `torch_geometric.data.HeteroData`.
    - Define Node Types: `['TF', 'Signaling', 'Enzyme', 'Metabolite']`.
    - Define Edge Types: `('TF', 'regulates', 'Enzyme')`, `('Enzyme', 'produces', 'Metabolite')`, etc.
    - **Feature Engineering**: Initialize node features using *node2vec* or random embeddings (initially), later replaced by expression profiles.

#### [NEW] Model Architecture (`src/model.py`)
- **Functionality**:
    - Implement **HAN** (Heterogeneous Graph Attention Network).
    - Use `torch_geometric.nn.HANConv`.
    - Define Semantic Attention (calculating importance of Metapaths).
    - **Metapaths**:
        1. TF $\to$ Enzyme $\to$ Metabolite
        2. Signaling $\to$ TF $\to$ Enzyme

#### [NEW] Training Pipeline (`src/trainer.py`)
- **Task**: Link Prediction (predict missing 'regulates' or 'associated_with' edges).
- **Loss**: Binary Cross Entropy with Negative Sampling.
- **Optimization**: AdamW.

## Verification Plan

### Automated Tests
- **Unit Tests**:
    - Verify Graph connectivity (no isolated nodes).
    - Test Matrix dimensions in HAN forward pass.
- **Integration Tests**:
    - Run a "dummy run" with a small synthetic subgraph to ensure the pipeline runs end-to-end without crashing.

### Manual Verification
- **Metapath Analysis**:
    - Extract attention weights from the trained HAN.
    - Check if the model prioritizes "biologically meaningful" paths (e.g., does it pay attention to the 'TF->Enzyme' path for metabolite prediction?).
- **Top-N Recommendation**:
    - Validate if known Ethylene regulators (EIN3, ERF) appear in the top-predicted list for Isoflavonoids.
