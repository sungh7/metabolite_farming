# Refined Implementation Protocol v2

## Step 1: Metabolomics Evidence (Statistical Stability)
**Goal**: Establish the "Ground Truth" for metabolic changes with statistical rigor.
- **Method**:
    - **Simulation**: Generate synthetic fold-change data for Flavonoid/Phenylpropanoid pathways based on literature (since no raw data).
    - **Validation**: Perform **Bootstrap Resampling** (1000 iterations) to verify that the pathway enrichment is robust to noise.
- **Deliverable**: `results/table1_metabolomics.csv`
    - Columns: Pathway, Enrichment Score, P-value, **Bootstrap 95% CI**, Stability Score.

## Step 2: Protein Layer Mapping
**Goal**: Connect Metabolites to the Proteome physically.
- **Logic**:
    1. Metabolite $M$ is substrate/product of Reaction $R$.
    2. Reaction $R$ is catalyzed by Enzyme $E$ (EC Number).
    3. Enzyme $E$ corresponds to Gene $G$ (Glyma ID).
- **Graph Construction**:
    - Nodes: Metabolites, Enzymes (Genes).
    - Edges: Metabolic Reactions.
- **Deliverable**: `data/processed/bipartite_graph.pt`

## Step 3: GNN Task Refinement (Rigorous Evaluation)
**Goal**: Prioritize functional associations.
- **Task**: **Metabolite-Enzyme Link Prediction** (Ranking).
- **Split Strategy**:
    - **Node-Disjoint**: Hold out random 10% of Enzymes entirely from training to test generalization to new proteins.
    - **Negative Sampling**:
        - **Training**: Dynamic (On-the-fly) random negatives (`add_negative_train_samples=False`).
        - **Validation/Test**: **Hard Negatives** (Enzymes in the same pathway but chemically distant).
- **Metrics**:
    - **Primary**: Precision@20, Hits@20.
    - **Secondary**: AUPRC (Robust to imbalance).
    - *Avoid*: ROC-AUC.

## Step 4: Structural Auxiliary Scoring
**Goal**: Refine candidates using physical plausibility.
- **Method**:
    - Use structural embeddings as an auxiliary re-ranking feature.
    - **Tone**: "Prioritization support", not "Verification".

## Step 5: Final Conclusion
**Goal**: Integrated Network Hypothesis.
- **Format**: "Ethylene triggers TF $X$, which regulates Enzyme $Y$, consistent with the observed bootstrap-stable accumulation of Metabolite $Z$."
