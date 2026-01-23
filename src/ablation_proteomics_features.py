"""
Proteomics Feature Ablation Study

Compares GNN performance with different node feature configurations:
1. Random (baseline): Learnable random embeddings
2. Proteomics: Log2FC + statistical features as initial node features
3. Combined: Proteomics features concatenated with learnable embeddings

Hypothesis: Incorporating expression data should improve predictions.
"""

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import os
import sys
sys.path.append(os.getcwd())

from torch_geometric.nn import HGTConv, Linear


def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


class HGTWithFeatures(nn.Module):
    """HGT model that can accept external node features."""

    def __init__(self, metadata, in_channels, hidden_channels, out_channels,
                 num_heads=4, num_layers=3, feature_dims=None):
        super().__init__()

        self.feature_dims = feature_dims or {}

        # Input projection: handle different feature dimensions per node type
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            feat_dim = self.feature_dims.get(node_type, in_channels)
            self.lin_dict[node_type] = Linear(feat_dim, hidden_channels)

        # HGT convolutions
        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HGTConv(hidden_channels, hidden_channels, metadata, heads=num_heads)
            self.convs.append(conv)

    def forward(self, x_dict, edge_index_dict):
        # Project to hidden dimension
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
        }

        # Message passing
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)

        return x_dict


class LinkPredictor(nn.Module):
    def __init__(self, in_channels):
        super().__init__()

    def forward(self, x_src, x_dst, edge_label_index):
        row, col = edge_label_index
        return (x_src[row] * x_dst[col]).sum(dim=-1)


def load_proteomics_features(proteomics_path, enzyme_mapping_path):
    """
    Load proteomics data and create feature matrix for enzymes.

    Returns:
        enzyme_features: dict mapping uniprot_id -> feature vector
        feature_names: list of feature names
    """
    # Load data
    prot_df = pd.read_csv(proteomics_path)
    mapping_df = pd.read_csv(enzyme_mapping_path)

    # Create STRING_ID to UniProt mapping
    string_to_uniprot = {}
    for _, row in mapping_df.iterrows():
        string_id = f"3847.{row['uniprot_id']}"
        string_to_uniprot[string_id] = row['uniprot_id']

    # Extract features
    feature_names = ['Log2FC', 'P_Value', 'Mean_Control', 'Mean_Ethylene']
    enzyme_features = {}

    for _, row in prot_df.iterrows():
        string_id = row['STRING_ID']
        if pd.isna(string_id):
            continue

        if string_id in string_to_uniprot:
            uniprot_id = string_to_uniprot[string_id]

            features = []
            for feat in feature_names:
                val = row[feat] if feat in row and not pd.isna(row[feat]) else 0.0
                features.append(val)

            # Add derived features
            features.append(abs(row['Log2FC']))  # Absolute fold change
            features.append(-np.log10(row['P_Value'] + 1e-10))  # -log10(p)
            features.append(1 if row['P_Value'] < 0.05 else 0)  # Significance flag

            enzyme_features[uniprot_id] = features

    feature_names_extended = feature_names + ['Abs_Log2FC', 'NegLog10P', 'Significant']

    return enzyme_features, feature_names_extended


def create_feature_matrix(data, enzyme_features, mapping_df, feature_dim=64, mode='random'):
    """
    Create node feature matrices based on mode.

    Modes:
    - 'random': Random learnable embeddings (baseline)
    - 'proteomics': Proteomics features only
    - 'combined': Proteomics + random embeddings
    """
    device = data['Enzyme'].x.device
    num_enzymes = data['Enzyme'].num_nodes

    if mode == 'random':
        # Keep original random features
        return data

    # Build mapping: enzyme_idx -> uniprot_id
    idx_to_uniprot = {}
    for _, row in mapping_df.iterrows():
        idx_to_uniprot[int(row['enzyme_idx'])] = row['uniprot_id']

    # Create feature matrix
    n_prot_features = 7  # Log2FC, P_Value, Mean_Control, Mean_Ethylene, Abs_FC, NegLogP, Sig

    if mode == 'proteomics':
        # Proteomics features only (padded to feature_dim)
        enzyme_x = torch.zeros(num_enzymes, feature_dim, device=device)

        matched = 0
        for idx in range(num_enzymes):
            if idx in idx_to_uniprot:
                uniprot_id = idx_to_uniprot[idx]
                if uniprot_id in enzyme_features:
                    feats = enzyme_features[uniprot_id]
                    enzyme_x[idx, :len(feats)] = torch.tensor(feats, dtype=torch.float32)
                    matched += 1

        print(f"  Proteomics features matched: {matched}/{num_enzymes} enzymes")

        # Normalize
        scaler_mean = enzyme_x.mean(dim=0)
        scaler_std = enzyme_x.std(dim=0) + 1e-8
        enzyme_x = (enzyme_x - scaler_mean) / scaler_std

        data['Enzyme'].x = enzyme_x

    elif mode == 'combined':
        # Concatenate proteomics + random
        random_x = data['Enzyme'].x  # Original 64-dim

        prot_x = torch.zeros(num_enzymes, n_prot_features, device=device)

        matched = 0
        for idx in range(num_enzymes):
            if idx in idx_to_uniprot:
                uniprot_id = idx_to_uniprot[idx]
                if uniprot_id in enzyme_features:
                    feats = enzyme_features[uniprot_id]
                    prot_x[idx, :len(feats)] = torch.tensor(feats, dtype=torch.float32)
                    matched += 1

        print(f"  Combined features matched: {matched}/{num_enzymes} enzymes")

        # Normalize proteomics features
        scaler_mean = prot_x.mean(dim=0)
        scaler_std = prot_x.std(dim=0) + 1e-8
        prot_x = (prot_x - scaler_mean) / scaler_std

        # Concatenate
        combined_x = torch.cat([random_x, prot_x], dim=1)
        data['Enzyme'].x = combined_x

    return data


def create_splits(num_nodes, seed=42):
    """Create train/val/test splits."""
    indices = np.arange(num_nodes)
    train_idx, temp_idx = train_test_split(indices, train_size=0.7, random_state=seed)
    val_idx, test_idx = train_test_split(temp_idx, train_size=0.5, random_state=seed)
    return train_idx, val_idx, test_idx


def train_epoch(model, predictor, optimizer, data, train_edges, device):
    """Single training epoch."""
    model.train()
    optimizer.zero_grad()

    x_dict = model(data.x_dict, data.edge_index_dict)

    num_pos = train_edges.size(1)
    num_metabolites = data['Metabolite'].num_nodes
    neg_src = train_edges[0]
    offset = torch.randint(1, 6, (num_pos,), device=device) * \
             (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
    neg_dst = (train_edges[1] + offset) % num_metabolites

    pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], train_edges)
    neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'],
                       torch.stack([neg_src, neg_dst]))

    loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - \
           torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()

    loss.backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
    optimizer.step()

    return loss.item()


def evaluate(model, predictor, data, edges, candidate_nodes, device):
    """Evaluate on edges."""
    model.eval()

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']

        src, dst = edges
        unique_mets = torch.unique(dst)

        hits_10, hits_20 = 0, 0
        reciprocal_ranks = []
        count = 0

        for met_idx in unique_mets:
            met_mask = (dst == met_idx)
            true_enzs = src[met_mask]

            num_cands = candidate_nodes.size(0)
            eval_edges = torch.stack([candidate_nodes, met_idx.repeat(num_cands)])

            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()
            is_true = torch.isin(candidate_nodes, true_enzs)

            if is_true.sum() == 0:
                continue

            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]

            if sorted_labels[:10].sum() > 0:
                hits_10 += 1
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1

            true_pos = torch.where(sorted_labels)[0]
            if len(true_pos) > 0:
                reciprocal_ranks.append(1.0 / (true_pos[0].item() + 1))

            count += 1

    return {
        'hits@10': hits_10 / count if count > 0 else 0.0,
        'hits@20': hits_20 / count if count > 0 else 0.0,
        'mrr': np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0
    }


def run_experiment(base_data, enzyme_features, mapping_df, mode, seed, device,
                   num_layers=3, epochs=50, patience=10):
    """Run single experiment."""
    set_seed(seed)

    # Clone and add features
    data = base_data.clone()
    data = create_feature_matrix(data, enzyme_features, mapping_df, mode=mode)

    num_enzymes = data['Enzyme'].num_nodes
    train_idx, val_idx, test_idx = create_splits(num_enzymes, seed=seed)

    train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    val_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    test_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_mask[train_idx] = True
    val_mask[val_idx] = True
    test_mask[test_idx] = True

    # Split edges
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]

    train_edges = edge_index[:, train_mask[src]]
    val_edges = edge_index[:, val_mask[src]]
    test_edges = edge_index[:, test_mask[src]]

    # Create training graph
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges

    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, train_mask[rev_index[1]]]

    # Determine feature dimensions
    feature_dims = {}
    for node_type in train_data.node_types:
        feature_dims[node_type] = train_data[node_type].x.shape[1]

    # Model
    model = HGTWithFeatures(
        train_data.metadata(),
        in_channels=64,
        hidden_channels=64,
        out_channels=64,
        num_heads=4,
        num_layers=num_layers,
        feature_dims=feature_dims
    ).to(device)

    predictor = LinkPredictor(64).to(device)

    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=0.01, weight_decay=1e-5
    )

    # Training with early stopping
    best_val_mrr = 0
    best_state = None

    train_nodes = torch.where(train_mask)[0]
    val_nodes = torch.where(val_mask)[0]
    test_nodes = torch.where(test_mask)[0]

    patience_counter = 0

    for epoch in range(1, epochs + 1):
        loss = train_epoch(model, predictor, optimizer, train_data, train_edges, device)

        if epoch % 5 == 0:
            val_metrics = evaluate(model, predictor, train_data, val_edges, val_nodes, device)

            if val_metrics['mrr'] > best_val_mrr:
                best_val_mrr = val_metrics['mrr']
                patience_counter = 0
                best_state = {
                    'model': {k: v.cpu().clone() for k, v in model.state_dict().items()},
                    'predictor': {k: v.cpu().clone() for k, v in predictor.state_dict().items()}
                }
            else:
                patience_counter += 1

            if patience_counter >= patience // 5:
                break

    # Load best and evaluate
    if best_state:
        model.load_state_dict({k: v.to(device) for k, v in best_state['model'].items()})
        predictor.load_state_dict({k: v.to(device) for k, v in best_state['predictor'].items()})

    test_metrics = evaluate(model, predictor, train_data, test_edges, test_nodes, device)

    return test_metrics


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Load data
    base_data = torch.load('data/processed/strict_bipartite_v2.pt')
    base_data = base_data.to(device)

    # Load proteomics features
    print("\nLoading proteomics features...")
    enzyme_features, feature_names = load_proteomics_features(
        'data/processed/pxd006989_mapped.csv',
        'data/processed/enzyme_string_mapping.csv'
    )
    print(f"  Features: {feature_names}")
    print(f"  Enzymes with proteomics data: {len(enzyme_features)}")

    # Load mapping
    mapping_df = pd.read_csv('data/processed/enzyme_string_mapping.csv')

    seeds = [42, 123, 456, 789, 1024]
    modes = ['random', 'proteomics', 'combined']

    all_results = {mode: [] for mode in modes}

    print("\n" + "=" * 80)
    print("PROTEOMICS FEATURE ABLATION")
    print("=" * 80)

    for mode in modes:
        print(f"\n--- Mode: {mode} ---")

        for seed in seeds:
            print(f"  Seed {seed}...", end=" ", flush=True)

            result = run_experiment(
                base_data, enzyme_features, mapping_df,
                mode=mode, seed=seed, device=device,
                num_layers=3, epochs=50
            )

            all_results[mode].append(result)
            print(f"H@20={result['hits@20']:.3f}, MRR={result['mrr']:.4f}")

    # Summary
    print("\n" + "=" * 80)
    print("PROTEOMICS FEATURE ABLATION RESULTS")
    print("=" * 80)

    print(f"\n{'Mode':<12} {'Hits@10':<18} {'Hits@20':<18} {'MRR':<18}")
    print("-" * 70)

    for mode in modes:
        results = all_results[mode]
        h10 = np.array([r['hits@10'] for r in results])
        h20 = np.array([r['hits@20'] for r in results])
        mrr = np.array([r['mrr'] for r in results])

        print(f"{mode:<12} "
              f"{h10.mean()*100:.2f}±{h10.std()*100:.2f}%    "
              f"{h20.mean()*100:.2f}±{h20.std()*100:.2f}%    "
              f"{mrr.mean():.4f}±{mrr.std():.4f}")

    # Improvement analysis
    print("\n" + "=" * 80)
    print("IMPROVEMENT vs RANDOM BASELINE")
    print("=" * 80)

    baseline_h20 = np.mean([r['hits@20'] for r in all_results['random']])
    baseline_mrr = np.mean([r['mrr'] for r in all_results['random']])

    for mode in ['proteomics', 'combined']:
        h20 = np.mean([r['hits@20'] for r in all_results[mode]])
        mrr = np.mean([r['mrr'] for r in all_results[mode]])

        h20_delta = (h20 - baseline_h20) / baseline_h20 * 100 if baseline_h20 > 0 else 0
        mrr_delta = (mrr - baseline_mrr) / baseline_mrr * 100 if baseline_mrr > 0 else 0

        print(f"{mode}: Hits@20 {h20_delta:+.1f}%, MRR {mrr_delta:+.1f}%")

    # Save results
    os.makedirs('results/gnn', exist_ok=True)

    with open('results/gnn/proteomics_feature_ablation.tsv', 'w') as f:
        f.write("Mode\tHits10_Mean\tHits10_Std\tHits20_Mean\tHits20_Std\tMRR_Mean\tMRR_Std\n")
        for mode in modes:
            results = all_results[mode]
            h10 = np.array([r['hits@10'] for r in results])
            h20 = np.array([r['hits@20'] for r in results])
            mrr = np.array([r['mrr'] for r in results])
            f.write(f"{mode}\t{h10.mean():.4f}\t{h10.std():.4f}\t"
                    f"{h20.mean():.4f}\t{h20.std():.4f}\t"
                    f"{mrr.mean():.4f}\t{mrr.std():.4f}\n")

    print("\nResults saved to: results/gnn/proteomics_feature_ablation.tsv")


if __name__ == "__main__":
    main()
