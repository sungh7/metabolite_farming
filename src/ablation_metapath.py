"""
Metapath Edge Ablation Study

Compares HGT performance with/without metapath shortcuts:
- Baseline: Original graph (TF→Enzyme→Metabolite = 2-hop)
- Metapath: Add TF→Metabolite edges computed from TF→Enzyme→Metabolite paths

Hypothesis: Metapath shortcuts should improve TF-Metabolite association learning.
"""

import torch
import numpy as np
from sklearn.model_selection import KFold
import os
import sys
sys.path.append(os.getcwd())

from src.model import HGT, LinkPredictor


def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def add_metapath_edges(data, device):
    """
    Add TF→Metabolite metapath edges.

    Computes: TF → Enzyme → Metabolite paths and creates direct edges.
    Edge type: ('TF', 'regulates', 'Metabolite')
    """
    # Get TF→Enzyme edges
    tf_enz_edges = data['TF', 'interacts', 'Enzyme'].edge_index

    # Get Enzyme→Metabolite edges
    enz_met_edges = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index

    # Build enzyme→metabolite mapping
    enz_to_mets = {}
    for i in range(enz_met_edges.shape[1]):
        enz = enz_met_edges[0, i].item()
        met = enz_met_edges[1, i].item()
        if enz not in enz_to_mets:
            enz_to_mets[enz] = []
        enz_to_mets[enz].append(met)

    # Compute TF→Metabolite edges via TF→Enzyme→Metabolite
    tf_met_edges = set()
    for i in range(tf_enz_edges.shape[1]):
        tf = tf_enz_edges[0, i].item()
        enz = tf_enz_edges[1, i].item()

        if enz in enz_to_mets:
            for met in enz_to_mets[enz]:
                tf_met_edges.add((tf, met))

    if len(tf_met_edges) == 0:
        print("  Warning: No TF→Metabolite metapath edges found")
        return data

    # Convert to tensor
    tf_met_list = list(tf_met_edges)
    tf_indices = torch.tensor([e[0] for e in tf_met_list], device=device)
    met_indices = torch.tensor([e[1] for e in tf_met_list], device=device)
    metapath_edge_index = torch.stack([tf_indices, met_indices])

    # Add new edge type
    data['TF', 'regulates', 'Metabolite'].edge_index = metapath_edge_index
    data['Metabolite', 'rev_regulates', 'TF'].edge_index = torch.stack([met_indices, tf_indices])

    print(f"  Added {len(tf_met_edges)} TF→Metabolite metapath edges")
    print(f"  Unique TFs: {len(set(tf_indices.tolist()))}")
    print(f"  Unique Metabolites: {len(set(met_indices.tolist()))}")

    return data


def train_and_evaluate(data, train_mask, test_mask, seed, device,
                       num_layers=2, epochs=30, use_metapath=False):
    """
    Train and evaluate HGT with optional metapath edges.
    """
    set_seed(seed)

    train_data = data.clone()

    # Add metapath edges if requested
    if use_metapath:
        train_data = add_metapath_edges(train_data, device)

    num_enzymes = data['Enzyme'].num_nodes

    # Mask edges for train/test split
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]

    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges

    # Handle reverse edges
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1]
    rev_mask = train_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]

    # Model with specified configuration
    model = HGT(
        train_data.metadata(),
        in_channels=64,
        hidden_channels=64,
        out_channels=64,
        num_heads=4,
        num_layers=num_layers
    ).to(device)

    predictor = LinkPredictor(64).to(device)

    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=0.01,
        weight_decay=1e-5
    )

    # Training
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()

        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        pos_edge_index = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index

        # Hard Negative Sampling
        num_pos = pos_edge_index.size(1)
        num_metabolites = train_data['Metabolite'].num_nodes
        neg_src = pos_edge_index[0]
        offset = torch.randint(1, 6, (num_pos,), device=device) * \
                 (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge_index[1] + offset) % num_metabolites

        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'],
                           torch.stack([neg_src, neg_dst]))

        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - \
               torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

    # Evaluation
    model.eval()
    test_enzymes = torch.where(~train_mask)[0]

    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']

        test_src = test_edges[0]
        test_dst = test_edges[1]
        unique_mets = torch.unique(test_dst)

        hits_10, hits_20, hits_50 = 0, 0, 0
        reciprocal_ranks = []
        count = 0

        for met_idx in unique_mets:
            met_mask = (test_dst == met_idx)
            true_enzs = test_src[met_mask]
            candidate_enzs = test_enzymes

            num_cands = candidate_enzs.size(0)
            eval_src = candidate_enzs
            eval_dst = met_idx.repeat(num_cands)
            eval_edges = torch.stack([eval_src, eval_dst])

            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()
            is_true = torch.isin(candidate_enzs, true_enzs)

            if is_true.sum() == 0:
                continue

            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]

            # Hits@K
            if sorted_labels[:10].sum() > 0:
                hits_10 += 1
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1
            if sorted_labels[:50].sum() > 0:
                hits_50 += 1

            # MRR
            true_positions = torch.where(sorted_labels)[0]
            if len(true_positions) > 0:
                first_true_rank = true_positions[0].item() + 1
                reciprocal_ranks.append(1.0 / first_true_rank)

            count += 1

    return {
        'hits10': hits_10 / count if count > 0 else 0.0,
        'hits20': hits_20 / count if count > 0 else 0.0,
        'hits50': hits_50 / count if count > 0 else 0.0,
        'mrr': np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0,
        'count': count
    }


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)

    print(f"Graph: {data.node_types}")
    print(f"Edge types: {len(data.edge_types)}")
    print()

    num_enzymes = data['Enzyme'].num_nodes
    seeds = [42, 123, 456]
    n_folds = 5

    # Configurations to test
    configs = [
        {'name': '2L-baseline', 'layers': 2, 'metapath': False},
        {'name': '2L-metapath', 'layers': 2, 'metapath': True},
        {'name': '3L-baseline', 'layers': 3, 'metapath': False},
        {'name': '3L-metapath', 'layers': 3, 'metapath': True},
    ]

    all_results = {cfg['name']: {'hits10': [], 'hits20': [], 'hits50': [], 'mrr': []}
                   for cfg in configs}

    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    enzyme_indices = np.arange(num_enzymes)

    print("=" * 80)
    print("METAPATH ABLATION STUDY")
    print("=" * 80)

    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(enzyme_indices)):
        train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
        train_mask[train_idx] = True
        test_mask = ~train_mask

        print(f"\n{'='*40}")
        print(f"Fold {fold_idx + 1}/{n_folds}")
        print(f"{'='*40}")

        for seed in seeds:
            for cfg in configs:
                print(f"\n  Config: {cfg['name']}, Seed: {seed}")

                # Reload data for each config (metapath modifies it)
                data_copy = torch.load('data/processed/strict_bipartite_v2.pt')
                data_copy = data_copy.to(device)

                result = train_and_evaluate(
                    data_copy, train_mask, test_mask, seed, device,
                    num_layers=cfg['layers'],
                    epochs=30,
                    use_metapath=cfg['metapath']
                )

                all_results[cfg['name']]['hits10'].append(result['hits10'])
                all_results[cfg['name']]['hits20'].append(result['hits20'])
                all_results[cfg['name']]['hits50'].append(result['hits50'])
                all_results[cfg['name']]['mrr'].append(result['mrr'])

                print(f"    H@10={result['hits10']:.3f}, H@20={result['hits20']:.3f}, "
                      f"MRR={result['mrr']:.4f}")

    # Summary
    print("\n" + "=" * 80)
    print("METAPATH ABLATION RESULTS")
    print("=" * 80)

    print(f"\n{'Config':<15} {'Hits@10':<18} {'Hits@20':<18} {'MRR':<18}")
    print("-" * 70)

    for cfg in configs:
        name = cfg['name']
        h10 = np.array(all_results[name]['hits10'])
        h20 = np.array(all_results[name]['hits20'])
        mrr = np.array(all_results[name]['mrr'])

        print(f"{name:<15} "
              f"{h10.mean():.4f} ± {h10.std():.4f}   "
              f"{h20.mean():.4f} ± {h20.std():.4f}   "
              f"{mrr.mean():.4f} ± {mrr.std():.4f}")

    # Comparison analysis
    print("\n" + "=" * 80)
    print("METAPATH EFFECT ANALYSIS")
    print("=" * 80)

    for base_layers in [2, 3]:
        base_name = f'{base_layers}L-baseline'
        meta_name = f'{base_layers}L-metapath'

        base_h20 = np.mean(all_results[base_name]['hits20'])
        meta_h20 = np.mean(all_results[meta_name]['hits20'])
        base_mrr = np.mean(all_results[base_name]['mrr'])
        meta_mrr = np.mean(all_results[meta_name]['mrr'])

        h20_delta = (meta_h20 - base_h20) / base_h20 * 100 if base_h20 > 0 else 0
        mrr_delta = (meta_mrr - base_mrr) / base_mrr * 100 if base_mrr > 0 else 0

        print(f"{base_layers}-layer: metapath vs baseline -> "
              f"Hits@20 {h20_delta:+.1f}%, MRR {mrr_delta:+.1f}%")

    # Save results
    os.makedirs('results/gnn', exist_ok=True)

    with open('results/gnn/metapath_ablation.tsv', 'w') as f:
        f.write("Config\tLayers\tMetapath\tHits10_Mean\tHits10_Std\t"
                "Hits20_Mean\tHits20_Std\tMRR_Mean\tMRR_Std\n")
        for cfg in configs:
            name = cfg['name']
            h10 = np.array(all_results[name]['hits10'])
            h20 = np.array(all_results[name]['hits20'])
            mrr = np.array(all_results[name]['mrr'])
            f.write(f"{name}\t{cfg['layers']}\t{cfg['metapath']}\t"
                    f"{h10.mean():.4f}\t{h10.std():.4f}\t"
                    f"{h20.mean():.4f}\t{h20.std():.4f}\t"
                    f"{mrr.mean():.4f}\t{mrr.std():.4f}\n")

    print("\nResults saved to: results/gnn/metapath_ablation.tsv")


if __name__ == "__main__":
    main()
