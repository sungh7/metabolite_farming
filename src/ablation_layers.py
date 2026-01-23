"""
GNN Layer Depth Ablation Study

Compares HGT performance with different numbers of layers:
- 2 layers (baseline): TF can reach Enzyme directly
- 3 layers: TF can reach Metabolite (2-hop)
- 4 layers: Deeper propagation

Hypothesis: 3+ layers should improve TF→Metabolite signal propagation.
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


def train_and_evaluate(data, train_mask, test_mask, seed, device,
                       num_layers=2, epochs=30):
    """
    Train and evaluate HGT with specified number of layers.
    """
    set_seed(seed)

    train_data = data.clone()
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

    # Model with specified layers
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
    print(f"Enzymes: {data['Enzyme'].num_nodes}")
    print(f"Metabolites: {data['Metabolite'].num_nodes}")
    print(f"Enzyme→Metabolite edges: {data['Enzyme', 'catalyzes', 'Metabolite'].edge_index.shape[1]}")
    print()

    num_enzymes = data['Enzyme'].num_nodes
    seeds = [42, 123, 456]
    n_folds = 5

    # Test different layer depths
    layer_configs = [2, 3, 4]
    all_results = {n: {'hits10': [], 'hits20': [], 'hits50': [], 'mrr': []}
                   for n in layer_configs}

    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    enzyme_indices = np.arange(num_enzymes)

    print("=" * 70)
    print("LAYER DEPTH ABLATION STUDY")
    print("=" * 70)

    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(enzyme_indices)):
        train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
        train_mask[train_idx] = True
        test_mask = ~train_mask

        print(f"\n--- Fold {fold_idx + 1}/{n_folds} ---")

        for seed in seeds:
            for num_layers in layer_configs:
                result = train_and_evaluate(
                    data, train_mask, test_mask, seed, device,
                    num_layers=num_layers, epochs=30
                )

                all_results[num_layers]['hits10'].append(result['hits10'])
                all_results[num_layers]['hits20'].append(result['hits20'])
                all_results[num_layers]['hits50'].append(result['hits50'])
                all_results[num_layers]['mrr'].append(result['mrr'])

                print(f"  Layers={num_layers}, Seed={seed}: "
                      f"H@10={result['hits10']:.3f}, H@20={result['hits20']:.3f}, "
                      f"MRR={result['mrr']:.4f}")

    # Summary
    print("\n" + "=" * 70)
    print("LAYER DEPTH ABLATION RESULTS")
    print("=" * 70)

    print(f"\n{'Layers':<10} {'Hits@10':<20} {'Hits@20':<20} {'Hits@50':<20} {'MRR':<20}")
    print("-" * 90)

    for num_layers in layer_configs:
        h10 = np.array(all_results[num_layers]['hits10'])
        h20 = np.array(all_results[num_layers]['hits20'])
        h50 = np.array(all_results[num_layers]['hits50'])
        mrr = np.array(all_results[num_layers]['mrr'])

        print(f"{num_layers:<10} "
              f"{h10.mean():.4f} ± {h10.std():.4f}    "
              f"{h20.mean():.4f} ± {h20.std():.4f}    "
              f"{h50.mean():.4f} ± {h50.std():.4f}    "
              f"{mrr.mean():.4f} ± {mrr.std():.4f}")

    # Improvement analysis
    print("\n" + "=" * 70)
    print("IMPROVEMENT ANALYSIS (vs 2-layer baseline)")
    print("=" * 70)

    baseline_h20 = np.mean(all_results[2]['hits20'])
    baseline_mrr = np.mean(all_results[2]['mrr'])

    for num_layers in [3, 4]:
        h20_mean = np.mean(all_results[num_layers]['hits20'])
        mrr_mean = np.mean(all_results[num_layers]['mrr'])

        h20_delta = (h20_mean - baseline_h20) / baseline_h20 * 100
        mrr_delta = (mrr_mean - baseline_mrr) / baseline_mrr * 100

        print(f"{num_layers}-layer vs 2-layer: "
              f"Hits@20 {h20_delta:+.1f}%, MRR {mrr_delta:+.1f}%")

    # Save results
    os.makedirs('results/gnn', exist_ok=True)

    with open('results/gnn/layer_ablation.tsv', 'w') as f:
        f.write("Layers\tHits10_Mean\tHits10_Std\tHits20_Mean\tHits20_Std\t"
                "Hits50_Mean\tHits50_Std\tMRR_Mean\tMRR_Std\n")
        for num_layers in layer_configs:
            h10 = np.array(all_results[num_layers]['hits10'])
            h20 = np.array(all_results[num_layers]['hits20'])
            h50 = np.array(all_results[num_layers]['hits50'])
            mrr = np.array(all_results[num_layers]['mrr'])
            f.write(f"{num_layers}\t{h10.mean():.4f}\t{h10.std():.4f}\t"
                    f"{h20.mean():.4f}\t{h20.std():.4f}\t"
                    f"{h50.mean():.4f}\t{h50.std():.4f}\t"
                    f"{mrr.mean():.4f}\t{mrr.std():.4f}\n")

    print("\nResults saved to: results/gnn/layer_ablation.tsv")


if __name__ == "__main__":
    main()
