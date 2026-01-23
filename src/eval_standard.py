"""
Standard Evaluation with Train/Val/Test Split

Produces performance tables in two formats:
1. Modern (Current Academic): Single test set, report mean±std across seeds
2. Traditional: Train/Val/Test metrics separately

Following conventions from:
- OGB (Open Graph Benchmark)
- PyG examples
- Knowledge Graph papers (TransE, RotatE, etc.)
"""

import torch
import numpy as np
from sklearn.model_selection import train_test_split
import os
import sys
sys.path.append(os.getcwd())

from src.model import HGT, LinkPredictor


def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def create_splits(num_nodes, train_ratio=0.7, val_ratio=0.15, test_ratio=0.15, seed=42):
    """Create train/val/test splits."""
    indices = np.arange(num_nodes)

    # First split: train vs (val+test)
    train_idx, temp_idx = train_test_split(
        indices, train_size=train_ratio, random_state=seed
    )

    # Second split: val vs test
    val_ratio_adjusted = val_ratio / (val_ratio + test_ratio)
    val_idx, test_idx = train_test_split(
        temp_idx, train_size=val_ratio_adjusted, random_state=seed
    )

    return train_idx, val_idx, test_idx


def evaluate_split(model, predictor, data, edge_index, candidate_nodes, device):
    """Evaluate on a specific split."""
    model.eval()

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']

        src = edge_index[0]
        dst = edge_index[1]
        unique_mets = torch.unique(dst)

        hits_1, hits_3, hits_10, hits_20, hits_50 = 0, 0, 0, 0, 0
        reciprocal_ranks = []
        count = 0

        for met_idx in unique_mets:
            met_mask = (dst == met_idx)
            true_enzs = src[met_mask]

            num_cands = candidate_nodes.size(0)
            eval_src = candidate_nodes
            eval_dst = met_idx.repeat(num_cands)
            eval_edges = torch.stack([eval_src, eval_dst])

            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()
            is_true = torch.isin(candidate_nodes, true_enzs)

            if is_true.sum() == 0:
                continue

            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]

            # Hits@K
            if sorted_labels[:1].sum() > 0:
                hits_1 += 1
            if sorted_labels[:3].sum() > 0:
                hits_3 += 1
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
        'hits@1': hits_1 / count if count > 0 else 0.0,
        'hits@3': hits_3 / count if count > 0 else 0.0,
        'hits@10': hits_10 / count if count > 0 else 0.0,
        'hits@20': hits_20 / count if count > 0 else 0.0,
        'hits@50': hits_50 / count if count > 0 else 0.0,
        'mrr': np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0,
        'count': count
    }


def train_epoch(model, predictor, optimizer, data, train_edges, device):
    """Single training epoch."""
    model.train()
    optimizer.zero_grad()

    x_dict = model(data.x_dict, data.edge_index_dict)

    # Hard Negative Sampling
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


def run_experiment(data, num_layers, seed, device, epochs=50, patience=10):
    """Run single experiment with early stopping."""
    set_seed(seed)

    num_enzymes = data['Enzyme'].num_nodes

    # Create node splits
    train_idx, val_idx, test_idx = create_splits(num_enzymes, seed=seed)

    train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    val_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    test_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)

    train_mask[train_idx] = True
    val_mask[val_idx] = True
    test_mask[test_idx] = True

    # Split edges based on source enzyme
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]

    train_edge_mask = train_mask[src]
    val_edge_mask = val_mask[src]
    test_edge_mask = test_mask[src]

    train_edges = edge_index[:, train_edge_mask]
    val_edges = edge_index[:, val_edge_mask]
    test_edges = edge_index[:, test_edge_mask]

    # Create training graph (only train edges)
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges

    # Handle reverse edges
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_mask = train_mask[rev_index[1]]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]

    # Model
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
        lr=0.01, weight_decay=1e-5
    )

    # Training with early stopping
    best_val_mrr = 0
    best_epoch = 0
    patience_counter = 0

    train_nodes = torch.where(train_mask)[0]
    val_nodes = torch.where(val_mask)[0]
    test_nodes = torch.where(test_mask)[0]

    for epoch in range(1, epochs + 1):
        loss = train_epoch(model, predictor, optimizer, train_data, train_edges, device)

        # Validate every 5 epochs
        if epoch % 5 == 0:
            val_metrics = evaluate_split(model, predictor, train_data, val_edges, val_nodes, device)

            if val_metrics['mrr'] > best_val_mrr:
                best_val_mrr = val_metrics['mrr']
                best_epoch = epoch
                patience_counter = 0
                # Save best model state
                best_state = {
                    'model': model.state_dict(),
                    'predictor': predictor.state_dict()
                }
            else:
                patience_counter += 1

            if patience_counter >= patience // 5:
                break

    # Load best model
    model.load_state_dict(best_state['model'])
    predictor.load_state_dict(best_state['predictor'])

    # Final evaluation on all splits
    train_metrics = evaluate_split(model, predictor, train_data, train_edges, train_nodes, device)
    val_metrics = evaluate_split(model, predictor, train_data, val_edges, val_nodes, device)
    test_metrics = evaluate_split(model, predictor, train_data, test_edges, test_nodes, device)

    return {
        'train': train_metrics,
        'val': val_metrics,
        'test': test_metrics,
        'best_epoch': best_epoch
    }


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)

    print(f"Enzymes: {data['Enzyme'].num_nodes}")
    print(f"Metabolites: {data['Metabolite'].num_nodes}")
    print(f"Edges: {data['Enzyme', 'catalyzes', 'Metabolite'].edge_index.shape[1]}")
    print()

    seeds = [42, 123, 456, 789, 1024]
    layer_configs = [2, 3]

    all_results = {}

    print("=" * 80)
    print("STANDARD EVALUATION (Train/Val/Test Split)")
    print("=" * 80)

    for num_layers in layer_configs:
        config_name = f'{num_layers}L-HGT'
        all_results[config_name] = {'train': [], 'val': [], 'test': []}

        print(f"\n--- {config_name} ---")

        for seed in seeds:
            print(f"  Seed {seed}...", end=" ", flush=True)

            result = run_experiment(data.clone(), num_layers, seed, device)

            all_results[config_name]['train'].append(result['train'])
            all_results[config_name]['val'].append(result['val'])
            all_results[config_name]['test'].append(result['test'])

            print(f"Test H@20={result['test']['hits@20']:.3f}, "
                  f"MRR={result['test']['mrr']:.4f} (epoch {result['best_epoch']})")

    # =========================================================================
    # Format 1: Modern Academic Style (OGB-like)
    # =========================================================================
    print("\n" + "=" * 80)
    print("FORMAT 1: Modern Academic Style (Test Set Only)")
    print("Following OGB/PyG conventions")
    print("=" * 80)

    print(f"\n{'Model':<12} {'Hits@1':<14} {'Hits@10':<14} {'Hits@20':<14} {'MRR':<14}")
    print("-" * 70)

    for config_name in all_results:
        test_results = all_results[config_name]['test']

        h1 = np.array([r['hits@1'] for r in test_results])
        h10 = np.array([r['hits@10'] for r in test_results])
        h20 = np.array([r['hits@20'] for r in test_results])
        mrr = np.array([r['mrr'] for r in test_results])

        print(f"{config_name:<12} "
              f"{h1.mean()*100:.2f}±{h1.std()*100:.2f}%    "
              f"{h10.mean()*100:.2f}±{h10.std()*100:.2f}%    "
              f"{h20.mean()*100:.2f}±{h20.std()*100:.2f}%    "
              f"{mrr.mean():.4f}±{mrr.std():.4f}")

    # =========================================================================
    # Format 2: Traditional Style (Train/Val/Test)
    # =========================================================================
    print("\n" + "=" * 80)
    print("FORMAT 2: Traditional Style (Train/Val/Test)")
    print("Following classic ML conventions")
    print("=" * 80)

    for config_name in all_results:
        print(f"\n### {config_name} ###")
        print(f"{'Split':<8} {'Hits@1':<12} {'Hits@10':<12} {'Hits@20':<12} {'MRR':<12}")
        print("-" * 60)

        for split in ['train', 'val', 'test']:
            results = all_results[config_name][split]

            h1 = np.array([r['hits@1'] for r in results])
            h10 = np.array([r['hits@10'] for r in results])
            h20 = np.array([r['hits@20'] for r in results])
            mrr = np.array([r['mrr'] for r in results])

            print(f"{split:<8} "
                  f"{h1.mean()*100:>5.2f}%       "
                  f"{h10.mean()*100:>5.2f}%       "
                  f"{h20.mean()*100:>5.2f}%       "
                  f"{mrr.mean():.4f}")

    # =========================================================================
    # Format 3: Paper-Ready Table (LaTeX-friendly)
    # =========================================================================
    print("\n" + "=" * 80)
    print("FORMAT 3: Paper-Ready Table")
    print("=" * 80)

    print("\n| Model | Split | Hits@1 | Hits@10 | Hits@20 | MRR |")
    print("|-------|-------|--------|---------|---------|-----|")

    for config_name in all_results:
        for split in ['train', 'val', 'test']:
            results = all_results[config_name][split]

            h1 = np.mean([r['hits@1'] for r in results]) * 100
            h10 = np.mean([r['hits@10'] for r in results]) * 100
            h20 = np.mean([r['hits@20'] for r in results]) * 100
            mrr = np.mean([r['mrr'] for r in results])

            print(f"| {config_name} | {split} | {h1:.2f}% | {h10:.2f}% | {h20:.2f}% | {mrr:.4f} |")

    # =========================================================================
    # Save Results
    # =========================================================================
    os.makedirs('results/gnn', exist_ok=True)

    # TSV format
    with open('results/gnn/standard_eval.tsv', 'w') as f:
        f.write("Model\tSplit\tHits@1\tHits@3\tHits@10\tHits@20\tHits@50\tMRR\tStd_H20\tStd_MRR\n")

        for config_name in all_results:
            for split in ['train', 'val', 'test']:
                results = all_results[config_name][split]

                h1 = np.array([r['hits@1'] for r in results])
                h3 = np.array([r['hits@3'] for r in results])
                h10 = np.array([r['hits@10'] for r in results])
                h20 = np.array([r['hits@20'] for r in results])
                h50 = np.array([r['hits@50'] for r in results])
                mrr = np.array([r['mrr'] for r in results])

                f.write(f"{config_name}\t{split}\t"
                        f"{h1.mean():.4f}\t{h3.mean():.4f}\t{h10.mean():.4f}\t"
                        f"{h20.mean():.4f}\t{h50.mean():.4f}\t{mrr.mean():.4f}\t"
                        f"{h20.std():.4f}\t{mrr.std():.4f}\n")

    print("\nResults saved to: results/gnn/standard_eval.tsv")

    # =========================================================================
    # Improvement Summary
    # =========================================================================
    print("\n" + "=" * 80)
    print("IMPROVEMENT SUMMARY (3L vs 2L)")
    print("=" * 80)

    for split in ['val', 'test']:
        h20_2l = np.mean([r['hits@20'] for r in all_results['2L-HGT'][split]])
        h20_3l = np.mean([r['hits@20'] for r in all_results['3L-HGT'][split]])
        mrr_2l = np.mean([r['mrr'] for r in all_results['2L-HGT'][split]])
        mrr_3l = np.mean([r['mrr'] for r in all_results['3L-HGT'][split]])

        h20_delta = (h20_3l - h20_2l) / h20_2l * 100
        mrr_delta = (mrr_3l - mrr_2l) / mrr_2l * 100

        print(f"{split.upper():>5}: Hits@20 {h20_delta:+.1f}%, MRR {mrr_delta:+.1f}%")


if __name__ == "__main__":
    main()
