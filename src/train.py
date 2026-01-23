"""
Unified Training Pipeline for Ethylene-Isoflavonoid GNN

This is the RECOMMENDED training script that integrates all
standardized components for reproducible experiments.

Features:
- Multi-seed evaluation support
- Configurable via centralized config
- Standardized negative sampling
- Proper train/val/test splits with early stopping
- Comprehensive metrics logging

Usage:
    # Single run
    python src/train.py --graph data/processed/strict_bipartite_v2.pt

    # Multi-seed evaluation
    python src/train.py --seeds 42,123,456,789,1024

    # Custom configuration
    python src/train.py --num-layers 3 --epochs 100 --patience 15
"""

import torch
import torch.nn as nn
import numpy as np
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import (
    DATA_PATHS, TRAINING_CONFIG, EVALUATION_CONFIG, SPLIT_CONFIG
)
from src.utils.seed import set_seed, get_device
from src.utils.negative_sampling import NegativeSampler
from src.utils.data_split import DataSplitter
from src.model import HGT, LinkPredictor


def train_epoch(
    model: nn.Module,
    predictor: nn.Module,
    optimizer: torch.optim.Optimizer,
    data,
    train_edges: torch.Tensor,
    negative_sampler: NegativeSampler,
    device: torch.device,
) -> float:
    """
    Execute a single training epoch.

    Args:
        model: GNN encoder
        predictor: Link predictor
        optimizer: Optimizer
        data: Training graph data
        train_edges: Training edge indices
        negative_sampler: Negative sampling strategy
        device: Compute device

    Returns:
        Training loss for this epoch
    """
    model.train()
    optimizer.zero_grad()

    # Forward pass
    x_dict = model(data.x_dict, data.edge_index_dict)

    # Positive edges
    pos_edge_index = train_edges

    # Negative sampling using standardized sampler
    num_metabolites = data['Metabolite'].num_nodes
    neg_src, neg_dst = negative_sampler.sample(
        pos_edge_index, num_metabolites, device
    )

    # Predictions
    pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
    neg_out = predictor(
        x_dict['Enzyme'], x_dict['Metabolite'],
        torch.stack([neg_src, neg_dst])
    )

    # BCE loss
    loss = (
        -torch.log(torch.sigmoid(pos_out) + 1e-15).mean()
        - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
    )

    loss.backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), TRAINING_CONFIG['max_grad_norm'])
    optimizer.step()

    return loss.item()


def evaluate(
    model: nn.Module,
    predictor: nn.Module,
    data,
    edge_index: torch.Tensor,
    candidate_nodes: torch.Tensor,
    device: torch.device,
    hits_at_k: List[int] = None,
) -> Dict[str, float]:
    """
    Evaluate model on a split.

    Args:
        model: GNN encoder
        predictor: Link predictor
        data: Graph data
        edge_index: Edges to evaluate
        candidate_nodes: Candidate nodes for ranking
        device: Compute device
        hits_at_k: List of K values for Hits@K metrics

    Returns:
        Dictionary of metrics
    """
    if hits_at_k is None:
        hits_at_k = EVALUATION_CONFIG['hits_at_k']

    model.eval()

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']

        src = edge_index[0]
        dst = edge_index[1]
        unique_mets = torch.unique(dst)

        # Initialize counters
        hits = {k: 0 for k in hits_at_k}
        reciprocal_ranks = []
        count = 0

        for met_idx in unique_mets:
            met_mask = (dst == met_idx)
            true_enzs = src[met_mask]

            # Score all candidates
            num_cands = candidate_nodes.size(0)
            eval_src = candidate_nodes
            eval_dst = met_idx.repeat(num_cands)
            eval_edges = torch.stack([eval_src, eval_dst])

            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()
            is_true = torch.isin(candidate_nodes, true_enzs)

            if is_true.sum() == 0:
                continue

            # Ranking
            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]

            # Hits@K
            for k in hits_at_k:
                if sorted_labels[:k].sum() > 0:
                    hits[k] += 1

            # MRR
            true_positions = torch.where(sorted_labels)[0]
            if len(true_positions) > 0:
                first_true_rank = true_positions[0].item() + 1
                reciprocal_ranks.append(1.0 / first_true_rank)

            count += 1

    # Compute metrics
    metrics = {f'hits@{k}': hits[k] / count if count > 0 else 0.0 for k in hits_at_k}
    metrics['mrr'] = np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0
    metrics['count'] = count

    return metrics


def run_experiment(
    graph_path: Path,
    num_layers: int,
    seed: int,
    device: torch.device,
    epochs: int = None,
    patience: int = None,
    neg_strategy: str = 'hard',
) -> Dict[str, Dict[str, float]]:
    """
    Run a single training experiment.

    Args:
        graph_path: Path to graph file
        num_layers: Number of GNN layers
        seed: Random seed
        device: Compute device
        epochs: Max training epochs
        patience: Early stopping patience
        neg_strategy: Negative sampling strategy

    Returns:
        Dictionary with train/val/test metrics
    """
    set_seed(seed)

    epochs = epochs or TRAINING_CONFIG['epochs']
    patience = patience or TRAINING_CONFIG['patience']

    # Load data
    data = torch.load(graph_path)
    data = data.to(device)

    # Split data
    splitter = DataSplitter(seed=seed)
    train_data, train_edges, val_edges, test_edges = splitter.split(data)

    # Get node masks
    train_mask, val_mask, test_mask = splitter.get_node_masks(train_data)
    train_nodes = torch.where(train_mask)[0]
    val_nodes = torch.where(val_mask)[0]
    test_nodes = torch.where(test_mask)[0]

    # Initialize model
    model = HGT(
        train_data.metadata(),
        in_channels=TRAINING_CONFIG['hidden_channels'],
        hidden_channels=TRAINING_CONFIG['hidden_channels'],
        out_channels=TRAINING_CONFIG['out_channels'],
        num_heads=TRAINING_CONFIG['num_heads'],
        num_layers=num_layers,
    ).to(device)

    predictor = LinkPredictor(TRAINING_CONFIG['out_channels']).to(device)

    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=TRAINING_CONFIG['learning_rate'],
        weight_decay=TRAINING_CONFIG['weight_decay'],
    )

    # Initialize negative sampler
    # For EC-class strategy, load EC mappings from the graph
    ec_to_indices = getattr(data['Metabolite'], 'ec_to_indices', None)
    met_to_ecs = getattr(data['Metabolite'], 'met_to_ecs', None)

    if neg_strategy == 'ec_class' and ec_to_indices and met_to_ecs:
        negative_sampler = NegativeSampler(
            strategy=neg_strategy,
            ec_to_indices=ec_to_indices,
            met_to_ecs=met_to_ecs,
        )
    else:
        negative_sampler = NegativeSampler(strategy=neg_strategy)

    # Training with early stopping
    best_val_mrr = 0
    best_epoch = 0
    patience_counter = 0
    best_state = None

    for epoch in range(1, epochs + 1):
        loss = train_epoch(
            model, predictor, optimizer, train_data,
            train_edges, negative_sampler, device
        )

        # Validate every 5 epochs
        if epoch % 5 == 0:
            val_metrics = evaluate(
                model, predictor, train_data, val_edges, val_nodes, device
            )

            if val_metrics['mrr'] > best_val_mrr:
                best_val_mrr = val_metrics['mrr']
                best_epoch = epoch
                patience_counter = 0
                best_state = {
                    'model': model.state_dict(),
                    'predictor': predictor.state_dict()
                }
            else:
                patience_counter += 1

            if patience_counter >= patience // 5:
                break

    # Load best model
    if best_state is not None:
        model.load_state_dict(best_state['model'])
        predictor.load_state_dict(best_state['predictor'])

    # Final evaluation
    train_metrics = evaluate(model, predictor, train_data, train_edges, train_nodes, device)
    val_metrics = evaluate(model, predictor, train_data, val_edges, val_nodes, device)
    test_metrics = evaluate(model, predictor, train_data, test_edges, test_nodes, device)

    return {
        'train': train_metrics,
        'val': val_metrics,
        'test': test_metrics,
        'best_epoch': best_epoch,
    }


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Unified training pipeline for Ethylene GNN"
    )
    parser.add_argument(
        '--graph', type=str,
        default=str(DATA_PATHS['bipartite_graph_v2']),
        help='Path to graph file'
    )
    parser.add_argument(
        '--seeds', type=str, default='42,123,456,789,1024',
        help='Comma-separated random seeds'
    )
    parser.add_argument(
        '--num-layers', type=int, default=2,
        help='Number of GNN layers'
    )
    parser.add_argument(
        '--epochs', type=int, default=TRAINING_CONFIG['epochs'],
        help='Max training epochs'
    )
    parser.add_argument(
        '--patience', type=int, default=TRAINING_CONFIG['patience'],
        help='Early stopping patience'
    )
    parser.add_argument(
        '--neg-strategy', type=str, default='ec_class',
        choices=['random', 'hard', 'mixed', 'ec_class'],
        help='Negative sampling strategy (ec_class recommended)'
    )
    parser.add_argument(
        '--output', type=str, default=None,
        help='Output file for results (default: results/gnn/train_results.tsv)'
    )

    args = parser.parse_args()

    # Parse seeds
    seeds = [int(s.strip()) for s in args.seeds.split(',')]

    device = get_device()
    print(f"Device: {device}")
    print(f"Graph: {args.graph}")
    print(f"Seeds: {seeds}")
    print(f"Layers: {args.num_layers}")
    print(f"Negative sampling: {args.neg_strategy}")
    print("=" * 60)

    # Load data info
    data = torch.load(args.graph)
    print(f"Enzymes: {data['Enzyme'].num_nodes}")
    print(f"Metabolites: {data['Metabolite'].num_nodes}")
    print(f"Edges: {data['Enzyme', 'catalyzes', 'Metabolite'].edge_index.shape[1]}")

    # Check EC mapping availability
    ec_to_indices = getattr(data['Metabolite'], 'ec_to_indices', None)
    met_to_ecs = getattr(data['Metabolite'], 'met_to_ecs', None)
    if ec_to_indices and met_to_ecs:
        print(f"EC mappings: {len(ec_to_indices)} EC classes, "
              f"{len(met_to_ecs)} metabolites with EC info")
    elif args.neg_strategy == 'ec_class':
        print("Warning: EC mappings not found in graph. "
              "Run data_pipeline.py to rebuild with EC mappings.")
    print()

    # Run experiments
    all_results = {'train': [], 'val': [], 'test': []}

    for seed in seeds:
        print(f"Seed {seed}...", end=" ", flush=True)

        result = run_experiment(
            Path(args.graph),
            args.num_layers,
            seed,
            device,
            args.epochs,
            args.patience,
            args.neg_strategy,
        )

        all_results['train'].append(result['train'])
        all_results['val'].append(result['val'])
        all_results['test'].append(result['test'])

        print(
            f"Test H@20={result['test']['hits@20']:.3f}, "
            f"MRR={result['test']['mrr']:.4f} (epoch {result['best_epoch']})"
        )

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY (Test Set)")
    print("=" * 60)

    test_results = all_results['test']
    for metric in ['hits@1', 'hits@10', 'hits@20', 'mrr']:
        values = np.array([r[metric] for r in test_results])
        if 'hits' in metric:
            print(f"{metric}: {values.mean()*100:.2f}% +/- {values.std()*100:.2f}%")
        else:
            print(f"{metric}: {values.mean():.4f} +/- {values.std():.4f}")

    # Save results
    output_path = args.output or str(DATA_PATHS['gnn_results'] / 'train_results.tsv')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, 'w') as f:
        f.write("seed\tsplit\thits@1\thits@3\thits@10\thits@20\thits@50\tmrr\n")
        for i, seed in enumerate(seeds):
            for split in ['train', 'val', 'test']:
                r = all_results[split][i]
                f.write(f"{seed}\t{split}\t")
                f.write(f"{r['hits@1']:.4f}\t{r['hits@3']:.4f}\t")
                f.write(f"{r['hits@10']:.4f}\t{r['hits@20']:.4f}\t")
                f.write(f"{r['hits@50']:.4f}\t{r['mrr']:.4f}\n")

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
