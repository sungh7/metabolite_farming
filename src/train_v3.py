"""
Training Pipeline v3 for Ethylene-Isoflavonoid GNN Project

학술적으로 방어 가능한 평가 (Academically Defensible Evaluation)

Key Changes from v2:
1. Filtered ranking evaluation (filter known positives)
2. Global→Local index mapping for candidate subsets
3. Ablation support (--no-rxn-neighbor, --no-tf-edges)
4. Two-tier evaluation: Full (306) and Experimental (10)

Usage:
    # Main training
    python src/train_v3.py --graph data/processed/graph_v3.pt --seeds 42,123,456

    # Ablation (no rxn_neighbor)
    python src/train_v3.py --graph data/processed/graph_v3.pt --no-rxn-neighbor
"""

import torch
import torch.nn as nn
import numpy as np
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from datetime import datetime
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import DATA_PATHS, TRAINING_CONFIG, EVALUATION_CONFIG, V3_CONFIG
from src.utils.seed import set_seed, get_device
from src.utils.negative_sampling import NegativeSampler
from src.utils.data_split import DataSplitter
from src.model_v3 import HGTv3, HGTv3Ablation, LinkPredictor, create_model_v3


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

    Supervision uses Tier-R edges only (catalyzes_R).

    Args:
        model: HGT v3 encoder
        predictor: Link predictor
        optimizer: Optimizer
        data: Training graph data
        train_edges: Training edge indices (Tier-R only)
        negative_sampler: Negative sampling strategy
        device: Compute device

    Returns:
        Training loss for this epoch
    """
    model.train()
    optimizer.zero_grad()

    # Forward pass
    x_dict = model(data)

    # Positive edges (Tier-R only)
    pos_edge_index = train_edges

    # Negative sampling
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


def evaluate_filtered(
    model: nn.Module,
    predictor: nn.Module,
    data,
    test_edges: torch.Tensor,
    tier_r_lookup: Dict[int, Set[int]],
    candidate_mets: List[int],
    device: torch.device,
    hits_at_k: List[int] = None,
) -> Dict[str, float]:
    """
    Filtered ranking evaluation.

    For each test edge (enzyme, metabolite):
    1. Score all candidate metabolites
    2. Filter (set to -inf) known Tier-R positives for that enzyme
    3. Compute rank of the true metabolite

    This prevents pessimistic rank estimates due to known positives
    ranking higher than the test positive.

    Args:
        model: HGT v3 encoder
        predictor: Link predictor
        data: Graph data
        test_edges: Test edges [2, num_test]
        tier_r_lookup: Dict[enzyme_idx, Set[metabolite_idx]] of known positives
        candidate_mets: List of candidate metabolite indices
        device: Compute device
        hits_at_k: List of K values for Hits@K

    Returns:
        Dictionary of metrics
    """
    if hits_at_k is None:
        hits_at_k = EVALUATION_CONFIG['hits_at_k']

    model.eval()

    with torch.no_grad():
        x_dict = model(data)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']

        # Global→Local mapping for candidate subset
        cand = list(candidate_mets)
        global_to_local = {m: i for i, m in enumerate(cand)}
        cand_tensor = torch.tensor(cand, device=device)

        # Get candidate metabolite embeddings
        cand_met_emb = met_emb[cand_tensor]  # [num_cand, dim]

        # Initialize metrics
        ranks = []
        hits = {k: 0 for k in hits_at_k}

        # Process each test edge
        src = test_edges[0].tolist()
        dst = test_edges[1].tolist()

        for enz_idx, met_idx in zip(src, dst):
            # Skip if metabolite not in candidate set
            if met_idx not in global_to_local:
                continue

            # Score all candidates: enz_emb[enz_idx] @ cand_met_emb.T
            enz_vec = enz_emb[enz_idx].unsqueeze(0)  # [1, dim]
            scores = torch.mm(enz_vec, cand_met_emb.t()).squeeze(0)  # [num_cand]

            # Filter known Tier-R positives (set to -inf)
            known_positives = tier_r_lookup.get(enz_idx, set())
            for known_met in known_positives:
                if known_met != met_idx and known_met in global_to_local:
                    local_idx = global_to_local[known_met]
                    scores[local_idx] = float('-inf')

            # Compute rank (1-based)
            dst_local = global_to_local[met_idx]
            dst_score = scores[dst_local]
            rank = (scores > dst_score).sum().item() + 1
            ranks.append(rank)

            # Hits@K
            for k in hits_at_k:
                if rank <= k:
                    hits[k] += 1

    # Compute final metrics
    count = len(ranks)
    if count == 0:
        return {f'hits@{k}': 0.0 for k in hits_at_k} | {'mrr': 0.0, 'count': 0}

    metrics = {f'hits@{k}': hits[k] / count for k in hits_at_k}
    metrics['mrr'] = np.mean([1.0 / r for r in ranks])
    metrics['count'] = count
    metrics['mean_rank'] = np.mean(ranks)

    return metrics


def evaluate_full(
    model: nn.Module,
    predictor: nn.Module,
    data,
    test_edges: torch.Tensor,
    tier_r_lookup: Dict[int, Set[int]],
    device: torch.device,
) -> Dict[str, float]:
    """
    Full evaluation on all metabolites (306).

    Args:
        model, predictor, data, test_edges, tier_r_lookup, device: As above

    Returns:
        Metrics dictionary
    """
    num_mets = data['Metabolite'].num_nodes
    all_mets = list(range(num_mets))
    return evaluate_filtered(
        model, predictor, data, test_edges,
        tier_r_lookup, all_mets, device
    )


def evaluate_experimental(
    model: nn.Module,
    predictor: nn.Module,
    data,
    test_edges: torch.Tensor,
    tier_r_lookup: Dict[int, Set[int]],
    device: torch.device,
) -> Dict[str, float]:
    """
    Evaluation on experimental metabolites only.

    Args:
        model, predictor, data, test_edges, tier_r_lookup, device: As above

    Returns:
        Metrics dictionary
    """
    # Get experimental metabolite indices
    is_exp = data['Metabolite'].is_experimental
    exp_mets = torch.where(is_exp)[0].tolist()

    if len(exp_mets) == 0:
        return {'hits@1': 0, 'mrr': 0, 'count': 0}

    # Filter test edges to only include experimental metabolites
    exp_set = set(exp_mets)
    mask = torch.tensor([dst.item() in exp_set for dst in test_edges[1]])
    filtered_edges = test_edges[:, mask]

    if filtered_edges.size(1) == 0:
        return {'hits@1': 0, 'mrr': 0, 'count': 0}

    return evaluate_filtered(
        model, predictor, data, filtered_edges,
        tier_r_lookup, exp_mets, device
    )


def run_experiment(
    graph_path: Path,
    num_layers: int,
    seed: int,
    device: torch.device,
    epochs: int = None,
    patience: int = None,
    neg_strategy: str = 'ec_class',
    use_rxn_neighbor: bool = True,
    use_tf_edges: bool = True,
) -> Dict[str, Dict[str, float]]:
    """
    Run a single training experiment.

    Args:
        graph_path: Path to v3 graph file
        num_layers: Number of GNN layers
        seed: Random seed
        device: Compute device
        epochs: Max training epochs
        patience: Early stopping patience
        neg_strategy: Negative sampling strategy
        use_rxn_neighbor: Whether to use rxn_neighbor edges
        use_tf_edges: Whether to use TF edges

    Returns:
        Dictionary with train/val/test metrics (full and experimental)
    """
    set_seed(seed)

    epochs = epochs or V3_CONFIG['epochs']
    patience = patience or V3_CONFIG['patience']

    # Load data
    data = torch.load(graph_path)
    data = data.to(device)

    # Get tier_r_lookup
    tier_r_lookup = getattr(data, 'tier_r_lookup', {})
    if not tier_r_lookup:
        print("Warning: tier_r_lookup not found. Filtered evaluation may be inaccurate.")

    # Split data using Tier-R edges
    # Get catalyzes_R edge index
    if ('Enzyme', 'catalyzes_R', 'Metabolite') in data.edge_types:
        tier_r_edges = data['Enzyme', 'catalyzes_R', 'Metabolite'].edge_index
    else:
        # Fallback to catalyzes
        tier_r_edges = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index

    # Manual split for supervision edges
    num_edges = tier_r_edges.size(1)
    perm = torch.randperm(num_edges)

    train_size = int(0.7 * num_edges)
    val_size = int(0.15 * num_edges)

    train_idx = perm[:train_size]
    val_idx = perm[train_size:train_size + val_size]
    test_idx = perm[train_size + val_size:]

    train_edges = tier_r_edges[:, train_idx]
    val_edges = tier_r_edges[:, val_idx]
    test_edges = tier_r_edges[:, test_idx]

    # Create model
    model, predictor = create_model_v3(
        data, num_layers,
        use_rxn_neighbor=use_rxn_neighbor,
        use_tf_edges=use_tf_edges,
    )
    model = model.to(device)
    predictor = predictor.to(device)

    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=V3_CONFIG['learning_rate'],
        weight_decay=V3_CONFIG['weight_decay'],
    )

    # Negative sampler
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
            model, predictor, optimizer, data,
            train_edges, negative_sampler, device
        )

        # Validate every 5 epochs
        if epoch % 5 == 0:
            val_metrics = evaluate_full(
                model, predictor, data, val_edges,
                tier_r_lookup, device
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

    # Final evaluation (filtered)
    train_full = evaluate_full(model, predictor, data, train_edges, tier_r_lookup, device)
    val_full = evaluate_full(model, predictor, data, val_edges, tier_r_lookup, device)
    test_full = evaluate_full(model, predictor, data, test_edges, tier_r_lookup, device)

    # Experimental subset evaluation
    train_exp = evaluate_experimental(model, predictor, data, train_edges, tier_r_lookup, device)
    val_exp = evaluate_experimental(model, predictor, data, val_edges, tier_r_lookup, device)
    test_exp = evaluate_experimental(model, predictor, data, test_edges, tier_r_lookup, device)

    return {
        'train_full': train_full,
        'val_full': val_full,
        'test_full': test_full,
        'train_exp': train_exp,
        'val_exp': val_exp,
        'test_exp': test_exp,
        'best_epoch': best_epoch,
    }


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="v3 Training Pipeline (Filtered Evaluation)"
    )
    parser.add_argument(
        '--graph', type=str,
        default=str(DATA_PATHS['graph_v3']),
        help='Path to v3 graph file'
    )
    parser.add_argument(
        '--seeds', type=str, default='42,123,456',
        help='Comma-separated random seeds'
    )
    parser.add_argument(
        '--num-layers', type=int, default=2,
        help='Number of GNN layers'
    )
    parser.add_argument(
        '--epochs', type=int, default=V3_CONFIG['epochs'],
        help='Max training epochs'
    )
    parser.add_argument(
        '--patience', type=int, default=V3_CONFIG['patience'],
        help='Early stopping patience'
    )
    parser.add_argument(
        '--neg-strategy', type=str, default='ec_class',
        choices=['random', 'hard', 'mixed', 'ec_class'],
        help='Negative sampling strategy'
    )
    parser.add_argument(
        '--no-rxn-neighbor', action='store_true',
        help='Ablation: disable rxn_neighbor edges'
    )
    parser.add_argument(
        '--no-tf-edges', action='store_true',
        help='Ablation: disable TF edges'
    )
    parser.add_argument(
        '--output', type=str, default=None,
        help='Output file for results'
    )

    args = parser.parse_args()

    # Parse seeds
    seeds = [int(s.strip()) for s in args.seeds.split(',')]

    device = get_device()
    use_rxn = not args.no_rxn_neighbor
    use_tf = not args.no_tf_edges

    print("=" * 70)
    print("v3 Training Pipeline (Filtered Evaluation)")
    print("=" * 70)
    print(f"Device: {device}")
    print(f"Graph: {args.graph}")
    print(f"Seeds: {seeds}")
    print(f"Layers: {args.num_layers}")
    print(f"rxn_neighbor: {'ON' if use_rxn else 'OFF (ablation)'}")
    print(f"TF edges: {'ON' if use_tf else 'OFF (ablation)'}")
    print("=" * 70)

    # Load data info
    if not Path(args.graph).exists():
        print(f"Error: Graph file not found: {args.graph}")
        print("Run data_pipeline_v3.py first to build the graph.")
        return

    data = torch.load(args.graph)
    print(f"Enzymes: {data['Enzyme'].num_nodes}")
    print(f"Metabolites: {data['Metabolite'].num_nodes}")

    if ('Enzyme', 'catalyzes_R', 'Metabolite') in data.edge_types:
        print(f"Tier-R edges: {data['Enzyme', 'catalyzes_R', 'Metabolite'].edge_index.shape[1]}")
    if ('Metabolite', 'rxn_neighbor', 'Metabolite') in data.edge_types:
        print(f"rxn_neighbor edges: {data['Metabolite', 'rxn_neighbor', 'Metabolite'].edge_index.shape[1]}")

    exp_count = data['Metabolite'].is_experimental.sum().item()
    print(f"Experimental metabolites: {exp_count}")
    print()

    # Run experiments
    all_results = {
        'test_full': [], 'test_exp': [],
        'val_full': [], 'val_exp': [],
    }

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
            use_rxn_neighbor=use_rxn,
            use_tf_edges=use_tf,
        )

        all_results['test_full'].append(result['test_full'])
        all_results['test_exp'].append(result['test_exp'])
        all_results['val_full'].append(result['val_full'])
        all_results['val_exp'].append(result['val_exp'])

        print(
            f"Full H@20={result['test_full']['hits@20']:.3f}, "
            f"MRR={result['test_full']['mrr']:.4f} | "
            f"Exp H@3={result['test_exp'].get('hits@3', 0):.3f} "
            f"(epoch {result['best_epoch']})"
        )

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY (Test Set - Filtered Evaluation)")
    print("=" * 70)

    print("\n--- Full Candidate Set ---")
    test_full = all_results['test_full']
    for metric in ['hits@1', 'hits@10', 'hits@20', 'mrr']:
        if metric in test_full[0]:
            values = np.array([r[metric] for r in test_full])
            if 'hits' in metric:
                print(f"{metric}: {values.mean()*100:.2f}% +/- {values.std()*100:.2f}%")
            else:
                print(f"{metric}: {values.mean():.4f} +/- {values.std():.4f}")

    print("\n--- Experimental Subset ---")
    test_exp = all_results['test_exp']
    for metric in ['hits@1', 'hits@3', 'mrr', 'count']:
        if metric in test_exp[0]:
            values = np.array([r[metric] for r in test_exp])
            if metric == 'count':
                print(f"{metric}: {values.mean():.1f}")
            elif 'hits' in metric:
                print(f"{metric}: {values.mean()*100:.2f}% +/- {values.std()*100:.2f}%")
            else:
                print(f"{metric}: {values.mean():.4f} +/- {values.std():.4f}")

    # Save results
    output_path = args.output or str(DATA_PATHS['gnn_results'] / 'train_v3_results.tsv')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    ablation_tag = ""
    if not use_rxn:
        ablation_tag += "_no_rxn"
    if not use_tf:
        ablation_tag += "_no_tf"

    with open(output_path, 'w') as f:
        f.write("seed\teval_type\thits@1\thits@3\thits@10\thits@20\thits@50\tmrr\tcount\tablation\n")
        for i, seed in enumerate(seeds):
            for eval_type in ['test_full', 'test_exp', 'val_full', 'val_exp']:
                r = all_results.get(eval_type, [{}])[i] if i < len(all_results.get(eval_type, [])) else {}
                f.write(f"{seed}\t{eval_type}\t")
                f.write(f"{r.get('hits@1', 0):.4f}\t{r.get('hits@3', 0):.4f}\t")
                f.write(f"{r.get('hits@10', 0):.4f}\t{r.get('hits@20', 0):.4f}\t")
                f.write(f"{r.get('hits@50', 0):.4f}\t{r.get('mrr', 0):.4f}\t")
                f.write(f"{r.get('count', 0)}\t{ablation_tag or 'full'}\n")

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
