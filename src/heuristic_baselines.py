"""
Heuristic Baselines for Link Prediction

Implements:
1. Adamic-Adar (Hetero): Type-aware neighbor overlap scoring
2. Resource Allocation (Hetero): Type-aware resource distribution

These are classic link prediction heuristics adapted for heterogeneous graphs.
"""

import torch
import numpy as np
import os
import sys
sys.path.append(os.getcwd())

def build_neighbor_dict(data, device):
    """Build neighbor dictionary for each node type."""
    neighbors = {}
    
    num_enzymes = data['Enzyme'].num_nodes
    num_metabolites = data['Metabolite'].num_nodes
    num_tfs = data['TF'].num_nodes if 'TF' in data.node_types else 0
    
    # Initialize
    neighbors['Enzyme'] = {i: set() for i in range(num_enzymes)}
    neighbors['Metabolite'] = {i: set() for i in range(num_metabolites)}
    if num_tfs > 0:
        neighbors['TF'] = {i: set() for i in range(num_tfs)}
    
    # Enzyme-Metabolite
    if ('Enzyme', 'catalyzes', 'Metabolite') in data.edge_types:
        edge = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        for i in range(edge.size(1)):
            enz = edge[0, i].item()
            met = edge[1, i].item()
            neighbors['Enzyme'][enz].add(('Metabolite', met))
            neighbors['Metabolite'][met].add(('Enzyme', enz))
    
    return neighbors

def adamic_adar_score(neighbors, src_type, src_idx, dst_type, dst_idx):
    """
    Adamic-Adar index: sum of 1/log(degree) for common neighbors.
    """
    src_neighbors = neighbors.get(src_type, {}).get(src_idx, set())
    dst_neighbors = neighbors.get(dst_type, {}).get(dst_idx, set())
    
    common = src_neighbors & dst_neighbors
    
    if len(common) == 0:
        return 0.0
    
    score = 0.0
    for neighbor_type, neighbor_idx in common:
        degree = len(neighbors.get(neighbor_type, {}).get(neighbor_idx, set()))
        if degree > 1:
            score += 1.0 / np.log(degree)
    
    return score

def resource_allocation_score(neighbors, src_type, src_idx, dst_type, dst_idx):
    """
    Resource Allocation: sum of 1/degree for common neighbors.
    """
    src_neighbors = neighbors.get(src_type, {}).get(src_idx, set())
    dst_neighbors = neighbors.get(dst_type, {}).get(dst_idx, set())
    
    common = src_neighbors & dst_neighbors
    
    if len(common) == 0:
        return 0.0
    
    score = 0.0
    for neighbor_type, neighbor_idx in common:
        degree = len(neighbors.get(neighbor_type, {}).get(neighbor_idx, set()))
        if degree > 0:
            score += 1.0 / degree
    
    return score

def evaluate_heuristic(data, heuristic_fn, device):
    """Evaluate heuristic baseline on Enzyme-Metabolite link prediction."""
    num_enzymes = data['Enzyme'].num_nodes
    
    # Build neighbor dict from training edges
    np.random.seed(42)
    indices = np.random.permutation(num_enzymes)
    split = int(0.9 * num_enzymes)
    train_enzymes = set(indices[:split])
    test_enzymes = set(indices[split:])
    
    # Create training graph
    train_data = data.clone()
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    mask = torch.tensor([int(src) in train_enzymes for src in edge_index[0].cpu().numpy()])
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index[:, mask]
    
    neighbors = build_neighbor_dict(train_data, device)
    
    # Get test edges
    test_mask = ~mask
    test_edges = edge_index[:, test_mask]
    test_src = test_edges[0].cpu().numpy()
    test_dst = test_edges[1]
    unique_mets = torch.unique(test_dst)
    
    hits_20 = 0
    reciprocal_ranks = []
    count = 0
    
    for met_idx in unique_mets:
        met_mask = (test_dst == met_idx)
        true_enzs = set(test_src[met_mask.cpu().numpy()])
        
        # Score all test enzymes
        test_enz_list = list(test_enzymes)
        scores = []
        for enz in test_enz_list:
            score = heuristic_fn(neighbors, 'Enzyme', enz, 'Metabolite', met_idx.item())
            scores.append(score)
        
        scores = np.array(scores)
        sorted_indices = np.argsort(-scores)
        
        is_true = np.array([test_enz_list[i] in true_enzs for i in sorted_indices])
        
        if is_true.sum() == 0:
            continue
        
        if is_true[:20].sum() > 0:
            hits_20 += 1
        
        true_positions = np.where(is_true)[0]
        if len(true_positions) > 0:
            reciprocal_ranks.append(1.0 / (true_positions[0] + 1))
        
        count += 1
    
    hits20_score = hits_20 / count if count > 0 else 0.0
    mrr_score = np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0
    
    return {'hits20': hits20_score, 'mrr': mrr_score}

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    heuristics = {
        'Adamic-Adar': adamic_adar_score,
        'Resource-Allocation': resource_allocation_score
    }
    
    results = {}
    
    for name, heuristic_fn in heuristics.items():
        print(f"\nEvaluating {name}...")
        result = evaluate_heuristic(data, heuristic_fn, device)
        results[name] = result
        print(f"  Hits@20: {result['hits20']:.4f}, MRR: {result['mrr']:.4f}")
    
    print("\n" + "="*50)
    print("HEURISTIC BASELINE RESULTS")
    print("="*50)
    for name, result in results.items():
        print(f"{name}: Hits@20={result['hits20']:.4f}, MRR={result['mrr']:.4f}")
    
    # Save
    os.makedirs('results/gnn', exist_ok=True)
    with open('results/gnn/heuristic_baselines.tsv', 'w') as f:
        f.write("Method\tHits20\tMRR\n")
        for name, result in results.items():
            f.write(f"{name}\t{result['hits20']:.4f}\t{result['mrr']:.4f}\n")
    print("\nResults saved to results/gnn/heuristic_baselines.tsv")

if __name__ == "__main__":
    main()
