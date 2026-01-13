"""
GNNShap: Edge-Level Shapley Value Explainability

Computes Shapley-based edge importance for link predictions.
Shows which edges in the neighborhood contribute most to predictions.
"""

import torch
import numpy as np
import json
import os
import sys
from itertools import combinations
sys.path.append(os.getcwd())

from src.model import HGT, LinkPredictor

def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)

def get_k_hop_edges(data, node_type, node_idx, k=2, max_edges=20):
    """
    Get edges within k-hops of a node.
    Returns list of (edge_type, edge_idx) tuples.
    """
    edges = []
    visited_nodes = {(node_type, node_idx)}
    frontier = {(node_type, node_idx)}
    
    for hop in range(k):
        new_frontier = set()
        for nt, ni in frontier:
            # Find edges connected to this node
            for edge_type in data.edge_types:
                edge_index = data[edge_type].edge_index
                
                if edge_type[0] == nt:
                    # Outgoing edges
                    mask = (edge_index[0] == ni)
                    for edge_idx in torch.where(mask)[0][:5]:  # Limit per edge type
                        edges.append((edge_type, edge_idx.item()))
                        neighbor = (edge_type[2], edge_index[1, edge_idx].item())
                        if neighbor not in visited_nodes:
                            new_frontier.add(neighbor)
                            visited_nodes.add(neighbor)
                
                if edge_type[2] == nt:
                    # Incoming edges
                    mask = (edge_index[1] == ni)
                    for edge_idx in torch.where(mask)[0][:5]:
                        edges.append((edge_type, edge_idx.item()))
                        neighbor = (edge_type[0], edge_index[0, edge_idx].item())
                        if neighbor not in visited_nodes:
                            new_frontier.add(neighbor)
                            visited_nodes.add(neighbor)
        
        frontier = new_frontier
        if len(edges) >= max_edges:
            break
    
    return edges[:max_edges]

def mask_edges(data, edges_to_mask):
    """
    Create a copy of data with specified edges masked (removed).
    """
    masked_data = data.clone()
    
    for edge_type, edge_idx in edges_to_mask:
        edge_index = data[edge_type].edge_index
        mask = torch.ones(edge_index.size(1), dtype=torch.bool, device=edge_index.device)
        mask[edge_idx] = False
        masked_data[edge_type].edge_index = edge_index[:, mask]
    
    return masked_data

def compute_shapley_values(data, model, predictor, src_type, src_idx, dst_type, dst_idx, 
                           edges, device, n_samples=50):
    """
    Compute Shapley values for each edge using sampling approximation.
    """
    model.eval()
    
    n_edges = len(edges)
    if n_edges == 0:
        return {}
    
    shapley_values = {i: 0.0 for i in range(n_edges)}
    
    # Monte Carlo sampling for Shapley approximation
    for _ in range(n_samples):
        # Random permutation
        perm = np.random.permutation(n_edges)
        
        prev_score = None
        for i, edge_idx in enumerate(perm):
            # Coalition before adding this edge
            coalition = [edges[perm[j]] for j in range(i)]
            
            # Score without this edge
            if i == 0:
                # Empty coalition - use all edges
                with torch.no_grad():
                    x_dict = model(data.x_dict, data.edge_index_dict)
                    edge = torch.tensor([[src_idx], [dst_idx]], device=device)
                    prev_score = predictor(x_dict[src_type], x_dict[dst_type], edge).sigmoid().item()
            
            # Score with this edge added
            edges_before = [edges[perm[j]] for j in range(i + 1)]
            edges_to_mask = [e for e in edges if e not in edges_before]
            
            masked_data = mask_edges(data, edges_to_mask)
            
            with torch.no_grad():
                x_dict = model(masked_data.x_dict, masked_data.edge_index_dict)
                edge = torch.tensor([[src_idx], [dst_idx]], device=device)
                curr_score = predictor(x_dict[src_type], x_dict[dst_type], edge).sigmoid().item()
            
            # Marginal contribution
            marginal = curr_score - prev_score
            shapley_values[edge_idx] += marginal
            prev_score = curr_score
    
    # Average
    for k in shapley_values:
        shapley_values[k] /= n_samples
    
    return shapley_values

def explain_prediction(data, model, predictor, src_type, src_idx, dst_type, dst_idx, device):
    """
    Generate edge-level explanation for a single prediction.
    """
    # Get edges in neighborhood
    src_edges = get_k_hop_edges(data, src_type, src_idx, k=2, max_edges=15)
    dst_edges = get_k_hop_edges(data, dst_type, dst_idx, k=2, max_edges=15)
    
    # Combine and deduplicate
    all_edges = list(set(src_edges + dst_edges))
    
    if len(all_edges) == 0:
        return {'edges': [], 'shapley_values': []}
    
    # Compute Shapley values
    shapley = compute_shapley_values(
        data, model, predictor, src_type, src_idx, dst_type, dst_idx,
        all_edges, device, n_samples=30
    )
    
    # Format results
    results = []
    for i, edge in enumerate(all_edges):
        results.append({
            'edge_type': str(edge[0]),
            'edge_idx': edge[1],
            'shapley_value': shapley.get(i, 0.0)
        })
    
    # Sort by absolute importance
    results.sort(key=lambda x: -abs(x['shapley_value']))
    
    return results

def compute_explanation_fidelity(data, model, predictor, src_type, src_idx, dst_type, dst_idx,
                                  top_edges, device):
    """
    Compute explanation fidelity: score drop when top edges are removed.
    """
    model.eval()
    
    # Original score
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        edge = torch.tensor([[src_idx], [dst_idx]], device=device)
        original_score = predictor(x_dict[src_type], x_dict[dst_type], edge).sigmoid().item()
    
    # Score with top edges removed
    edges_to_mask = [(tuple(e['edge_type'].strip("()' ").replace("'", "").split(", ")), e['edge_idx']) 
                     for e in top_edges[:5]]
    
    # Skip if edge parsing fails
    valid_edges = []
    for et, ei in edges_to_mask:
        if len(et) == 3:
            valid_edges.append((tuple(et), ei))
    
    if not valid_edges:
        return {'original': original_score, 'after_removal': original_score, 'fidelity': 0.0}
    
    masked_data = mask_edges(data, valid_edges)
    
    with torch.no_grad():
        x_dict = model(masked_data.x_dict, masked_data.edge_index_dict)
        edge = torch.tensor([[src_idx], [dst_idx]], device=device)
        perturbed_score = predictor(x_dict[src_type], x_dict[dst_type], edge).sigmoid().item()
    
    fidelity = original_score - perturbed_score
    
    return {
        'original': original_score,
        'after_removal': perturbed_score,
        'fidelity': fidelity
    }

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    set_seed(42)
    
    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    # Train model
    print("Training model...")
    model = HGT(data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    for epoch in range(20):
        model.train()
        optimizer.zero_grad()
        
        x_dict = model(data.x_dict, data.edge_index_dict)
        pos_edge = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        
        num_pos = pos_edge.size(1)
        num_met = data['Metabolite'].num_nodes
        neg_src = pos_edge[0]
        offset = torch.randint(1, 6, (num_pos,), device=device) * (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge[1] + offset) % num_met
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
    
    # Get top-5 predictions
    print("\nGenerating GNNShap explanations...")
    model.eval()
    
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
    
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes, device=device)
    test_enzymes = indices[int(0.9 * num_enzymes):][:10]
    
    # Find high-scoring predictions
    all_preds = []
    for enz_idx in test_enzymes:
        for met_idx in range(min(30, data['Metabolite'].num_nodes)):
            edge = torch.stack([enz_idx.unsqueeze(0), torch.tensor([met_idx], device=device)])
            score = predictor(x_dict['Enzyme'], x_dict['Metabolite'], edge).sigmoid().item()
            all_preds.append({'enzyme': enz_idx.item(), 'metabolite': met_idx, 'score': score})
    
    all_preds.sort(key=lambda x: -x['score'])
    top_preds = all_preds[:5]
    
    # Generate explanations
    explanations = []
    
    for pred in top_preds:
        print(f"\nExplaining Enzyme[{pred['enzyme']}] → Metabolite[{pred['metabolite']}] (Score: {pred['score']:.4f})")
        
        edge_importance = explain_prediction(
            data, model, predictor, 'Enzyme', pred['enzyme'], 'Metabolite', pred['metabolite'], device
        )
        
        # Compute fidelity
        fidelity = compute_explanation_fidelity(
            data, model, predictor, 'Enzyme', pred['enzyme'], 'Metabolite', pred['metabolite'],
            edge_importance, device
        )
        
        explanations.append({
            'prediction': pred,
            'top_edges': edge_importance[:10],
            'fidelity': fidelity
        })
        
        print(f"  Fidelity: {fidelity['fidelity']:.4f} (Original: {fidelity['original']:.4f} → After: {fidelity['after_removal']:.4f})")
        print("  Top-3 edges:")
        for i, e in enumerate(edge_importance[:3], 1):
            print(f"    {i}. {e['edge_type']} (Shapley: {e['shapley_value']:.4f})")
    
    # Save results
    os.makedirs('results/explainability', exist_ok=True)
    
    with open('results/explainability/gnnshap_explanations.json', 'w') as f:
        json.dump(explanations, f, indent=2, default=float)
    
    # Summary table
    with open('results/explainability/gnnshap_summary.tsv', 'w') as f:
        f.write("Rank\tEnzyme\tMetabolite\tScore\tFidelity\tTop_Edge\tShapley\n")
        for i, exp in enumerate(explanations, 1):
            pred = exp['prediction']
            fid = exp['fidelity']['fidelity']
            top_e = exp['top_edges'][0] if exp['top_edges'] else {'edge_type': 'N/A', 'shapley_value': 0}
            f.write(f"{i}\t{pred['enzyme']}\t{pred['metabolite']}\t{pred['score']:.4f}\t{fid:.4f}\t{top_e['edge_type']}\t{top_e['shapley_value']:.4f}\n")
    
    print("\n" + "="*50)
    print("GNNShap COMPLETE")
    print("="*50)
    print("Results saved to results/explainability/gnnshap_*.{json,tsv}")

if __name__ == "__main__":
    main()
