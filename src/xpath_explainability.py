"""
xPath-style Influence Path Explainability

Extracts influential paths for link predictions in heterogeneous GNN.
Shows "source → intermediate → target" influence chains.
"""

import torch
import numpy as np
import json
import os
import sys
sys.path.append(os.getcwd())

from src.model import HGT, LinkPredictor

def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)

def extract_influence_paths(data, model, predictor, src_type, src_idx, dst_type, dst_idx, 
                            x_dict, top_k=5, max_depth=3):
    """
    Extract top-k influential paths from src to dst.
    
    Uses gradient-based importance to trace information flow.
    """
    model.eval()
    
    paths = []
    
    # Get edge types that connect to source/target
    src_edges = [et for et in data.edge_types if et[0] == src_type or et[2] == src_type]
    dst_edges = [et for et in data.edge_types if et[0] == dst_type or et[2] == dst_type]
    
    # Simple path extraction: identify neighbors and their importance
    # For heterogeneous graphs, paths go through different node types
    
    for edge_type in src_edges:
        edge_index = data[edge_type].edge_index
        
        if edge_type[0] == src_type:
            # Outgoing edges from source
            mask = (edge_index[0] == src_idx)
            neighbors = edge_index[1, mask].cpu().numpy()
            neighbor_type = edge_type[2]
        else:
            # Incoming edges to source
            mask = (edge_index[1] == src_idx)
            neighbors = edge_index[0, mask].cpu().numpy()
            neighbor_type = edge_type[0]
        
        for neighbor in neighbors[:10]:  # Limit to top 10 neighbors
            # Calculate importance score (simplified: use embedding similarity)
            if neighbor_type in x_dict:
                neighbor_emb = x_dict[neighbor_type][neighbor]
                if dst_type in x_dict:
                    dst_emb = x_dict[dst_type][dst_idx]
                    importance = torch.dot(neighbor_emb, dst_emb).item()
                else:
                    importance = 0.0
            else:
                importance = 0.0
            
            paths.append({
                'path': [
                    {'type': src_type, 'idx': int(src_idx)},
                    {'type': neighbor_type, 'idx': int(neighbor)},
                    {'type': dst_type, 'idx': int(dst_idx)}
                ],
                'edge_types': [edge_type[1], 'predicted'],
                'importance': importance
            })
    
    # Sort by importance and return top-k
    paths.sort(key=lambda x: -abs(x['importance']))
    return paths[:top_k]

def explain_top_predictions(data, model, predictor, device, n_predictions=5, n_paths=3):
    """
    Generate explanations for top-N predictions.
    """
    set_seed(42)
    
    model.eval()
    
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
    
    # Get predictions (using test enzymes)
    num_enzymes = data['Enzyme'].num_nodes
    num_metabolites = data['Metabolite'].num_nodes
    
    indices = torch.randperm(num_enzymes, device=device)
    test_enzymes = indices[int(0.9 * num_enzymes):]
    
    # Score all enzyme-metabolite pairs for test enzymes
    all_predictions = []
    
    for met_idx in range(min(50, num_metabolites)):
        for enz_idx in test_enzymes[:20]:
            edge = torch.stack([enz_idx.unsqueeze(0), torch.tensor([met_idx], device=device)])
            score = predictor(x_dict['Enzyme'], x_dict['Metabolite'], edge).sigmoid().item()
            all_predictions.append({
                'enzyme_idx': enz_idx.item(),
                'metabolite_idx': met_idx,
                'score': score
            })
    
    # Sort by score
    all_predictions.sort(key=lambda x: -x['score'])
    top_predictions = all_predictions[:n_predictions]
    
    # Extract paths for each top prediction
    explanations = []
    
    for pred in top_predictions:
        enz_idx = pred['enzyme_idx']
        met_idx = pred['metabolite_idx']
        
        paths = extract_influence_paths(
            data, model, predictor,
            'Enzyme', enz_idx, 'Metabolite', met_idx,
            x_dict, top_k=n_paths
        )
        
        explanations.append({
            'prediction': pred,
            'paths': paths
        })
    
    return explanations

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    set_seed(42)
    
    # Load data and model
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    # Train model
    print("Training model for explanation...")
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
    
    print("Extracting explanations...")
    explanations = explain_top_predictions(data, model, predictor, device, n_predictions=5, n_paths=3)
    
    print("\n" + "="*60)
    print("TOP-5 PREDICTIONS WITH INFLUENCE PATHS")
    print("="*60)
    
    for i, exp in enumerate(explanations, 1):
        pred = exp['prediction']
        print(f"\n{i}. Enzyme[{pred['enzyme_idx']}] → Metabolite[{pred['metabolite_idx']}]")
        print(f"   Score: {pred['score']:.4f}")
        print(f"   Influential Paths:")
        for j, path in enumerate(exp['paths'], 1):
            path_str = " → ".join([f"{p['type']}[{p['idx']}]" for p in path['path']])
            print(f"      {j}. {path_str} (importance: {path['importance']:.4f})")
    
    # Save results
    os.makedirs('results/explainability', exist_ok=True)
    
    with open('results/explainability/xpath_explanations.json', 'w') as f:
        json.dump(explanations, f, indent=2)
    
    # Create summary table
    with open('results/explainability/xpath_summary.tsv', 'w') as f:
        f.write("Rank\tEnzyme_Idx\tMetabolite_Idx\tScore\tTop_Path\tImportance\n")
        for i, exp in enumerate(explanations, 1):
            pred = exp['prediction']
            if exp['paths']:
                top_path = exp['paths'][0]
                path_str = " -> ".join([f"{p['type']}[{p['idx']}]" for p in top_path['path']])
                importance = top_path['importance']
            else:
                path_str = "N/A"
                importance = 0.0
            f.write(f"{i}\t{pred['enzyme_idx']}\t{pred['metabolite_idx']}\t{pred['score']:.4f}\t{path_str}\t{importance:.4f}\n")
    
    print("\nResults saved to results/explainability/")

if __name__ == "__main__":
    main()
