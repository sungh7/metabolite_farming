"""
Distance-based Performance Evaluation

Evaluates Hits@20 stratified by Reaction Distance (3, 4, ≥5).
Shows how GNN generalizes to increasingly distant metabolite-enzyme pairs.
"""

import torch
import numpy as np
from sklearn.metrics import average_precision_score
import os
import sys
sys.path.append(os.getcwd())

from src.model import HGT, LinkPredictor

def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

def get_reaction_distance(enz_idx, met_idx, enz_pathways, met_pathways):
    """
    Simulate reaction distance based on pathway membership.
    In real implementation, this would use KEGG reaction graph BFS.
    
    Distance simulation:
    - Same pathway, adjacent index: 1-2
    - Same pathway, distant index: 3-4
    - Different pathway: ≥5
    """
    enz_path = enz_pathways.get(enz_idx, 'Other')
    met_path = met_pathways.get(met_idx, 'Other')
    
    if enz_path == met_path and enz_path != 'Other':
        # Same pathway - distance based on index proximity
        # Metabolites 0-19: Phenylpropanoid, 20-49: Flavonoid
        if enz_path == 'Phenylpropanoid':
            base = 0
            size = 20
        else:
            base = 20
            size = 30
        
        rel_pos = (met_idx - base) % size
        if rel_pos < 5:
            return 3  # Close in pathway
        elif rel_pos < 15:
            return 4  # Medium distance
        else:
            return 5  # Far in pathway
    else:
        return 6  # Cross-pathway (very distant)

def evaluate_by_distance(data, model, predictor, test_edges, test_enzymes, device):
    """Evaluate Hits@20 stratified by reaction distance."""
    num_enzymes = data['Enzyme'].num_nodes
    num_metabolites = data['Metabolite'].num_nodes
    
    # Pathway assignments (same as bipartite_builder)
    np.random.seed(42)
    met_pathways = {}
    for i in range(20): met_pathways[i] = 'Phenylpropanoid'
    for i in range(20, 50): met_pathways[i] = 'Flavonoid'
    for i in range(50, num_metabolites): met_pathways[i] = 'Other'
    
    enz_pathways = {}
    for i in range(num_enzymes):
        r = np.random.rand()
        if r < 0.05: p = 'Phenylpropanoid'
        elif r < 0.15: p = 'Flavonoid'
        else: p = 'Other'
        enz_pathways[i] = p
    
    model.eval()
    test_enzymes_np = test_enzymes.cpu().numpy()
    
    # Group test edges by distance
    distance_groups = {3: [], 4: [], 5: []}  # 5 includes ≥5
    
    test_src = test_edges[0].cpu().numpy()
    test_dst = test_edges[1].cpu().numpy()
    
    for i in range(len(test_src)):
        enz_idx = int(test_src[i])
        met_idx = int(test_dst[i])
        dist = get_reaction_distance(enz_idx, met_idx, enz_pathways, met_pathways)
        
        if dist == 3:
            distance_groups[3].append((enz_idx, met_idx))
        elif dist == 4:
            distance_groups[4].append((enz_idx, met_idx))
        else:  # ≥5
            distance_groups[5].append((enz_idx, met_idx))
    
    results = {}
    
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        for dist, pairs in distance_groups.items():
            if len(pairs) == 0:
                results[dist] = (0.0, 0)
                continue
            
            hits_20 = 0
            count = 0
            
            # Group by metabolite
            met_to_enzs = {}
            for enz_idx, met_idx in pairs:
                if met_idx not in met_to_enzs:
                    met_to_enzs[met_idx] = []
                met_to_enzs[met_idx].append(enz_idx)
            
            for met_idx, true_enzs in met_to_enzs.items():
                true_enzs_set = set(true_enzs)
                
                # Score all test enzymes
                candidate_enzs = test_enzymes
                num_cands = candidate_enzs.size(0)
                eval_src = candidate_enzs
                eval_dst = torch.full((num_cands,), met_idx, device=device)
                eval_edges = torch.stack([eval_src, eval_dst])
                
                scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()
                is_true = torch.tensor([int(e) in true_enzs_set for e in candidate_enzs.cpu().numpy()], device=device)
                
                if is_true.sum() == 0:
                    continue
                
                sorted_indices = torch.argsort(scores, descending=True)
                sorted_labels = is_true[sorted_indices]
                
                if sorted_labels[:20].sum() > 0:
                    hits_20 += 1
                count += 1
            
            results[dist] = (hits_20 / count if count > 0 else 0.0, count)
    
    return results

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    set_seed(42)
    
    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    num_enzymes = data['Enzyme'].num_nodes
    
    # Same split as training (90/10)
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    test_enz_mask = ~train_enz_mask
    
    # Get train data for model
    train_data = data.clone()
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_enz_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1]
    rev_mask = train_enz_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    test_enzymes = indices[split:]
    
    # Train model
    print("Training model...")
    model = HGT(train_data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    for epoch in range(20):
        model.train()
        optimizer.zero_grad()
        
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        pos_edge_index = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        
        num_pos = pos_edge_index.size(1)
        num_metabolites = train_data['Metabolite'].num_nodes
        neg_src = pos_edge_index[0]
        offset = torch.randint(1, 6, (num_pos,), device=device) * (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge_index[1] + offset) % num_metabolites
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
    
    print("Evaluating by distance...")
    results = evaluate_by_distance(train_data, model, predictor, test_edges, test_enzymes, device)
    
    print("\n" + "="*50)
    print("DISTANCE-BASED PERFORMANCE")
    print("="*50)
    for dist in [3, 4, 5]:
        hits20, count = results[dist]
        label = f"≥{dist}" if dist == 5 else str(dist)
        print(f"Distance {label}: Hits@20 = {hits20:.4f} (n={count})")
    print("="*50)
    
    # Save results
    os.makedirs('results/gnn', exist_ok=True)
    with open('results/gnn/distance_analysis.txt', 'w') as f:
        f.write("Distance-based Performance Analysis\n")
        f.write("="*40 + "\n")
        for dist in [3, 4, 5]:
            hits20, count = results[dist]
            label = f"≥{dist}" if dist == 5 else str(dist)
            f.write(f"Distance {label}: Hits@20 = {hits20:.4f} (n={count})\n")
    print("Results saved to results/gnn/distance_analysis.txt")

if __name__ == "__main__":
    main()
