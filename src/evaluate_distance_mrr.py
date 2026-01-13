"""
Enhanced Distance-based Evaluation with MRR

Adds Mean Reciprocal Rank (MRR) to complement Hits@20 for small sample sizes.
MRR is less sensitive to sample size and provides finer-grained ranking quality.
"""

import torch
import numpy as np
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
    """Simulate reaction distance based on pathway membership."""
    enz_path = enz_pathways.get(enz_idx, 'Other')
    met_path = met_pathways.get(met_idx, 'Other')
    
    if enz_path == met_path and enz_path != 'Other':
        if enz_path == 'Phenylpropanoid':
            base, size = 0, 20
        else:
            base, size = 20, 30
        
        rel_pos = (met_idx - base) % size
        if rel_pos < 5:
            return 3
        elif rel_pos < 15:
            return 4
        else:
            return 5
    else:
        return 6

def evaluate_with_mrr(data, model, predictor, test_edges, test_enzymes, device):
    """Evaluate Hits@20 and MRR stratified by reaction distance."""
    num_enzymes = data['Enzyme'].num_nodes
    num_metabolites = data['Metabolite'].num_nodes
    
    # Pathway assignments
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
    distance_groups = {3: [], 4: [], 5: []}
    
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
        else:
            distance_groups[5].append((enz_idx, met_idx))
    
    results = {}
    
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        for dist, pairs in distance_groups.items():
            if len(pairs) == 0:
                results[dist] = {'hits20': 0.0, 'mrr': 0.0, 'n': 0}
                continue
            
            hits_20 = 0
            reciprocal_ranks = []
            count = 0
            
            met_to_enzs = {}
            for enz_idx, met_idx in pairs:
                if met_idx not in met_to_enzs:
                    met_to_enzs[met_idx] = []
                met_to_enzs[met_idx].append(enz_idx)
            
            for met_idx, true_enzs in met_to_enzs.items():
                true_enzs_set = set(true_enzs)
                
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
                
                # Hits@20
                if sorted_labels[:20].sum() > 0:
                    hits_20 += 1
                
                # MRR: find rank of first true positive
                true_positions = torch.where(sorted_labels)[0]
                if len(true_positions) > 0:
                    first_true_rank = true_positions[0].item() + 1  # 1-indexed
                    reciprocal_ranks.append(1.0 / first_true_rank)
                
                count += 1
            
            hits20_score = hits_20 / count if count > 0 else 0.0
            mrr_score = np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0
            
            results[dist] = {
                'hits20': hits20_score,
                'mrr': mrr_score,
                'n': count
            }
    
    return results

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    set_seed(42)
    
    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    num_enzymes = data['Enzyme'].num_nodes
    
    # Split
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    
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
        num_met = train_data['Metabolite'].num_nodes
        neg_src = pos_edge_index[0]
        offset = torch.randint(1, 6, (num_pos,), device=device) * (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge_index[1] + offset) % num_met
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
    
    print("Evaluating with MRR...")
    results = evaluate_with_mrr(train_data, model, predictor, test_edges, test_enzymes, device)
    
    print("\n" + "="*60)
    print("DISTANCE-BASED PERFORMANCE (with MRR)")
    print("="*60)
    print(f"{'Distance':<12} {'Hits@20':<12} {'MRR':<12} {'N':<8}")
    print("-"*44)
    for dist in [3, 4, 5]:
        label = f"≥{dist}" if dist == 5 else str(dist)
        r = results[dist]
        print(f"{label:<12} {r['hits20']:.4f}       {r['mrr']:.4f}       {r['n']}")
    print("="*60)
    
    # Save results
    os.makedirs('results/gnn', exist_ok=True)
    with open('results/gnn/distance_mrr_analysis.txt', 'w') as f:
        f.write("Distance-based Performance with MRR\n")
        f.write("="*50 + "\n\n")
        f.write(f"{'Distance':<12} {'Hits@20':<12} {'MRR':<12} {'N':<8}\n")
        f.write("-"*44 + "\n")
        for dist in [3, 4, 5]:
            label = f"≥{dist}" if dist == 5 else str(dist)
            r = results[dist]
            f.write(f"{label:<12} {r['hits20']:.4f}       {r['mrr']:.4f}       {r['n']}\n")
    print("Results saved to results/gnn/distance_mrr_analysis.txt")

if __name__ == "__main__":
    main()
