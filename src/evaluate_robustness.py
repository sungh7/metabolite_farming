"""
Robustness Evaluation Script (5-Fold CV × 3 Seeds = 15 Runs)

Produces mean ± std and 95% CI for Hits@20 under Hard Negative conditions.
"""

import torch
import numpy as np
from sklearn.model_selection import KFold
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

def train_one_fold(data, train_mask, test_mask, seed, device, epochs=20):
    """Train and evaluate on one fold with given seed."""
    set_seed(seed)
    
    # Clone data
    train_data = data.clone()
    
    # Get enzyme indices
    num_enzymes = data['Enzyme'].num_nodes
    train_enz_mask = train_mask
    test_enz_mask = test_mask
    
    # Mask edges
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_enz_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    
    # Handle reverse edges
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1]
    rev_mask = train_enz_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    # Model
    model = HGT(train_data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    # Get edge weights for MSI-based loss weighting
    edge_weight = data['Enzyme', 'catalyzes', 'Metabolite'].get('edge_weight', None)
    if edge_weight is not None:
        train_weights = edge_weight[mask].to(device)
    else:
        train_weights = None
    
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
        offset = torch.randint(1, 6, (num_pos,), device=device) * (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge_index[1] + offset) % num_metabolites
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        # MSI-Weighted Loss
        if train_weights is not None:
            pos_loss = -(train_weights * torch.log(torch.sigmoid(pos_out) + 1e-15)).mean()
        else:
            pos_loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean()
        neg_loss = -torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss = pos_loss + neg_loss
        
        loss.backward()
        optimizer.step()
    
    # Evaluation
    model.eval()
    test_enzymes = torch.where(test_enz_mask)[0]
    
    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        test_src = test_edges[0]
        test_dst = test_edges[1]
        unique_mets = torch.unique(test_dst)
        
        hits_20 = 0
        count = 0
        
        for met_idx in unique_mets:
            mask = (test_dst == met_idx)
            true_enzs = test_src[mask]
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
            
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1
            count += 1
    
    return hits_20 / count if count > 0 else 0.0

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    num_enzymes = data['Enzyme'].num_nodes
    seeds = [42, 123, 456]
    n_folds = 5
    
    results = []
    
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    enzyme_indices = np.arange(num_enzymes)
    
    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(enzyme_indices)):
        train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
        train_mask[train_idx] = True
        test_mask = ~train_mask
        
        for seed in seeds:
            hits20 = train_one_fold(data, train_mask, test_mask, seed, device)
            results.append(hits20)
            print(f"Fold {fold_idx+1}/5, Seed {seed}: Hits@20 = {hits20:.4f}")
    
    results = np.array(results)
    mean_hits = results.mean()
    std_hits = results.std()
    ci_low = np.percentile(results, 2.5)
    ci_high = np.percentile(results, 97.5)
    
    print("\n" + "="*50)
    print(f"FINAL RESULTS (n={len(results)} runs)")
    print(f"Hits@20: {mean_hits:.4f} ± {std_hits:.4f}")
    print(f"95% CI: [{ci_low:.4f}, {ci_high:.4f}]")
    print("="*50)
    
    # Save results
    os.makedirs('results/gnn', exist_ok=True)
    with open('results/gnn/robustness_results.txt', 'w') as f:
        f.write(f"Robustness Evaluation (5-Fold CV × 3 Seeds = {len(results)} Runs)\n")
        f.write(f"Hits@20: {mean_hits:.4f} ± {std_hits:.4f}\n")
        f.write(f"95% CI: [{ci_low:.4f}, {ci_high:.4f}]\n")
        f.write(f"\nIndividual runs:\n")
        for i, r in enumerate(results):
            f.write(f"  Run {i+1}: {r:.4f}\n")
    print("Results saved to results/gnn/robustness_results.txt")

if __name__ == "__main__":
    main()
