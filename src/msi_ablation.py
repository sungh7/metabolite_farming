"""
MSI Ablation Study

Compares performance with and without MSI-weighted loss.
Demonstrates the value of incorporating annotation confidence into training.
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

def train_and_evaluate(data, train_mask, test_mask, seed, device, use_msi_weights=True, epochs=20):
    """Train and evaluate with or without MSI weights."""
    set_seed(seed)
    
    train_data = data.clone()
    num_enzymes = data['Enzyme'].num_nodes
    
    # Mask edges
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1]
    rev_mask = train_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    # Model
    model = HGT(train_data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    # Get edge weights for MSI-based loss weighting
    edge_weight = data['Enzyme', 'catalyzes', 'Metabolite'].get('edge_weight', None)
    if edge_weight is not None and use_msi_weights:
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
        
        # MSI-Weighted or Uniform Loss
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
    test_enzymes = torch.where(~train_mask)[0]
    
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
    
    results_with_msi = []
    results_without_msi = []
    
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    enzyme_indices = np.arange(num_enzymes)
    
    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(enzyme_indices)):
        train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
        train_mask[train_idx] = True
        test_mask = ~train_mask
        
        for seed in seeds:
            # With MSI weights
            hits20_with = train_and_evaluate(data, train_mask, test_mask, seed, device, use_msi_weights=True)
            results_with_msi.append(hits20_with)
            
            # Without MSI weights
            hits20_without = train_and_evaluate(data, train_mask, test_mask, seed, device, use_msi_weights=False)
            results_without_msi.append(hits20_without)
            
            print(f"Fold {fold_idx+1}/5, Seed {seed}: With MSI = {hits20_with:.4f}, Without MSI = {hits20_without:.4f}")
    
    results_with_msi = np.array(results_with_msi)
    results_without_msi = np.array(results_without_msi)
    
    print("\n" + "="*50)
    print("MSI ABLATION RESULTS (n=15 runs)")
    print("="*50)
    print(f"With MSI Weights:    {results_with_msi.mean():.4f} ± {results_with_msi.std():.4f}")
    print(f"Without MSI Weights: {results_without_msi.mean():.4f} ± {results_without_msi.std():.4f}")
    
    # Statistical comparison
    diff = results_with_msi - results_without_msi
    print(f"\nDifference (With - Without): {diff.mean():.4f} ± {diff.std():.4f}")
    print("="*50)
    
    # Save results
    os.makedirs('results/gnn', exist_ok=True)
    with open('results/gnn/msi_ablation.txt', 'w') as f:
        f.write("MSI Ablation Study (5-Fold CV × 3 Seeds = 15 Runs)\n")
        f.write("="*50 + "\n")
        f.write(f"With MSI Weights:    {results_with_msi.mean():.4f} ± {results_with_msi.std():.4f}\n")
        f.write(f"Without MSI Weights: {results_without_msi.mean():.4f} ± {results_without_msi.std():.4f}\n")
        f.write(f"\nDifference (With - Without): {diff.mean():.4f} ± {diff.std():.4f}\n")
    print("Results saved to results/gnn/msi_ablation.txt")

if __name__ == "__main__":
    main()
