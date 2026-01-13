"""
Inductive Split Evaluation

Implements complete node removal (enzyme-holdout) for cold-start evaluation.
Tests if model can predict for truly unseen enzymes.
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

def create_inductive_split(data, holdout_fraction=0.1, seed=42):
    """
    Create truly inductive split by removing nodes from graph.
    
    Unlike transductive split (edges masked), inductive split
    removes nodes entirely - testing cold-start generalization.
    """
    np.random.seed(seed)
    
    num_enzymes = data['Enzyme'].num_nodes
    n_holdout = int(num_enzymes * holdout_fraction)
    
    # Randomly select holdout enzymes
    all_enzymes = np.arange(num_enzymes)
    np.random.shuffle(all_enzymes)
    holdout_enzymes = set(all_enzymes[:n_holdout])
    train_enzymes = set(all_enzymes[n_holdout:])
    
    return holdout_enzymes, train_enzymes

def remove_nodes_from_graph(data, holdout_enzymes, device):
    """
    Create graph with holdout enzymes completely removed.
    """
    train_data = data.clone()
    num_enzymes = data['Enzyme'].num_nodes
    
    # Create mapping from old to new enzyme indices
    train_enzyme_list = sorted([i for i in range(num_enzymes) if i not in holdout_enzymes])
    old_to_new = {old: new for new, old in enumerate(train_enzyme_list)}
    
    # Filter catalyzes edges
    cat_edge = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    mask = torch.tensor([int(src) not in holdout_enzymes for src in cat_edge[0].cpu().numpy()])
    filtered_cat = cat_edge[:, mask]
    # Remap enzyme indices
    new_src = torch.tensor([old_to_new[int(src)] for src in filtered_cat[0].cpu().numpy()])
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = torch.stack([
        new_src.to(device), filtered_cat[1]
    ])
    
    # Filter rev_catalyzes edges
    rev_edge = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    mask = torch.tensor([int(dst) not in holdout_enzymes for dst in rev_edge[1].cpu().numpy()])
    filtered_rev = rev_edge[:, mask]
    new_dst = torch.tensor([old_to_new[int(dst)] for dst in filtered_rev[1].cpu().numpy()])
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        filtered_rev[0], new_dst.to(device)
    ])
    
    # Update enzyme features
    train_data['Enzyme'].x = data['Enzyme'].x[train_enzyme_list]
    train_data['Enzyme'].num_nodes = len(train_enzyme_list)
    
    return train_data, old_to_new, train_enzyme_list

def train_and_evaluate_inductive(data, seed, device, epochs=20):
    """Train on reduced graph, evaluate on held-out enzymes."""
    set_seed(seed)
    
    # Create inductive split
    holdout_enzymes, train_enzymes = create_inductive_split(data, holdout_fraction=0.1, seed=seed)
    
    # Create training graph (holdout nodes removed)
    train_data, old_to_new, train_enzyme_list = remove_nodes_from_graph(data, holdout_enzymes, device)
    
    # Train model
    model = HGT(train_data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        pos_edge_index = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        
        if pos_edge_index.size(1) == 0:
            continue
        
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
    
    # Evaluate on held-out enzymes
    # For truly inductive: we need to embed held-out enzymes WITHOUT their edges
    # This simulates "new enzyme" scenario
    model.eval()
    
    # Get held-out enzyme edges from original graph
    cat_edge = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    holdout_mask = torch.tensor([int(src) in holdout_enzymes for src in cat_edge[0].cpu().numpy()])
    test_edges = cat_edge[:, holdout_mask]
    
    if test_edges.size(1) == 0:
        return {'hits20': 0.0, 'mrr': 0.0, 'n_holdout': len(holdout_enzymes)}
    
    with torch.no_grad():
        # Get embeddings from training graph
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        met_emb = x_dict['Metabolite']
        
        # For cold-start: use random embeddings for held-out enzymes
        holdout_enz_emb = torch.randn(len(holdout_enzymes), 64).to(device)
        
        test_src = test_edges[0].cpu().numpy()
        test_dst = test_edges[1]
        unique_mets = torch.unique(test_dst)
        
        hits_20 = 0
        reciprocal_ranks = []
        count = 0
        
        holdout_list = sorted(list(holdout_enzymes))
        holdout_to_idx = {old: new for new, old in enumerate(holdout_list)}
        
        for met_idx in unique_mets:
            mask = (test_dst == met_idx)
            true_enzs = [holdout_to_idx[int(e)] for e in test_src[mask.cpu().numpy()] if int(e) in holdout_to_idx]
            
            if len(true_enzs) == 0:
                continue
            
            # Score all holdout enzymes
            num_cands = len(holdout_list)
            eval_src = torch.arange(num_cands, device=device)
            eval_dst = met_idx.repeat(num_cands)
            eval_edges = torch.stack([eval_src, eval_dst])
            
            scores = predictor(holdout_enz_emb, met_emb, eval_edges).sigmoid()
            is_true = torch.zeros(num_cands, dtype=torch.bool, device=device)
            for t in true_enzs:
                is_true[t] = True
            
            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]
            
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1
            
            true_positions = torch.where(sorted_labels)[0]
            if len(true_positions) > 0:
                reciprocal_ranks.append(1.0 / (true_positions[0].item() + 1))
            
            count += 1
    
    hits20_score = hits_20 / count if count > 0 else 0.0
    mrr_score = np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0
    
    return {'hits20': hits20_score, 'mrr': mrr_score, 'n_holdout': len(holdout_enzymes)}

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    seeds = [42, 123, 456]
    results = []
    
    print("Running Inductive (Cold-Start) Evaluation...")
    for seed in seeds:
        result = train_and_evaluate_inductive(data, seed, device)
        results.append(result)
        print(f"Seed {seed}: Hits@20={result['hits20']:.4f}, MRR={result['mrr']:.4f} (n_holdout={result['n_holdout']})")
    
    hits20_arr = np.array([r['hits20'] for r in results])
    mrr_arr = np.array([r['mrr'] for r in results])
    
    print("\n" + "="*50)
    print("INDUCTIVE SPLIT RESULTS")
    print("="*50)
    print(f"Hits@20: {hits20_arr.mean():.4f} ± {hits20_arr.std():.4f}")
    print(f"MRR:     {mrr_arr.mean():.4f} ± {mrr_arr.std():.4f}")
    print("="*50)
    
    # Save
    os.makedirs('results/gnn', exist_ok=True)
    with open('results/gnn/inductive_performance.tsv', 'w') as f:
        f.write("Split_Type\tHits20_Mean\tHits20_Std\tMRR_Mean\tMRR_Std\n")
        f.write(f"inductive\t{hits20_arr.mean():.4f}\t{hits20_arr.std():.4f}\t{mrr_arr.mean():.4f}\t{mrr_arr.std():.4f}\n")
    print("Results saved to results/gnn/inductive_performance.tsv")

if __name__ == "__main__":
    main()
