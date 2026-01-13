"""
All-Constant Features Ablation

Tests if node identity signal is necessary by replacing all node features
with identical constant vectors. If performance drops significantly,
it confirms that node-distinguishing information is essential.
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

def train_and_evaluate(data, train_mask, test_mask, seed, device, 
                       feature_mode='learnable', epochs=20):
    """
    Train and evaluate with different feature modes.
    
    feature_mode: 'learnable', 'constant', 'shuffled'
    """
    set_seed(seed)
    
    train_data = data.clone()
    num_enzymes = data['Enzyme'].num_nodes
    
    # Handle feature modes
    if feature_mode == 'constant':
        # All nodes get identical constant vector
        for node_type in train_data.node_types:
            if hasattr(train_data[node_type], 'x') and train_data[node_type].x is not None:
                n_nodes = train_data[node_type].num_nodes
                dim = train_data[node_type].x.size(1)
                # Same vector for all nodes (no identity signal)
                train_data[node_type].x = torch.ones(n_nodes, dim, device=device)
    
    elif feature_mode == 'shuffled':
        # Shuffle features across nodes (break ID-feature correspondence)
        for node_type in train_data.node_types:
            if hasattr(train_data[node_type], 'x') and train_data[node_type].x is not None:
                n_nodes = train_data[node_type].num_nodes
                perm = torch.randperm(n_nodes)
                train_data[node_type].x = train_data[node_type].x[perm]
    
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
    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()), lr=0.01
    )
    
    # Training
    for epoch in range(epochs):
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
        reciprocal_ranks = []
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
            
            true_positions = torch.where(sorted_labels)[0]
            if len(true_positions) > 0:
                reciprocal_ranks.append(1.0 / (true_positions[0].item() + 1))
            
            count += 1
    
    hits20_score = hits_20 / count if count > 0 else 0.0
    mrr_score = np.mean(reciprocal_ranks) if reciprocal_ranks else 0.0
    
    return {'hits20': hits20_score, 'mrr': mrr_score}

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    num_enzymes = data['Enzyme'].num_nodes
    seeds = [42, 123, 456]
    n_folds = 5
    
    modes = ['learnable', 'constant', 'shuffled']
    all_results = {mode: {'hits20': [], 'mrr': []} for mode in modes}
    
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    enzyme_indices = np.arange(num_enzymes)
    
    print("Running Feature Mode Ablation...")
    
    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(enzyme_indices)):
        train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
        train_mask[train_idx] = True
        test_mask = ~train_mask
        
        for seed in seeds:
            for mode in modes:
                result = train_and_evaluate(
                    data, train_mask, test_mask, seed, device,
                    feature_mode=mode
                )
                all_results[mode]['hits20'].append(result['hits20'])
                all_results[mode]['mrr'].append(result['mrr'])
                print(f"Fold {fold_idx+1}, Seed {seed}, {mode}: Hits@20={result['hits20']:.4f}")
    
    print("\n" + "="*60)
    print("FEATURE MODE ABLATION RESULTS")
    print("="*60)
    
    for mode in modes:
        hits20 = np.array(all_results[mode]['hits20'])
        mrr = np.array(all_results[mode]['mrr'])
        print(f"\n{mode.upper()}:")
        print(f"  Hits@20: {hits20.mean():.4f} ± {hits20.std():.4f}")
        print(f"  MRR:     {mrr.mean():.4f} ± {mrr.std():.4f}")
    
    # Save results
    os.makedirs('results/gnn', exist_ok=True)
    
    with open('results/gnn/constant_feature_ablation.tsv', 'w') as f:
        f.write("Mode\tHits20_Mean\tHits20_Std\tMRR_Mean\tMRR_Std\n")
        for mode in modes:
            hits20 = np.array(all_results[mode]['hits20'])
            mrr = np.array(all_results[mode]['mrr'])
            f.write(f"{mode}\t{hits20.mean():.4f}\t{hits20.std():.4f}\t{mrr.mean():.4f}\t{mrr.std():.4f}\n")
    
    print("\nResults saved to results/gnn/constant_feature_ablation.tsv")

if __name__ == "__main__":
    main()
