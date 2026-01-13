"""
Tier-P Ablation Experiment
Compares performance across different edge configurations:
1. Tier-R only (reaction-grounded)
2. Tier-R + Tier-P (current)
3. Tier-P only (pathway-supported)
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import average_precision_score
import sys
import os

sys.path.append(os.getcwd())
from src.model import HGT, LinkPredictor

def train_and_evaluate(data, device, seed=42):
    """Train HGT and return Hits@20."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_enz_mask[src]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1]
    rev_mask = train_enz_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
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
    
    # Evaluate
    model.eval()
    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        test_enzymes = indices[split:]
        
        unique_mets = torch.unique(test_edges[1])
        hits_20 = 0
        count = 0
        
        for met_idx in unique_mets:
            mask = (test_edges[1] == met_idx)
            true_enzs = test_edges[0][mask]
            
            num_cands = test_enzymes.size(0)
            eval_edges = torch.stack([test_enzymes, met_idx.repeat(num_cands)])
            scores = predictor(x_dict['Enzyme'], x_dict['Metabolite'], eval_edges).sigmoid()
            
            is_true = torch.isin(test_enzymes, true_enzs)
            if is_true.sum() == 0:
                continue
            
            sorted_indices = torch.argsort(scores, descending=True)
            sorted_labels = is_true[sorted_indices]
            
            if sorted_labels[:20].sum() > 0:
                hits_20 += 1
            count += 1
    
    return hits_20 / count if count > 0 else 0, train_edges.size(1), test_edges.size(1)

def build_tier_graph(base_graph_path, full_edges_path, tier_config):
    """Build graph with specified tier configuration."""
    data = torch.load(base_graph_path).clone()
    full_edges = pd.read_csv(full_edges_path, sep='\t')
    
    # Load metabolite list from expanded graph
    expanded = torch.load('data/processed/expanded_bipartite_graph.pt')
    met_list = expanded['Metabolite'].compound_ids
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    
    num_enzymes = data['Enzyme'].num_nodes
    n_metabolites = len(met_list)
    
    # EC to enzyme mapping (from original expanded builder)
    np.random.seed(42)
    ec_to_enz = {}
    for ec in full_edges['enzyme_ec'].unique():
        n_enz = np.random.randint(5, 20)
        ec_to_enz[ec] = np.random.choice(num_enzymes, n_enz, replace=False).tolist()
    
    # MTBLS531 covered metabolites (Tier-R)
    covered_mtbls = {'C00062', 'C00078', 'C00858', 'C02495', 'C10216', 'C01177', 'C06037', 'C04079', 'C19865', 'C01004'}
    
    tier_r_src, tier_r_dst = [], []
    tier_p_src, tier_p_dst = [], []
    
    for _, row in full_edges.iterrows():
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']
        
        if met_id not in met_to_idx:
            continue
        
        met_idx = met_to_idx[met_id]
        enz_indices = ec_to_enz.get(ec, [])
        is_mtbls = met_id in covered_mtbls
        
        for enz_idx in enz_indices:
            if is_mtbls:
                tier_r_src.append(enz_idx)
                tier_r_dst.append(met_idx)
            else:
                tier_p_src.append(enz_idx)
                tier_p_dst.append(met_idx)
    
    # Apply tier configuration
    if tier_config == 'tier_r_only':
        all_src, all_dst = tier_r_src, tier_r_dst
        weights = [1.0] * len(tier_r_src)
    elif tier_config == 'tier_p_only':
        all_src, all_dst = tier_p_src, tier_p_dst
        weights = [0.5] * len(tier_p_src)
    else:  # tier_r_plus_p
        all_src = tier_r_src + tier_p_src
        all_dst = tier_r_dst + tier_p_dst
        weights = [1.0] * len(tier_r_src) + [0.5] * len(tier_p_src)
    
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = torch.randn(n_metabolites, 64)
    
    edge_index = torch.tensor([all_src, all_dst], dtype=torch.long)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = torch.tensor(weights, dtype=torch.float)
    
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(all_dst), torch.tensor(all_src)
    ])
    
    return data, len(tier_r_src), len(tier_p_src)

def main():
    print("=" * 60)
    print("Tier-P Ablation Experiment")
    print("=" * 60)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    base_graph = 'data/processed/strict_graph.pt'
    full_edges = 'data/kegg/full_enzyme_metabolite_edges.tsv'
    
    results = []
    
    for tier_config in ['tier_r_only', 'tier_r_plus_p', 'tier_p_only']:
        print(f"\n--- {tier_config} ---")
        
        data, n_tier_r, n_tier_p = build_tier_graph(base_graph, full_edges, tier_config)
        data = data.to(device)
        
        total_edges = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index.size(1)
        print(f"Edges: {total_edges} (Tier-R: {n_tier_r if tier_config != 'tier_p_only' else 0}, Tier-P: {n_tier_p if tier_config != 'tier_r_only' else 0})")
        
        # Run 3 seeds
        hits_list = []
        for seed in [42, 123, 456]:
            hits, train_n, test_n = train_and_evaluate(data, device, seed)
            hits_list.append(hits)
            print(f"  Seed {seed}: Hits@20 = {hits:.4f}")
        
        mean_hits = np.mean(hits_list)
        std_hits = np.std(hits_list)
        
        results.append({
            'config': tier_config,
            'edges': total_edges,
            'tier_r': n_tier_r if tier_config != 'tier_p_only' else 0,
            'tier_p': n_tier_p if tier_config != 'tier_r_only' else 0,
            'hits_20_mean': mean_hits,
            'hits_20_std': std_hits
        })
        
        print(f"  Mean: {mean_hits:.4f} ± {std_hits:.4f}")
    
    # Save results
    df = pd.DataFrame(results)
    df.to_csv('results/tier_ablation.csv', index=False)
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(df.to_string(index=False))
    print(f"\nSaved to: results/tier_ablation.csv")

if __name__ == "__main__":
    main()
