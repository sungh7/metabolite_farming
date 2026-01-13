"""
Train and evaluate HGT on ethylene-conditioned graph.
Tests whether ethylene encoding improves isoflavonoid module prioritization.
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import hypergeom
import sys
import os

sys.path.append(os.getcwd())

# Isoflavonoid metabolites
ISOFLAVONOID_METS = {'C02495', 'C00858', 'C10216'}

class HGT_ET(torch.nn.Module):
    """HGT adapted for ethylene-conditioned features (dim 72)."""
    def __init__(self, metadata, in_channels, hidden_channels, out_channels, num_heads=4, num_layers=2):
        super().__init__()
        from torch_geometric.nn import HGTConv, Linear
        
        self.lin_dict = torch.nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = Linear(in_channels, hidden_channels)
        
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            conv = HGTConv(hidden_channels, hidden_channels, metadata, num_heads)
            self.convs.append(conv)
        
        self.lin_out = Linear(hidden_channels, out_channels)
    
    def forward(self, x_dict, edge_index_dict):
        for node_type, x in x_dict.items():
            x_dict[node_type] = self.lin_dict[node_type](x).relu()
        
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
        
        for node_type, x in x_dict.items():
            x_dict[node_type] = self.lin_out(x)
        
        return x_dict

class LinkPredictor(torch.nn.Module):
    def __init__(self, in_channels):
        super().__init__()
        self.lin1 = torch.nn.Linear(in_channels * 2, in_channels)
        self.lin2 = torch.nn.Linear(in_channels, 1)
    
    def forward(self, x_src, x_dst, edge_index):
        src_feat = x_src[edge_index[0]]
        dst_feat = x_dst[edge_index[1]]
        x = torch.cat([src_feat, dst_feat], dim=-1)
        x = self.lin1(x).relu()
        return self.lin2(x).squeeze(-1)

def train_and_evaluate(data, device, seed=42):
    """Train HGT and evaluate isoflavonoid enrichment."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    num_enzymes = data['Enzyme'].num_nodes
    in_dim = data['Enzyme'].x.shape[1]  # 72 for ethylene-conditioned
    
    # Split
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    edge_weight = data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight
    
    mask = train_enz_mask[edge_index[0]]
    train_edges = edge_index[:, mask]
    train_weights = edge_weight[mask]
    test_edges = edge_index[:, ~mask]
    test_enzymes = indices[split:]
    
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = train_weights
    
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_mask = train_enz_mask[rev_index[1]]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    # Model
    model = HGT_ET(train_data.metadata(), in_dim, 64, 64, num_heads=4, num_layers=2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    # Train with weighted loss
    for epoch in range(20):
        model.train()
        optimizer.zero_grad()
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        pos_edges = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        pos_weights = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight
        
        num_pos = pos_edges.size(1)
        num_met = train_data['Metabolite'].num_nodes
        neg_dst = torch.randint(0, num_met, (num_pos,), device=device)
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edges)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], 
                           torch.stack([pos_edges[0], neg_dst]))
        
        # Weighted BCE - higher weight edges contribute more
        pos_loss = -(pos_weights * torch.log(torch.sigmoid(pos_out) + 1e-15)).mean()
        neg_loss = -torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss = pos_loss + neg_loss
        
        loss.backward()
        optimizer.step()
    
    # Evaluate: source-conditioned ranking
    model.eval()
    met_list = data['Metabolite'].compound_ids
    
    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        # Use high-degree enzymes as source (proxy for signaling hubs)
        degrees = torch.bincount(train_edges[0], minlength=num_enzymes)
        hub_nodes = torch.argsort(degrees, descending=True)[:20]
        
        # Compute metabolite scores from hubs
        num_metabolites = data['Metabolite'].num_nodes
        scores = torch.zeros(num_metabolites, device=device)
        
        for enz_idx in hub_nodes:
            enz_feat = x_dict['Enzyme'][enz_idx].unsqueeze(0).expand(num_metabolites, -1)
            score = (enz_feat * x_dict['Metabolite']).sum(dim=1)
            scores += score
        
        ranked_indices = torch.argsort(scores, descending=True).cpu().numpy()
        
        # Isoflavonoid enrichment
        iso_indices = set()
        for i, met_id in enumerate(met_list):
            if met_id in ISOFLAVONOID_METS:
                iso_indices.add(i)
        
        top_20_set = set(ranked_indices[:20])
        hits = len(top_20_set & iso_indices)
        
        pvalue = 1 - hypergeom.cdf(hits - 1, num_metabolites, len(iso_indices), 20)
        
        # Also compute standard Hits@20
        unique_mets = torch.unique(test_edges[1])
        hits_20, count = 0, 0
        
        for met_idx in unique_mets:
            mask = (test_edges[1] == met_idx)
            true_enzs = test_edges[0][mask]
            
            num_cands = test_enzymes.size(0)
            eval_edges = torch.stack([test_enzymes, met_idx.repeat(num_cands)])
            scores = predictor(x_dict['Enzyme'], x_dict['Metabolite'], eval_edges).sigmoid()
            
            is_true = torch.isin(test_enzymes, true_enzs)
            if is_true.sum() == 0:
                continue
            
            sorted_idx = torch.argsort(scores, descending=True)
            if is_true[sorted_idx][:20].sum() > 0:
                hits_20 += 1
            count += 1
        
        hits_20_rate = hits_20 / count if count > 0 else 0
    
    return {
        'iso_hits': hits,
        'iso_pvalue': pvalue,
        'hits_20': hits_20_rate,
        'top_10_mets': [met_list[i] for i in ranked_indices[:10]]
    }

def main():
    print("=" * 60)
    print("Ethylene-Conditioned Graph Evaluation")
    print("=" * 60)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    results = []
    
    # Compare: Standard graph vs Ethylene-conditioned graph
    for graph_name, graph_path in [
        ('Standard', 'data/processed/expanded_bipartite_graph.pt'),
        ('ET-Conditioned', 'data/processed/ethylene_conditioned_graph.pt')
    ]:
        print(f"\n--- {graph_name} Graph ---")
        
        data = torch.load(graph_path).to(device)
        print(f"Feature dim: {data['Enzyme'].x.shape[1]}")
        
        # Run 3 seeds
        iso_hits_list, hits_20_list, pval_list = [], [], []
        
        for seed in [42, 123, 456]:
            res = train_and_evaluate(data, device, seed)
            iso_hits_list.append(res['iso_hits'])
            hits_20_list.append(res['hits_20'])
            pval_list.append(res['iso_pvalue'])
            print(f"  Seed {seed}: Iso={res['iso_hits']}/3 (p={res['iso_pvalue']:.4f}), "
                  f"Hits@20={res['hits_20']:.3f}")
        
        results.append({
            'graph': graph_name,
            'iso_hits_mean': np.mean(iso_hits_list),
            'iso_pvalue_mean': np.mean(pval_list),
            'hits_20_mean': np.mean(hits_20_list),
            'hits_20_std': np.std(hits_20_list)
        })
        
        print(f"  Mean: Iso={np.mean(iso_hits_list):.1f}, Hits@20={np.mean(hits_20_list):.3f}")
    
    # Save
    df = pd.DataFrame(results)
    df.to_csv('results/ethylene_graph_comparison.csv', index=False)
    
    print("\n" + "=" * 60)
    print("Comparison Summary")
    print("=" * 60)
    print(df.to_string(index=False))
    
    # Interpretation
    std_res = df[df['graph'] == 'Standard'].iloc[0]
    et_res = df[df['graph'] == 'ET-Conditioned'].iloc[0]
    
    if et_res['iso_hits_mean'] > std_res['iso_hits_mean']:
        print(f"\n✓ Ethylene conditioning IMPROVES isoflavonoid enrichment!")
        print(f"  Standard: {std_res['iso_hits_mean']:.1f}/3")
        print(f"  ET-Conditioned: {et_res['iso_hits_mean']:.1f}/3")
    else:
        print(f"\n× No improvement from ethylene conditioning")
    
    print(f"\nSaved to: results/ethylene_graph_comparison.csv")

if __name__ == "__main__":
    main()
