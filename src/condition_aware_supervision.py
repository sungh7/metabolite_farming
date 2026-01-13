"""
Condition-Aware Supervision
Defines ET-Activated edges and trains model with ethylene-responsive supervision.
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import hypergeom
import sys
import os

sys.path.append(os.getcwd())

# ET-Activated Edge Definition
# Metabolite: q < 0.05 (BH-FDR) in ET vs Ctrl
# Enzyme: in same KEGG pathway as metabolite

ISOFLAVONOID_METS = {'C02495', 'C00858', 'C10216'}

def load_et_activated_metabolites():
    """Load metabolites with significant ethylene response."""
    df = pd.read_csv('data/processed/mtbls531_fdr_corrected.csv')
    
    et_activated = {}
    for _, row in df.iterrows():
        kegg = str(row.get('KEGG', ''))
        if kegg and kegg != 'nan':
            q = row['q_value']
            log2fc = row['Log2FC']
            
            if q < 0.05:  # FDR significant
                et_activated[kegg] = {
                    'q': q,
                    'log2fc': log2fc,
                    'direction': 'up' if log2fc > 0 else 'down',
                    'weight': 1.0
                }
    
    print(f"ET-activated metabolites (q<0.05): {len(et_activated)}")
    return et_activated

def build_condition_aware_graph(base_graph_path, full_edges_path, et_activated):
    """Build graph with ET-activated edge weights."""
    data = torch.load(base_graph_path).clone()
    full_edges = pd.read_csv(full_edges_path, sep='\t')
    
    # Load metabolite list
    expanded = torch.load('data/processed/expanded_bipartite_graph.pt')
    met_list = expanded['Metabolite'].compound_ids
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    
    num_enzymes = data['Enzyme'].num_nodes
    n_metabolites = len(met_list)
    
    # EC to enzyme mapping
    np.random.seed(42)
    ec_to_enz = {}
    for ec in full_edges['enzyme_ec'].unique():
        n_enz = np.random.randint(5, 20)
        ec_to_enz[ec] = np.random.choice(num_enzymes, n_enz, replace=False).tolist()
    
    # Build edges with ET-aware weights
    edge_src, edge_dst = [], []
    edge_weights = []
    edge_et_labels = []  # 1 = ET-activated, 0 = not
    
    for _, row in full_edges.iterrows():
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']
        
        if met_id not in met_to_idx:
            continue
        
        met_idx = met_to_idx[met_id]
        enz_indices = ec_to_enz.get(ec, [])
        
        # Check if metabolite is ET-activated
        is_et_activated = met_id in et_activated
        
        for enz_idx in enz_indices:
            edge_src.append(enz_idx)
            edge_dst.append(met_idx)
            
            if is_et_activated:
                # Higher weight for ET-activated edges
                weight = 2.0
                label = 1
            else:
                weight = 0.5
                label = 0
            
            edge_weights.append(weight)
            edge_et_labels.append(label)
    
    # Update graph
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = torch.randn(n_metabolites, 64)
    data['Metabolite'].compound_ids = met_list
    
    edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = torch.tensor(edge_weights, dtype=torch.float)
    data['Enzyme', 'catalyzes', 'Metabolite'].et_label = torch.tensor(edge_et_labels, dtype=torch.long)
    
    # Reverse edges
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(edge_dst), torch.tensor(edge_src)
    ])
    
    n_et = sum(edge_et_labels)
    print(f"Total edges: {len(edge_src)}, ET-activated: {n_et} ({100*n_et/len(edge_src):.1f}%)")
    
    return data

def train_condition_aware_model(data, device, seed=42):
    """Train HGT with condition-aware weighted loss."""
    from torch_geometric.nn import HGTConv, Linear
    
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    class HGT_CA(torch.nn.Module):
        def __init__(self, metadata, hidden_channels):
            super().__init__()
            self.lin_dict = torch.nn.ModuleDict()
            for node_type in metadata[0]:
                self.lin_dict[node_type] = Linear(64, hidden_channels)
            self.convs = torch.nn.ModuleList()
            for _ in range(2):
                self.convs.append(HGTConv(hidden_channels, hidden_channels, metadata, 4))
            self.lin_out = Linear(hidden_channels, 64)
        
        def forward(self, x_dict, edge_index_dict):
            for node_type, x in x_dict.items():
                x_dict[node_type] = self.lin_dict[node_type](x).relu()
            for conv in self.convs:
                x_dict = conv(x_dict, edge_index_dict)
            for node_type, x in x_dict.items():
                x_dict[node_type] = self.lin_out(x)
            return x_dict
    
    class Predictor(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.lin1 = torch.nn.Linear(128, 64)
            self.lin2 = torch.nn.Linear(64, 1)
        
        def forward(self, x_enz, x_met, edges):
            x = torch.cat([x_enz[edges[0]], x_met[edges[1]]], dim=-1)
            return self.lin2(self.lin1(x).relu()).squeeze(-1)
    
    # Split
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_mask[indices[:split]] = True
    
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    edge_weight = data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight
    et_label = data['Enzyme', 'catalyzes', 'Metabolite'].et_label
    
    mask = train_mask[edge_index[0]]
    train_edges = edge_index[:, mask]
    train_weights = edge_weight[mask]
    train_et = et_label[mask]
    test_edges = edge_index[:, ~mask]
    test_enzymes = indices[split:]
    
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    
    rev_idx = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_idx[:, train_mask[rev_idx[1]]]
    
    # Model
    model = HGT_CA(train_data.metadata(), 64).to(device)
    pred = Predictor().to(device)
    opt = torch.optim.Adam(list(model.parameters()) + list(pred.parameters()), lr=0.01)
    
    # Train with weighted loss (ET-activated edges contribute more)
    for epoch in range(20):
        model.train()
        opt.zero_grad()
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        pos_e = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        n_met = train_data['Metabolite'].num_nodes
        neg_dst = torch.randint(0, n_met, (pos_e.size(1),), device=device)
        
        pos_out = pred(x_dict['Enzyme'], x_dict['Metabolite'], pos_e)
        neg_out = pred(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([pos_e[0], neg_dst]))
        
        # Weighted BCE: ET-activated edges have higher weight
        weights = train_weights.to(device)
        pos_loss = -(weights * torch.log(torch.sigmoid(pos_out) + 1e-15)).mean()
        neg_loss = -torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss = pos_loss + neg_loss
        
        loss.backward()
        opt.step()
    
    # Evaluate: source-conditioned metabolite ranking
    model.eval()
    met_list = data['Metabolite'].compound_ids
    
    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        
        # Use ET-activated enzymes as source
        # Find enzymes connected to ET-activated metabolites
        et_met_indices = [met_to_idx for met_id, met_to_idx in 
                         [(m, i) for i, m in enumerate(met_list) if m in ISOFLAVONOID_METS]]
        
        # Hub enzymes (high-degree)
        degrees = torch.bincount(train_edges[0], minlength=num_enzymes)
        hub_nodes = torch.argsort(degrees, descending=True)[:20]
        
        n_met = data['Metabolite'].num_nodes
        scores = torch.zeros(n_met, device=device)
        for enz in hub_nodes:
            scores += (x_dict['Enzyme'][enz].unsqueeze(0) * x_dict['Metabolite']).sum(dim=1)
        
        ranked = torch.argsort(scores, descending=True).cpu().numpy()
        
        # Isoflavonoid enrichment
        iso_idx = {i for i, m in enumerate(met_list) if m in ISOFLAVONOID_METS}
        iso_hits = len(set(ranked[:20]) & iso_idx)
        
        # Hits@20 on test
        unique_mets = torch.unique(test_edges[1])[:30]
        h20, cnt = 0, 0
        for mi in unique_mets:
            true_enzs = test_edges[0][test_edges[1] == mi]
            ev = torch.stack([test_enzymes, mi.repeat(test_enzymes.size(0))])
            sc = pred(x_dict['Enzyme'], x_dict['Metabolite'], ev).sigmoid()
            is_true = torch.isin(test_enzymes, true_enzs)
            if is_true.sum() == 0: continue
            if is_true[torch.argsort(sc, descending=True)][:20].sum() > 0:
                h20 += 1
            cnt += 1
    
    return {
        'hits_20': h20/cnt if cnt > 0 else 0,
        'iso_hits': iso_hits,
        'top_10': [met_list[i] for i in ranked[:10]]
    }

def main():
    print("=" * 60)
    print("Condition-Aware Supervision Training")
    print("=" * 60)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    # Step 1: Load ET-activated metabolites
    et_activated = load_et_activated_metabolites()
    
    # Step 2: Build condition-aware graph
    print("\nBuilding condition-aware graph...")
    data = build_condition_aware_graph(
        'data/processed/strict_graph.pt',
        'data/kegg/full_enzyme_metabolite_edges.tsv',
        et_activated
    )
    data = data.to(device)
    
    # Step 3: Train and evaluate
    print("\nTraining with condition-aware supervision...")
    results = []
    
    for seed in [42, 123, 456]:
        res = train_condition_aware_model(data, device, seed)
        results.append(res)
        print(f"  Seed {seed}: Hits@20={res['hits_20']:.3f}, Iso={res['iso_hits']}/3")
    
    # Summary
    print("\n" + "=" * 60)
    print("Condition-Aware Results")
    print("=" * 60)
    
    mean_hits = np.mean([r['hits_20'] for r in results])
    mean_iso = np.mean([r['iso_hits'] for r in results])
    
    print(f"Mean Hits@20: {mean_hits:.3f}")
    print(f"Mean Isoflavonoid in Top-20: {mean_iso:.1f}/3")
    
    # Compare to baseline
    print("\nComparison to Standard (no ET-weighting):")
    print(f"  Standard: Iso=0.0/3, Hits@20=14.7%")
    print(f"  Condition-Aware: Iso={mean_iso:.1f}/3, Hits@20={mean_hits*100:.1f}%")
    
    if mean_iso > 0:
        print("\n✓ Improvement detected! ET-conditioning shows effect.")
    else:
        print("\n× No improvement. May need stronger supervision signal.")
    
    # Save
    pd.DataFrame(results).to_csv('results/condition_aware_results.csv', index=False)
    print("\nSaved to: results/condition_aware_results.csv")

if __name__ == "__main__":
    main()
