"""
Minimum Defense Package
1. Permutation test for ET improvement
2. Weight-only / Feature-only ablation
3. GSEA isoflavonoid module
"""

import torch
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, ttest_rel, mannwhitneyu
from scipy.stats import false_discovery_control
from torch_geometric.nn import HGTConv, Linear
import warnings
warnings.filterwarnings('ignore')

ISOFLAVONOID_METS = {'C02495', 'C00858', 'C10216'}

class HGT_Auto(torch.nn.Module):
    def __init__(self, metadata, x_dict, hidden_channels, out_channels, num_heads=4, num_layers=2):
        super().__init__()
        self.lin_dict = torch.nn.ModuleDict()
        for node_type in metadata[0]:
            in_ch = x_dict[node_type].shape[1]
            self.lin_dict[node_type] = Linear(in_ch, hidden_channels)
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            self.convs.append(HGTConv(hidden_channels, hidden_channels, metadata, num_heads))
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
        return self.lin2(self.lin1(torch.cat([src_feat, dst_feat], dim=-1)).relu()).squeeze(-1)

def quick_train_eval(data, device, seed=42):
    """Quick train and eval for ablation."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_mask[indices[:split]] = True
    
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    mask = train_mask[edge_index[0]]
    train_edges = edge_index[:, mask]
    test_edges = edge_index[:, ~mask]
    test_enzymes = indices[split:]
    
    train_data = data.clone()
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, train_mask[rev_index[1]]]
    
    model = HGT_Auto(train_data.metadata(), train_data.x_dict, 64, 64, 4, 2).to(device)
    pred = LinkPredictor(64).to(device)
    opt = torch.optim.Adam(list(model.parameters()) + list(pred.parameters()), lr=0.01)
    
    for _ in range(20):
        model.train()
        opt.zero_grad()
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        pos_e = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        n_met = train_data['Metabolite'].num_nodes
        neg_dst = torch.randint(0, n_met, (pos_e.size(1),), device=device)
        pos_out = pred(x_dict['Enzyme'], x_dict['Metabolite'], pos_e)
        neg_out = pred(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([pos_e[0], neg_dst]))
        loss = -torch.log(torch.sigmoid(pos_out)+1e-15).mean() - torch.log(1-torch.sigmoid(neg_out)+1e-15).mean()
        loss.backward()
        opt.step()
    
    model.eval()
    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        met_list = data['Metabolite'].compound_ids
        
        # Isoflavonoid ranking
        degrees = torch.bincount(train_edges[0], minlength=num_enzymes)
        hub_nodes = torch.argsort(degrees, descending=True)[:20]
        n_met = data['Metabolite'].num_nodes
        scores = torch.zeros(n_met, device=device)
        for enz in hub_nodes:
            scores += (x_dict['Enzyme'][enz].unsqueeze(0) * x_dict['Metabolite']).sum(dim=1)
        ranked = torch.argsort(scores, descending=True).cpu().numpy()
        iso_idx = {i for i, m in enumerate(met_list) if m in ISOFLAVONOID_METS}
        iso_hits = len(set(ranked[:20]) & iso_idx)
        
        # Hits@20
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
    
    return h20/cnt if cnt > 0 else 0, iso_hits

def main():
    print("=" * 60)
    print("Minimum Defense Package")
    print("=" * 60)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Load graphs
    std_graph = torch.load('data/processed/expanded_bipartite_graph.pt').to(device)
    et_graph = torch.load('data/processed/ethylene_conditioned_graph.pt').to(device)
    
    seeds = [42, 123, 456, 789, 999]
    
    # === 1. Paired comparison ===
    print("\n[1] Paired Comparison: Standard vs ET-Conditioned")
    std_results, et_results = [], []
    
    for seed in seeds:
        std_h20, std_iso = quick_train_eval(std_graph, device, seed)
        et_h20, et_iso = quick_train_eval(et_graph, device, seed)
        std_results.append(std_h20)
        et_results.append(et_h20)
        print(f"  Seed {seed}: Std={std_h20:.3f}, ET={et_h20:.3f}, Δ={et_h20-std_h20:+.3f}")
    
    delta = np.array(et_results) - np.array(std_results)
    mean_delta = np.mean(delta)
    ci_95 = 1.96 * np.std(delta) / np.sqrt(len(delta))
    
    print(f"\n  ΔHits@20: {mean_delta:+.3f} [{mean_delta-ci_95:+.3f}, {mean_delta+ci_95:+.3f}] (95% CI)")
    
    # Paired t-test
    t_stat, p_val = ttest_rel(et_results, std_results)
    print(f"  Paired t-test: t={t_stat:.3f}, p={p_val:.4f}")
    
    # === 2. Permutation test ===
    print("\n[2] Permutation Test: Is ET improvement significant?")
    observed_delta = np.mean(et_results) - np.mean(std_results)
    
    n_perm = 100
    perm_deltas = []
    
    for i in range(n_perm):
        # Shuffle ET encoding by permuting enzyme/metabolite indices
        perm_graph = et_graph.clone()
        num_enz = perm_graph['Enzyme'].num_nodes
        perm_idx = torch.randperm(num_enz)
        perm_graph['Enzyme'].x = perm_graph['Enzyme'].x[perm_idx]
        
        perm_h20, _ = quick_train_eval(perm_graph, device, seed=42)
        std_h20, _ = quick_train_eval(std_graph, device, seed=42)
        perm_deltas.append(perm_h20 - std_h20)
        
        if (i+1) % 20 == 0:
            print(f"  Permutation {i+1}/{n_perm}")
    
    perm_p = np.mean(np.array(perm_deltas) >= observed_delta)
    print(f"  Observed Δ = {observed_delta:.3f}")
    print(f"  Permutation p-value = {perm_p:.3f} (n={n_perm})")
    
    # === 3. Ablation: Weight-only vs Feature-only ===
    print("\n[3] Ablation: Weight-only vs Feature-only")
    
    # Weight-only: use ET edge weights but standard 64-dim features
    weight_only_graph = std_graph.clone()
    weight_only_graph['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = \
        et_graph['Enzyme', 'catalyzes', 'Metabolite'].edge_weight[:std_graph['Enzyme', 'catalyzes', 'Metabolite'].edge_index.size(1)]
    
    wo_results = []
    for seed in seeds[:3]:
        h20, _ = quick_train_eval(weight_only_graph, device, seed)
        wo_results.append(h20)
    print(f"  Weight-only: {np.mean(wo_results):.3f} ± {np.std(wo_results):.3f}")
    
    # Feature-only: use ET features but standard edge weights
    feat_only_graph = et_graph.clone()
    # Reset edge weights to 1.0
    n_edges = feat_only_graph['Enzyme', 'catalyzes', 'Metabolite'].edge_index.size(1)
    feat_only_graph['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = torch.ones(n_edges)
    
    fo_results = []
    for seed in seeds[:3]:
        h20, _ = quick_train_eval(feat_only_graph, device, seed)
        fo_results.append(h20)
    print(f"  Feature-only: {np.mean(fo_results):.3f} ± {np.std(fo_results):.3f}")
    
    # === 4. GSEA-like module analysis ===
    print("\n[4] GSEA-like Module Enrichment")
    
    # Use factorial analysis results
    factorial_df = pd.read_csv('results/ethylene_main_effect.csv')
    
    print("  Isoflavonoid targets from factorial analysis:")
    for _, row in factorial_df.iterrows():
        if row['kegg_id'] in ['C02495', 'C00858', 'C10216']:
            print(f"    {row['name']}: q={row['q_ethylene']:.2e}")
    
    # Module-level summary
    iso_targets = factorial_df[factorial_df['kegg_id'].isin(['C02495', 'C00858', 'C10216'])]
    mean_log2fc = iso_targets['log2fc_ethylene'].mean()
    
    # Fisher's method for combining p-values
    from scipy.stats import combine_pvalues
    pvals = iso_targets['p_ethylene'].values
    combined_stat, combined_p = combine_pvalues(pvals, method='fisher')
    
    print(f"\n  Isoflavonoid Module (n=3):")
    print(f"    Mean Log2FC: {mean_log2fc:.3f}")
    print(f"    Combined p-value (Fisher): {combined_p:.2e}")
    
    # === Summary ===
    print("\n" + "=" * 60)
    print("SUMMARY FOR PAPER")
    print("=" * 60)
    
    results_summary = {
        'delta_hits20': f"{mean_delta:+.3f} [{mean_delta-ci_95:+.3f}, {mean_delta+ci_95:+.3f}]",
        'paired_p': p_val,
        'permutation_p': perm_p,
        'weight_only': f"{np.mean(wo_results):.3f}",
        'feature_only': f"{np.mean(fo_results):.3f}",
        'iso_module_p': combined_p
    }
    
    for k, v in results_summary.items():
        print(f"  {k}: {v}")
    
    # Save
    pd.DataFrame([results_summary]).to_csv('results/defense_package.csv', index=False)
    print("\nSaved to: results/defense_package.csv")

if __name__ == "__main__":
    main()
