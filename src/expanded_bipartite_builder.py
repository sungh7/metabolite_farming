"""
Expanded Tiered Bipartite Graph Builder
Uses full KEGG data (5,302 edges) instead of sample.
"""

import torch
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
import os
from pathlib import Path

CURRENCY_METABOLITES = {
    'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006',
    'C00008', 'C00009', 'C00010', 'C00011', 'C00013', 'C00014',
    'C00020', 'C00027', 'C00044', 'C00080'
}

def build_expanded_graph(graph_path: str, output_path: str, kegg_dir: str = "data/kegg"):
    """Build graph using full KEGG edge data."""
    print("=" * 60)
    print("Building Expanded Bipartite Graph")
    print("=" * 60)
    
    kegg_dir = Path(kegg_dir)
    
    # 1. Load Base PPI Graph
    base_graph = torch.load(graph_path)
    data = base_graph.clone()
    num_enzymes = data['Enzyme'].num_nodes
    print(f"Base graph: {num_enzymes} enzymes")
    
    # 2. Load full KEGG edges
    edges_path = kegg_dir / "full_enzyme_metabolite_edges.tsv"
    if not edges_path.exists():
        print("Error: full_enzyme_metabolite_edges.tsv not found")
        return
    
    edges_df = pd.read_csv(edges_path, sep='\t')
    print(f"Full KEGG edges: {len(edges_df)}")
    
    # 3. Load MTBLS531 metabolites
    mtbls_path = kegg_dir / "metabolites.csv"
    mtbls_df = pd.read_csv(mtbls_path)
    target_mets = set(mtbls_df['compound_id'])
    
    # 4. Get all unique metabolites from KEGG edges
    all_kegg_mets = set(edges_df['metabolite_id'].unique())
    
    # Filter out currency metabolites
    all_kegg_mets = all_kegg_mets - CURRENCY_METABOLITES
    
    # Covered MTBLS531 metabolites
    covered_mtbls = target_mets & all_kegg_mets
    print(f"MTBLS531 coverage: {len(covered_mtbls)}/{len(target_mets)}")
    print(f"Covered: {sorted(covered_mtbls)}")
    
    # 5. Create metabolite index (only use covered MTBLS531 + some extra for network)
    # Strategy: Include all MTBLS531 mets + sample from KEGG for structure
    
    # All MTBLS531 mets as primary
    met_list = sorted(covered_mtbls)
    
    # Add more from KEGG to enrich network (up to 200 total)
    uncovered_mtbls = target_mets - covered_mtbls
    remaining_kegg = list(all_kegg_mets - covered_mtbls)[:200 - len(met_list)]
    met_list.extend(remaining_kegg)
    
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    n_metabolites = len(met_list)
    print(f"Total metabolites in graph: {n_metabolites}")
    
    # 6. Load gene-EC mapping to connect ECs to enzyme indices
    gene_ec_path = kegg_dir / "gene_ec_mapping.tsv"
    gene_ec_df = pd.read_csv(gene_ec_path, sep='\t')
    
    # Create EC to enzyme indices mapping (simplified: random assignment)
    # In production: map gene_id to STRING ID
    ec_to_enz = {}
    np.random.seed(42)
    
    for ec in edges_df['enzyme_ec'].unique():
        # Assign 5-20 random enzymes to each EC
        n_enz = np.random.randint(5, 20)
        ec_to_enz[ec] = np.random.choice(num_enzymes, n_enz, replace=False).tolist()
    
    # 7. Build edges
    tier_r_src, tier_r_dst = [], []  # Reaction-grounded for MTBLS531 mets
    tier_p_src, tier_p_dst = [], []  # Others
    
    for _, row in edges_df.iterrows():
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
    
    print(f"Tier-R edges (MTBLS531 mets): {len(tier_r_src)}")
    print(f"Tier-P edges (KEGG network): {len(tier_p_src)}")
    
    # 8. Combine with weights
    all_src = tier_r_src + tier_p_src
    all_dst = tier_r_dst + tier_p_dst
    all_weights = [1.0] * len(tier_r_src) + [0.5] * len(tier_p_src)
    
    # 9. Update HeteroData
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = torch.randn(n_metabolites, 64)
    data['Metabolite'].compound_ids = met_list
    
    edge_index = torch.tensor([all_src, all_dst], dtype=torch.long)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = torch.tensor(all_weights, dtype=torch.float)
    
    # Reverse edges
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(all_dst), torch.tensor(all_src)
    ])
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_weight = torch.tensor(all_weights, dtype=torch.float)
    
    # 10. Summary
    print("\n" + "=" * 60)
    print("Graph Summary")
    print("=" * 60)
    print(f"Nodes: Enzymes={num_enzymes}, Metabolites={n_metabolites}")
    print(f"Edges: Total={len(all_src)}, Tier-R={len(tier_r_src)}, Tier-P={len(tier_p_src)}")
    print(f"Edge weight mean: {np.mean(all_weights):.3f}")
    
    # Check core metabolites
    core_mets = {'C00062', 'C00858', 'C02495'}
    covered_core = core_mets & set(met_list)
    print(f"Core metabolites: {len(covered_core)}/3")
    
    # 11. Save
    torch.save(data, output_path)
    print(f"\nSaved to: {output_path}")
    
    return data

if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    
    build_expanded_graph(
        'data/processed/strict_graph.pt',
        'data/processed/expanded_bipartite_graph.pt'
    )
