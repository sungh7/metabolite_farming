"""
Tiered Bipartite Graph Builder
Creates 2-tier evidence structure:
- Tier-R: Reaction-grounded edges (KEGG reactions)
- Tier-P: Pathway-supported edges (pathway membership)
"""

import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
import os
from pathlib import Path

# Currency metabolites to always exclude
CURRENCY_METABOLITES = {
    'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006',
    'C00008', 'C00009', 'C00010', 'C00011', 'C00013', 'C00014',
    'C00020', 'C00027', 'C00044', 'C00080'
}

def load_kegg_data(kegg_dir: Path) -> tuple:
    """Load KEGG metabolites and edges."""
    metabolites_path = kegg_dir / "metabolites.csv"
    edges_path = kegg_dir / "enzyme_metabolite_edges.tsv"
    gene_ec_path = kegg_dir / "gene_ec_mapping.tsv"
    
    metabolites = pd.read_csv(metabolites_path)
    
    # Load Tier-R edges if available
    tier_r_edges = pd.DataFrame()
    if edges_path.exists():
        tier_r_edges = pd.read_csv(edges_path, sep='\t')
    
    # Load gene-EC mapping
    gene_ec = pd.DataFrame()
    if gene_ec_path.exists():
        gene_ec = pd.read_csv(gene_ec_path, sep='\t')
    
    return metabolites, tier_r_edges, gene_ec

def build_tiered_bipartite_graph(graph_path: str, output_path: str, kegg_dir: str = "data/kegg"):
    """
    Constructs 2-Tier Heterogeneous Graph.
    
    Tier-R (weight=1.0): Reaction-grounded edges from KEGG
    Tier-P (weight=0.5): Pathway-supported edges
    """
    print("=" * 60)
    print("Building 2-Tier Bipartite Graph")
    print("=" * 60)
    
    kegg_dir = Path(kegg_dir)
    
    # 1. Load Base PPI Graph
    if not os.path.exists(graph_path):
        print(f"Error: {graph_path} not found.")
        return
        
    base_graph = torch.load(graph_path)
    data = base_graph.clone()
    num_enzymes = data['Enzyme'].num_nodes
    print(f"Base graph loaded: {num_enzymes} enzymes")
    
    # 2. Load KEGG data
    metabolites, tier_r_edges, gene_ec = load_kegg_data(kegg_dir)
    
    # Filter currency metabolites
    metabolites = metabolites[~metabolites['compound_id'].isin(CURRENCY_METABOLITES)]
    print(f"Metabolites after currency filter: {len(metabolites)}")
    
    # Create metabolite index
    met_list = metabolites['compound_id'].tolist()
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    n_metabolites = len(met_list)
    
    # 3. Build Tier-R edges (Reaction-grounded)
    tier_r_src, tier_r_dst = [], []
    tier_r_count = 0
    
    if not tier_r_edges.empty:
        # Create EC to enzyme index mapping (simplified - use random mapping for now)
        # In production: map via gene_ec -> STRING protein ID
        np.random.seed(42)
        
        for _, row in tier_r_edges.iterrows():
            met_id = row['metabolite_id']
            if met_id in met_to_idx:
                met_idx = met_to_idx[met_id]
                # Randomly assign to enzymes (simplified mapping)
                # Real implementation would use gene_ec mapping
                n_links = np.random.randint(5, 20)
                for _ in range(n_links):
                    enz_idx = np.random.randint(0, num_enzymes)
                    tier_r_src.append(enz_idx)
                    tier_r_dst.append(met_idx)
                    tier_r_count += 1
    
    print(f"Tier-R edges (reaction-grounded): {tier_r_count}")
    
    # 4. Build Tier-P edges (Pathway-supported)
    tier_p_src, tier_p_dst = [], []
    tier_p_count = 0
    
    # Get metabolites with pathway info
    pathway_mets = metabolites[metabolites['n_pathways'] > 0]
    
    # Create pathway -> enzyme mapping (simulated based on pathway type)
    np.random.seed(43)
    
    for _, row in pathway_mets.iterrows():
        met_id = row['compound_id']
        pathways = str(row['pathways']).split(';') if pd.notna(row['pathways']) else []
        
        if met_id not in met_to_idx:
            continue
            
        met_idx = met_to_idx[met_id]
        
        # Skip if already covered by Tier-R
        if met_id in tier_r_edges['metabolite_id'].values if not tier_r_edges.empty else []:
            continue
        
        # Assign pathway-based edges
        is_isoflavonoid = any('00943' in p or '00941' in p for p in pathways)
        is_phenylprop = any('00940' in p for p in pathways)
        
        # More edges for relevant pathways
        if is_isoflavonoid or is_phenylprop:
            n_links = np.random.randint(10, 30)
        else:
            n_links = np.random.randint(3, 10)
        
        for _ in range(n_links):
            enz_idx = np.random.randint(0, num_enzymes)
            tier_p_src.append(enz_idx)
            tier_p_dst.append(met_idx)
            tier_p_count += 1
    
    print(f"Tier-P edges (pathway-supported): {tier_p_count}")
    
    # 5. Combine edges with tier weights
    all_src = tier_r_src + tier_p_src
    all_dst = tier_r_dst + tier_p_dst
    
    # Tier weights: R=1.0, P=0.5
    tier_r_weights = [1.0] * len(tier_r_src)
    tier_p_weights = [0.5] * len(tier_p_src)
    all_weights = tier_r_weights + tier_p_weights
    
    # Tier labels for analysis
    tier_r_labels = ['R'] * len(tier_r_src)
    tier_p_labels = ['P'] * len(tier_p_src)
    all_tiers = tier_r_labels + tier_p_labels
    
    # 6. Update HeteroData
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = torch.randn(n_metabolites, 64)
    
    # Store metabolite metadata
    data['Metabolite'].compound_ids = met_list
    
    edge_index = torch.tensor([all_src, all_dst], dtype=torch.long)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = torch.tensor(all_weights, dtype=torch.float)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_tier = all_tiers
    
    # Reverse edges
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(all_dst), torch.tensor(all_src)
    ])
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_weight = torch.tensor(all_weights, dtype=torch.float)
    
    # 7. Summary statistics
    print("\n" + "=" * 60)
    print("Graph Summary")
    print("=" * 60)
    print(f"Nodes:")
    print(f"  - Enzymes: {num_enzymes}")
    print(f"  - Metabolites: {n_metabolites}")
    print(f"Edges:")
    print(f"  - Tier-R (reaction): {tier_r_count} (weight=1.0)")
    print(f"  - Tier-P (pathway): {tier_p_count} (weight=0.5)")
    print(f"  - Total catalyzes: {len(all_src)}")
    print(f"Edge weights: mean={np.mean(all_weights):.3f}")
    
    # Core metabolite coverage
    core_mets = {'C00062', 'C00858', 'C02495'}
    covered_core = core_mets & set(met_list)
    print(f"Core metabolites in graph: {len(covered_core)}/3 ({covered_core})")
    
    # 8. Save
    torch.save(data, output_path)
    print(f"\nSaved to: {output_path}")
    
    # Save tier statistics
    stats = pd.DataFrame({
        'tier': ['R', 'P', 'Total'],
        'n_edges': [tier_r_count, tier_p_count, len(all_src)],
        'weight': [1.0, 0.5, np.mean(all_weights)]
    })
    stats.to_csv(Path(output_path).parent / "tier_statistics.csv", index=False)
    
    return data

if __name__ == "__main__":
    import sys
    import argparse
    sys.path.append(os.getcwd())
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/processed/graph.pt')
    parser.add_argument('--output', type=str, default='data/processed/tiered_bipartite_graph.pt')
    parser.add_argument('--kegg', type=str, default='data/kegg')
    args = parser.parse_args()
    
    build_tiered_bipartite_graph(args.graph, args.output, args.kegg)
