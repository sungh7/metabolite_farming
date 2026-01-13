"""
Proper Bipartite Graph Builder with Real Gene Mappings (Fixed)
Uses enzyme_string_mapping.csv to get correct enzyme local indices
"""

import torch
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
import os
from pathlib import Path
from tqdm import tqdm

CURRENCY_METABOLITES = {
    'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006',
    'C00008', 'C00009', 'C00010', 'C00011', 'C00013', 'C00014',
    'C00020', 'C00027', 'C00044', 'C00080'
}


def build_proper_graph(graph_path: str, output_path: str, kegg_dir: str = "data/kegg"):
    """Build graph with proper gene ID mapping using enzyme_string_mapping.csv."""
    print("=" * 60)
    print("Building Proper Bipartite Graph with Gene Mappings (Fixed)")
    print("=" * 60)
    
    kegg_dir = Path(kegg_dir)
    processed_dir = Path("data/processed")
    
    # 1. Load Base PPI Graph
    base_graph = torch.load(graph_path)
    data = base_graph.clone()
    num_enzymes = data['Enzyme'].num_nodes
    print(f"Base graph: {num_enzymes} enzymes")
    
    # 2. Load enzyme_string_mapping.csv (enzyme_idx, string_id, uniprot_id)
    mapping_path = processed_dir / "enzyme_string_mapping.csv"
    if not mapping_path.exists():
        raise FileNotFoundError(f"enzyme_string_mapping.csv not found. Run graph_builder.py first.")
    
    enzyme_df = pd.read_csv(mapping_path)
    uniprot_to_enzyme_idx = dict(zip(enzyme_df['uniprot_id'], enzyme_df['enzyme_idx']))
    print(f"Enzyme mappings: {len(uniprot_to_enzyme_idx)}")
    
    # 3. Load KEGG gene-EC mapping
    gene_ec_path = kegg_dir / "gene_ec_mapping.tsv"
    gene_ec_df = pd.read_csv(gene_ec_path, sep='\t')
    print(f"KEGG gene-EC entries: {len(gene_ec_df)}")
    
    # 4. Load KEGG -> UniProt mapping (cached)
    mapping_cache = kegg_dir / "kegg_uniprot_mapping.csv"
    if not mapping_cache.exists():
        raise FileNotFoundError("kegg_uniprot_mapping.csv not found. Run proper_bipartite_builder.py first to fetch.")
    
    kegg_uniprot_df = pd.read_csv(mapping_cache)
    kegg_to_uniprot = dict(zip(kegg_uniprot_df['kegg_gene'].astype(str), kegg_uniprot_df['uniprot']))
    print(f"KEGG-UniProt mappings: {len(kegg_to_uniprot)}")
    
    # 5. Create EC -> Enzyme indices mapping (using real gene IDs)
    ec_to_enzyme_indices = {}
    matched = 0
    unmatched = 0
    
    for _, row in gene_ec_df.iterrows():
        gene_id = str(row['gene_id'])
        ec = row['ec']
        
        # Get UniProt ID from KEGG gene
        uniprot = kegg_to_uniprot.get(gene_id)
        if uniprot and uniprot in uniprot_to_enzyme_idx:
            enzyme_idx = uniprot_to_enzyme_idx[uniprot]
            if ec not in ec_to_enzyme_indices:
                ec_to_enzyme_indices[ec] = []
            if enzyme_idx not in ec_to_enzyme_indices[ec]:
                ec_to_enzyme_indices[ec].append(enzyme_idx)
            matched += 1
        else:
            unmatched += 1
    
    print(f"Gene-EC matches: {matched}, unmatched: {unmatched}")
    print(f"ECs with enzyme mapping: {len(ec_to_enzyme_indices)}")
    
    # Count unique enzymes mapped
    unique_enzymes = set()
    for enz_list in ec_to_enzyme_indices.values():
        unique_enzymes.update(enz_list)
    print(f"Unique enzymes with EC: {len(unique_enzymes)}")
    
    # 6. Load full KEGG enzyme-metabolite edges
    edges_path = kegg_dir / "full_enzyme_metabolite_edges.tsv"
    edges_df = pd.read_csv(edges_path, sep='\t')
    print(f"Full KEGG edges: {len(edges_df)}")
    
    # 7. Create metabolite index
    all_kegg_mets = set(edges_df['metabolite_id'].unique()) - CURRENCY_METABOLITES
    met_list = sorted(all_kegg_mets)[:200]  # Limit to 200
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    n_metabolites = len(met_list)
    print(f"Metabolites: {n_metabolites}")
    
    # 8. Build edges with real mapping
    edge_src, edge_dst = [], []
    
    for _, row in tqdm(edges_df.iterrows(), total=len(edges_df), desc="Building edges"):
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']
        
        if met_id not in met_to_idx:
            continue
        if ec not in ec_to_enzyme_indices:
            continue
            
        met_idx = met_to_idx[met_id]
        
        for enz_idx in ec_to_enzyme_indices[ec]:
            edge_src.append(enz_idx)
            edge_dst.append(met_idx)
    
    print(f"Total edges: {len(edge_src)}")
    
    # 9. Remove duplicates
    edge_set = set(zip(edge_src, edge_dst))
    edge_src, edge_dst = zip(*edge_set) if edge_set else ([], [])
    edge_src, edge_dst = list(edge_src), list(edge_dst)
    print(f"Unique edges: {len(edge_src)}")
    
    # 10. Update HeteroData
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = torch.randn(n_metabolites, 64)
    data['Metabolite'].compound_ids = met_list
    
    edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long)
    edge_weights = torch.ones(len(edge_src), dtype=torch.float)
    
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = edge_weights
    
    # Reverse edges
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(edge_dst), torch.tensor(edge_src)
    ])
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_weight = edge_weights
    
    # 11. Verify indices
    max_enz_idx = max(edge_src)
    print(f"Max enzyme index: {max_enz_idx}, Enzyme num_nodes: {num_enzymes}")
    assert max_enz_idx < num_enzymes, f"Index out of bounds: {max_enz_idx} >= {num_enzymes}"
    
    # 12. Summary
    print("\n" + "=" * 60)
    print("Graph Summary")
    print("=" * 60)
    print(f"Enzymes: {num_enzymes}")
    print(f"Metabolites: {n_metabolites}")
    print(f"Edges: {len(edge_src)}")
    
    torch.save(data, output_path)
    print(f"\nSaved to: {output_path}")
    
    return data


if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    
    build_proper_graph(
        'data/processed/strict_graph.pt',
        'data/processed/proper_bipartite_graph.pt'
    )
