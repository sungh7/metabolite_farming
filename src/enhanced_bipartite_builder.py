"""
Enhanced Bipartite Graph Builder
- Solution 1: NCBI Gene ID direct mapping (via STRING aliases)
- Solution 2: Tier-P pathway-based edges
"""

import torch
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
import os
import gzip
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict

CURRENCY_METABOLITES = {
    'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006',
    'C00008', 'C00009', 'C00010', 'C00011', 'C00013', 'C00014',
    'C00020', 'C00027', 'C00044', 'C00080'
}


def load_string_ncbi_mapping():
    """Load NCBI Gene ID to STRING protein mapping from aliases file."""
    aliases_path = Path("data/raw/3847.protein.aliases.v12.0.txt.gz")
    
    if not aliases_path.exists():
        print(f"Warning: {aliases_path} not found. Falling back to UniProt only.")
        return {}
    
    ncbi_to_string = {}
    print(f"Loading STRING aliases from {aliases_path}...")
    
    with gzip.open(aliases_path, 'rt') as f:
        next(f)  # Skip header
        for line in tqdm(f, desc="Parsing aliases"):
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                string_id = parts[0]  # e.g., "3847.A0A075W8S1"
                alias = parts[1]
                source = parts[2]
                
                # Look for NCBI/Entrez gene IDs
                if 'GeneID' in source or 'NCBI' in source or 'Entrez' in source:
                    ncbi_to_string[alias] = string_id
    
    print(f"NCBI to STRING mappings: {len(ncbi_to_string)}")
    return ncbi_to_string


def build_enhanced_graph(graph_path: str, output_path: str, kegg_dir: str = "data/kegg"):
    """Build graph with NCBI mapping and Tier-P edges."""
    print("=" * 60)
    print("Building Enhanced Bipartite Graph")
    print("Solution 1: NCBI Gene ID mapping")
    print("Solution 2: Tier-P pathway edges")
    print("=" * 60)
    
    kegg_dir = Path(kegg_dir)
    processed_dir = Path("data/processed")
    
    # 1. Load Base PPI Graph
    base_graph = torch.load(graph_path)
    data = base_graph.clone()
    num_enzymes = data['Enzyme'].num_nodes
    print(f"Base graph: {num_enzymes} enzymes")
    
    # 2. Load enzyme_string_mapping.csv
    enzyme_df = pd.read_csv(processed_dir / "enzyme_string_mapping.csv")
    uniprot_to_enzyme_idx = dict(zip(enzyme_df['uniprot_id'], enzyme_df['enzyme_idx']))
    string_to_enzyme_idx = dict(zip(enzyme_df['string_id'], enzyme_df['enzyme_idx']))
    print(f"Enzyme mappings (UniProt): {len(uniprot_to_enzyme_idx)}")
    
    # 3. Solution 1: Load NCBI to STRING mapping
    ncbi_to_string = load_string_ncbi_mapping()
    
    # Also use KEGG gene ID directly if matches STRING format
    # Create NCBI -> Enzyme index
    ncbi_to_enzyme_idx = {}
    for ncbi, string_id in ncbi_to_string.items():
        if string_id in string_to_enzyme_idx:
            ncbi_to_enzyme_idx[ncbi] = string_to_enzyme_idx[string_id]
    print(f"NCBI to Enzyme mappings: {len(ncbi_to_enzyme_idx)}")
    
    # 4. Load KEGG gene-EC mapping
    gene_ec_df = pd.read_csv(kegg_dir / "gene_ec_mapping.tsv", sep='\t')
    print(f"KEGG gene-EC entries: {len(gene_ec_df)}")
    
    # 5. Load KEGG -> UniProt mapping
    kegg_uniprot_df = pd.read_csv(kegg_dir / "kegg_uniprot_mapping.csv")
    kegg_to_uniprot = dict(zip(kegg_uniprot_df['kegg_gene'].astype(str), kegg_uniprot_df['uniprot']))
    
    # 6. Create EC -> Enzyme indices (combining UniProt and NCBI mappings)
    ec_to_enzyme_indices = defaultdict(set)
    matched_uniprot = 0
    matched_ncbi = 0
    unmatched = 0
    
    for _, row in gene_ec_df.iterrows():
        gene_id = str(row['gene_id'])
        ec = row['ec']
        
        enzyme_idx = None
        
        # Try UniProt mapping first
        uniprot = kegg_to_uniprot.get(gene_id)
        if uniprot and uniprot in uniprot_to_enzyme_idx:
            enzyme_idx = uniprot_to_enzyme_idx[uniprot]
            matched_uniprot += 1
        # Fall back to NCBI direct mapping
        elif gene_id in ncbi_to_enzyme_idx:
            enzyme_idx = ncbi_to_enzyme_idx[gene_id]
            matched_ncbi += 1
        else:
            unmatched += 1
        
        if enzyme_idx is not None:
            ec_to_enzyme_indices[ec].add(enzyme_idx)
    
    print(f"Matched via UniProt: {matched_uniprot}")
    print(f"Matched via NCBI: {matched_ncbi}")
    print(f"Unmatched: {unmatched}")
    print(f"ECs with enzyme mapping: {len(ec_to_enzyme_indices)}")
    
    # Count unique enzymes
    unique_enzymes = set()
    for enz_set in ec_to_enzyme_indices.values():
        unique_enzymes.update(enz_set)
    print(f"Unique enzymes with EC: {len(unique_enzymes)}")
    
    # 7. Load full KEGG enzyme-metabolite edges
    edges_df = pd.read_csv(kegg_dir / "full_enzyme_metabolite_edges.tsv", sep='\t')
    print(f"Full KEGG edges: {len(edges_df)}")
    
    # 8. Create metabolite index
    # PRIORITY: Include experimental metabolites
    exp_met_path = Path("data/kegg/metabolites.csv")
    exp_mets = set()
    if exp_met_path.exists():
        exp_df = pd.read_csv(exp_met_path)
        exp_mets = set(exp_df['compound_id'].unique())
        print(f"Experimental metabolites: {len(exp_mets)}")
    
    # Base KEGG metabolites (excluding currency)
    all_kegg_mets = set(edges_df['metabolite_id'].unique()) - CURRENCY_METABOLITES
    
    # Union: Experimental + Top 200 Context
    # Ensure experimental ones correspond to valid KEGG IDs in edges_df? 
    # Not necessarily. If they are in edges_df, good. If not, they might be isolated unless we have edges.
    # We should only include those that HAVE edges in edges_df to avoid zero-degree nodes?
    # Or maybe we need Tier-P edges for them even if no Direct edges.
    
    # Filter exp_mets to those present in edges_df (otherwise we can't link them)
    # Actually, edges_df has all KEGG connections. If exp met is not in edges_df, it has no known enzyme links in KEGG.
    # But Tier-P might link them if we assume class based? No, Tier-P is Enzyme->Metabolite via EC.
    # If Met is not in edges_df, we don't know its ECs.
    
    valid_exp_mets = exp_mets.intersection(set(edges_df['metabolite_id'].unique()))
    print(f"Experimental mets with KEGG edges: {len(valid_exp_mets)} / {len(exp_mets)}")
    
    # Fill remaining slots with top degree metabolites from KEGG
    # Calculate degree for all mets
    met_counts = edges_df['metabolite_id'].value_counts()
    top_context_mets = set(met_counts.index[:300]) - CURRENCY_METABOLITES
    
    final_met_set = valid_exp_mets.union(top_context_mets)
    met_list = sorted(final_met_set)
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    n_metabolites = len(met_list)
    print(f"Final Graph Metabolites: {n_metabolites} (Experimental: {len(valid_exp_mets)})")
    
    # 9. Solution 2: Build pathway membership for Tier-P
    # Group metabolites by EC (shared EC = potential pathway connection)
    ec_to_metabolites = defaultdict(set)
    for _, row in edges_df.iterrows():
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']
        if met_id in met_to_idx:
            ec_to_metabolites[ec].add(met_to_idx[met_id])
    
    # 10. Build edges
    tier_r_edges = set()  # Reaction-grounded (direct EC link)
    tier_p_edges = set()  # Pathway-supported (shared EC class)
    
    for _, row in tqdm(edges_df.iterrows(), total=len(edges_df), desc="Building edges"):
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']
        
        if met_id not in met_to_idx:
            continue
        if ec not in ec_to_enzyme_indices:
            continue
            
        met_idx = met_to_idx[met_id]
        
        # Tier-R: Direct reaction link
        for enz_idx in ec_to_enzyme_indices[ec]:
            tier_r_edges.add((enz_idx, met_idx))
        
        # Tier-P: Enzymes with same EC class (first 3 digits)
        # e.g., EC 1.1.1.x enzymes all relate to similar metabolites
        ec_class = '.'.join(ec.split('.')[:3])  # e.g., "1.1.1"
        for other_ec, other_met_set in ec_to_metabolites.items():
            if other_ec.startswith(ec_class) and other_ec != ec:
                for other_met_idx in other_met_set:
                    for enz_idx in ec_to_enzyme_indices[ec]:
                        tier_p_edges.add((enz_idx, other_met_idx))
    
    print(f"Tier-R edges: {len(tier_r_edges)}")
    print(f"Tier-P edges (raw): {len(tier_p_edges)}")
    
    # Remove Tier-P edges that overlap with Tier-R
    tier_p_edges = tier_p_edges - tier_r_edges
    print(f"Tier-P edges (unique): {len(tier_p_edges)}")
    
    # 11. Combine with weights
    all_edges = []
    all_weights = []
    
    for src, dst in tier_r_edges:
        all_edges.append((src, dst))
        all_weights.append(1.0)  # Tier-R weight
    
    for src, dst in tier_p_edges:
        all_edges.append((src, dst))
        all_weights.append(0.5)  # Tier-P weight
    
    edge_src = [e[0] for e in all_edges]
    edge_dst = [e[1] for e in all_edges]
    
    print(f"Total edges: {len(all_edges)}")
    
    # 12. Update HeteroData
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = torch.randn(n_metabolites, 64)
    data['Metabolite'].compound_ids = met_list
    
    edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long)
    edge_weights_tensor = torch.tensor(all_weights, dtype=torch.float)
    
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = edge_weights_tensor
    
    # Reverse edges
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(edge_dst), torch.tensor(edge_src)
    ])
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_weight = edge_weights_tensor
    
    # 13. Verify indices
    if edge_src:
        max_enz_idx = max(edge_src)
        print(f"Max enzyme index: {max_enz_idx}, Enzyme num_nodes: {num_enzymes}")
        assert max_enz_idx < num_enzymes, f"Index out of bounds"
    
    # 14. Summary
    print("\n" + "=" * 60)
    print("Graph Summary")
    print("=" * 60)
    print(f"Enzymes: {num_enzymes}")
    print(f"Metabolites: {n_metabolites}")
    print(f"Tier-R edges: {len(tier_r_edges)}")
    print(f"Tier-P edges: {len(tier_p_edges)}")
    print(f"Total edges: {len(all_edges)}")
    print(f"Mean edge weight: {np.mean(all_weights):.3f}")
    
    torch.save(data, output_path)
    print(f"\nSaved to: {output_path}")
    
    return data


if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    
    build_enhanced_graph(
        'data/processed/strict_graph.pt',
        'data/processed/enhanced_bipartite_graph.pt'
    )
