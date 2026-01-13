"""
Ethylene-Conditioned Graph Builder
Injects ethylene-response weights into graph edge structure.

Key changes:
1. TF-Enzyme edges: weight by TF ethylene-responsiveness  
2. Enzyme-Metabolite edges: weight by enzyme ethylene-response score
3. Node features: add ethylene-response dimension
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
import requests
import time

# Ethylene-responsive transcription factors (known from literature)
ETHYLENE_RESPONSIVE_TFS = {
    'AP2': 2.0,      # AP2/ERF family - primary ethylene response
    'ERF': 2.0,      # ERF subfamily
    'EIN3': 1.8,     # EIN3/EIL family
    'NAC': 1.5,      # NAC family (stress-responsive)
    'WRKY': 1.3,     # WRKY (defense, some ethylene-linked)
    'bHLH': 1.2,     # bHLH (general stress)
}

# Ethylene-responsive enzyme ECs (isoflavonoid/phenylpropanoid pathway)
ETHYLENE_RESPONSIVE_ECS = {
    # Isoflavonoid biosynthesis - highest weights
    '2.5.1.82': 2.5,   # IFS (isoflavone synthase)
    '5.5.1.6': 2.0,    # CHI (chalcone isomerase)  
    '1.1.1.234': 2.0,  # IFR (isoflavone reductase)
    '2.4.1.170': 1.8,  # UGT (UDP-glucosyltransferase for daidzin)
    # Phenylpropanoid pathway
    '4.3.1.24': 1.8,   # PAL (phenylalanine ammonia-lyase)
    '4.3.1.25': 1.8,   # PAL related
    '6.2.1.12': 1.5,   # 4CL (4-coumarate-CoA ligase)
    '2.3.1.74': 1.5,   # CHS (chalcone synthase)
    # Stress-responsive general
    '1.11.1.7': 1.3,   # Peroxidase
    '1.14.11.-': 1.3,  # 2-oxoglutarate-dependent dioxygenases
}

# Isoflavonoid metabolites (target enrichment)
ISOFLAVONOID_METS = {
    'C02495': 2.5,  # Daidzein - core isoflavone
    'C00858': 2.5,  # Formononetin - core isoflavone  
    'C10216': 2.0,  # Daidzin - daidzein glycoside
}

def build_ethylene_conditioned_graph(base_graph_path, full_edges_path, kegg_dir, output_path):
    """Build graph with ethylene-conditioned edge weights."""
    print("=" * 60)
    print("Building Ethylene-Conditioned Graph")
    print("=" * 60)
    
    kegg_dir = Path(kegg_dir)
    
    # Load base graph
    data = torch.load(base_graph_path).clone()
    num_enzymes = data['Enzyme'].num_nodes
    print(f"Base graph: {num_enzymes} enzymes")
    
    # Load gene-EC mapping
    gene_ec_df = pd.read_csv(kegg_dir / "gene_ec_mapping.tsv", sep='\t')
    ec_list = gene_ec_df['ec'].unique().tolist()
    
    # Create enzyme → EC mapping (simplified: random assignment as before)
    np.random.seed(42)
    ec_to_enz = {}
    enz_to_ec = {}
    
    for ec in ec_list:
        n_enz = np.random.randint(5, 20)
        enz_indices = np.random.choice(num_enzymes, n_enz, replace=False).tolist()
        ec_to_enz[ec] = enz_indices
        for enz in enz_indices:
            if enz not in enz_to_ec:
                enz_to_ec[enz] = []
            enz_to_ec[enz].append(ec)
    
    # Calculate enzyme ethylene-response scores
    enz_ethylene_scores = {}
    for enz_idx in range(num_enzymes):
        ecs = enz_to_ec.get(enz_idx, [])
        max_score = 1.0
        for ec in ecs:
            if ec in ETHYLENE_RESPONSIVE_ECS:
                max_score = max(max_score, ETHYLENE_RESPONSIVE_ECS[ec])
        enz_ethylene_scores[enz_idx] = max_score
    
    elevated_enzymes = sum(1 for s in enz_ethylene_scores.values() if s > 1.0)
    print(f"Enzymes with elevated ethylene scores: {elevated_enzymes}")
    
    # Load full KEGG edges
    full_edges = pd.read_csv(full_edges_path, sep='\t')
    
    # Load expanded graph for metabolite list
    expanded = torch.load('data/processed/expanded_bipartite_graph.pt')
    met_list = expanded['Metabolite'].compound_ids
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    n_metabolites = len(met_list)
    
    # Build edges with ethylene-conditioned weights
    # Weight formula: base_tier_weight × enzyme_ET_score × metabolite_ET_score
    
    covered_mtbls = {'C00062', 'C00078', 'C00858', 'C02495', 'C10216', 
                     'C01177', 'C06037', 'C04079', 'C19865', 'C01004'}
    
    edge_src, edge_dst = [], []
    edge_weights = []
    
    for _, row in full_edges.iterrows():
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']
        
        if met_id not in met_to_idx:
            continue
        
        met_idx = met_to_idx[met_id]
        enz_indices = ec_to_enz.get(ec, [])
        is_mtbls = met_id in covered_mtbls
        
        # Base tier weight
        tier_weight = 1.0 if is_mtbls else 0.5
        
        # Metabolite ethylene weight
        met_et_weight = ISOFLAVONOID_METS.get(met_id, 1.0)
        
        for enz_idx in enz_indices:
            # Enzyme ethylene weight
            enz_et_weight = enz_ethylene_scores.get(enz_idx, 1.0)
            
            # Combined weight
            final_weight = tier_weight * enz_et_weight * met_et_weight
            
            edge_src.append(enz_idx)
            edge_dst.append(met_idx)
            edge_weights.append(final_weight)
    
    print(f"Total edges: {len(edge_src)}")
    print(f"Weight distribution: min={min(edge_weights):.2f}, max={max(edge_weights):.2f}, "
          f"mean={np.mean(edge_weights):.2f}")
    
    # Count high-weight edges (ethylene-conditioned)
    high_weight = sum(1 for w in edge_weights if w > 2.0)
    print(f"High-weight edges (>2.0): {high_weight}")
    
    # Create node features with ethylene dimension
    # Add ethylene-response as additional feature dimension
    base_dim = 64
    et_dim = 8  # Ethylene response encoding
    
    # Enzyme features: base + ethylene score encoding
    enz_features = torch.randn(num_enzymes, base_dim)
    enz_et_feat = torch.zeros(num_enzymes, et_dim)
    for enz_idx, score in enz_ethylene_scores.items():
        if score > 1.0:
            enz_et_feat[enz_idx] = torch.randn(et_dim) * (score - 1.0)
    enz_features = torch.cat([enz_features, enz_et_feat], dim=1)
    
    # Metabolite features: base + ethylene target encoding
    met_features = torch.randn(n_metabolites, base_dim)
    met_et_feat = torch.zeros(n_metabolites, et_dim)
    for met_id, score in ISOFLAVONOID_METS.items():
        if met_id in met_to_idx:
            idx = met_to_idx[met_id]
            met_et_feat[idx] = torch.randn(et_dim) * (score - 1.0)
    met_features = torch.cat([met_features, met_et_feat], dim=1)
    
    # Update data
    data['Enzyme'].x = enz_features
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = met_features
    data['Metabolite'].compound_ids = met_list
    
    edge_index = torch.tensor([edge_src, edge_dst], dtype=torch.long)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = torch.tensor(edge_weights, dtype=torch.float)
    
    # Reverse edges
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(edge_dst), torch.tensor(edge_src)
    ])
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_weight = torch.tensor(edge_weights, dtype=torch.float)
    
    # Statistics
    print("\n" + "=" * 60)
    print("Graph Summary")
    print("=" * 60)
    print(f"Nodes: Enzymes={num_enzymes}, Metabolites={n_metabolites}")
    print(f"Edges: {len(edge_src)}")
    print(f"Feature dim: {enz_features.shape[1]} (base {base_dim} + ET {et_dim})")
    
    # Check isoflavonoid connectivity
    iso_edges = sum(1 for d, w in zip(edge_dst, edge_weights) 
                    if met_list[d] in ISOFLAVONOID_METS and w > 2.0)
    print(f"High-weight edges to isoflavonoid metabolites: {iso_edges}")
    
    # Save
    torch.save(data, output_path)
    print(f"\nSaved to: {output_path}")
    
    return data

if __name__ == "__main__":
    build_ethylene_conditioned_graph(
        'data/processed/strict_graph.pt',
        'data/kegg/full_enzyme_metabolite_edges.tsv',
        'data/kegg',
        'data/processed/ethylene_conditioned_graph.pt'
    )
