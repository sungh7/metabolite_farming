"""
DEPRECATED: Simulated Bipartite Graph Builder

This module uses simulated enzyme-metabolite edges and is kept for
backward compatibility. For production use, please use:
    - src/data_pipeline.py (recommended)
    - src/enhanced_bipartite_builder.py (real KEGG data)

The simulated edges in this module do not reflect actual biochemical
relationships and should only be used for initial testing.
"""

import warnings
import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
import os
from src.dataloader import StringDBLoader

# Issue deprecation warning on import
warnings.warn(
    "bipartite_builder.py is deprecated. Use data_pipeline.py or "
    "enhanced_bipartite_builder.py instead for real KEGG data.",
    DeprecationWarning,
    stacklevel=2
)

def build_bipartite_graph(graph_path, output_path):
    """
    Constructs Heterogeneous Graph (PPI + Metabolite-Enzyme).
    Args:
        graph_path (str): Path to the base PPI graph (Strict or Full).
        output_path (str): Path to save the merged graph.
    """
    print(f"Building Bipartite Graph from {graph_path}...")
    
    # 1. Load Base PPI Graph
    if not os.path.exists(graph_path):
        print(f"Error: {graph_path} not found.")
        return
        
    base_graph = torch.load(graph_path)
    
    # Clone to avoid modifying original?
    data = base_graph.clone()
    
    # 2. Load Metabolites (Metadata)
    # We use the counts from our sim/data.
    # Note: If real integration, we should map real metabolites.
    # For MVP, we stick to the simulated/mapped logic for 'Metabolite' nodes,
    # assuming we have a mapping of Enzyme ID -> Reactant ID.
    
    # For this MVP code, we are using the simulated bipartite logic from before
    # BUT we must respect the Enzyme Nodes existing in base_graph.
    
    num_enzymes = data['Enzyme'].num_nodes
    
    # Re-generate simulated metabolite logic consistent with previous steps
    n_metabolites = 200
    compounds = [f"Metabolite_{i}" for i in range(n_metabolites)]
    met_map = {name: i for i, name in enumerate(compounds)}
    
    # Assign Pathways (Simulated)
    met_pathways = {}
    for i in range(20): met_pathways[compounds[i]] = 'Phenylpropanoid'
    for i in range(20, 50): met_pathways[compounds[i]] = 'Flavonoid'
    for i in range(50, n_metabolites): met_pathways[compounds[i]] = 'Other'
    
    # 3. Create Edges
    edges_src = [] # Enzyme
    edges_dst = [] # Metabolite
    
    np.random.seed(42)
    
    enz_pathways = {}
    for i in range(num_enzymes):
        r = np.random.rand()
        if r < 0.05: p = 'Phenylpropanoid'
        elif r < 0.15: p = 'Flavonoid'
        else: p = 'Other'
        enz_pathways[i] = p
        
    for enz_idx, enz_path in enz_pathways.items():
        for met_name, met_path in met_pathways.items():
            if enz_path == met_path and met_path != 'Other':
                if np.random.rand() < 0.1: 
                   edges_src.append(enz_idx)
                   edges_dst.append(met_map[met_name])
                   
    # 4. Update HeteroData
    # Add Metabolite Nodes
    data['Metabolite'].num_nodes = n_metabolites
    # Initialize features (Random or Chem descriptors)
    data['Metabolite'].x = torch.randn(n_metabolites, 64)
    
    # Add Bipartite Edges with MSI Level-based Weights
    edge_index = torch.tensor([edges_src, edges_dst], dtype=torch.long)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index
    
    # MSI Level -> Weight Mapping (Phase 5 Evidence Weighting)
    # Level 2 (Strong evidence, e.g., Phenylpropanoid/Flavonoid pathways): 1.0
    # Level 3 (Moderate evidence): 0.7
    # Level 4 (Unknown/Other): 0.4
    def get_msi_weight(met_idx):
        met_name = compounds[met_idx]
        pathway = met_pathways.get(met_name, 'Other')
        if pathway in ['Phenylpropanoid', 'Flavonoid']:
            return 1.0  # Level 2 equivalent
        else:
            return 0.4  # Level 3-4 equivalent
    
    edge_weights = torch.tensor([get_msi_weight(dst) for dst in edges_dst], dtype=torch.float)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = edge_weights
    print(f"Edge weights assigned: Mean={edge_weights.mean():.3f}, Min={edge_weights.min():.3f}, Max={edge_weights.max():.3f}")
    
    # Reverse edge (Metabolite -> Enzyme) useful for message passing
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([torch.tensor(edges_dst), torch.tensor(edges_src)])
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_weight = edge_weights  # Same weights
    
    print(f"Merged Graph: {data}")
    torch.save(data, output_path)
    print(f"Saved to {output_path}")

if __name__ == "__main__":
    import sys
    import argparse
    sys.path.append(os.getcwd())
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/processed/graph.pt')
    parser.add_argument('--output', type=str, default='data/processed/bipartite_graph.pt')
    args = parser.parse_args()
    
    build_bipartite_graph(args.graph, args.output)
