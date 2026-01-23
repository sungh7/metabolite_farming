"""
Legacy Bipartite Graph Builder (Original Implementation)

This is the original bipartite_builder.py preserved for backward compatibility.
It uses SIMULATED enzyme-metabolite edges for initial prototyping.

DO NOT USE FOR PRODUCTION - Use src/data_pipeline.py instead.

Original behavior:
- Simulates 200 metabolites
- Creates random pathway assignments
- Generates edges based on pathway matching with 10% probability
- MSI-level based edge weights
"""

import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
import os


def build_bipartite_graph(graph_path, output_path):
    """
    Constructs Heterogeneous Graph (PPI + Metabolite-Enzyme).

    LEGACY VERSION - Uses simulated data.

    Args:
        graph_path (str): Path to the base PPI graph (Strict or Full).
        output_path (str): Path to save the merged graph.
    """
    print("[LEGACY] Building Bipartite Graph from simulated data...")
    print(f"Loading base graph from {graph_path}...")

    # 1. Load Base PPI Graph
    if not os.path.exists(graph_path):
        print(f"Error: {graph_path} not found.")
        return

    base_graph = torch.load(graph_path)
    data = base_graph.clone()
    num_enzymes = data['Enzyme'].num_nodes

    # 2. Load Metabolites (Simulated)
    n_metabolites = 200
    compounds = [f"Metabolite_{i}" for i in range(n_metabolites)]
    met_map = {name: i for i, name in enumerate(compounds)}

    # Assign Pathways (Simulated)
    met_pathways = {}
    for i in range(20):
        met_pathways[compounds[i]] = 'Phenylpropanoid'
    for i in range(20, 50):
        met_pathways[compounds[i]] = 'Flavonoid'
    for i in range(50, n_metabolites):
        met_pathways[compounds[i]] = 'Other'

    # 3. Create Edges
    edges_src = []  # Enzyme
    edges_dst = []  # Metabolite

    np.random.seed(42)

    enz_pathways = {}
    for i in range(num_enzymes):
        r = np.random.rand()
        if r < 0.05:
            p = 'Phenylpropanoid'
        elif r < 0.15:
            p = 'Flavonoid'
        else:
            p = 'Other'
        enz_pathways[i] = p

    for enz_idx, enz_path in enz_pathways.items():
        for met_name, met_path in met_pathways.items():
            if enz_path == met_path and met_path != 'Other':
                if np.random.rand() < 0.1:
                    edges_src.append(enz_idx)
                    edges_dst.append(met_map[met_name])

    # 4. Update HeteroData
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].x = torch.randn(n_metabolites, 64)

    # Add Bipartite Edges with MSI Level-based Weights
    edge_index = torch.tensor([edges_src, edges_dst], dtype=torch.long)
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = edge_index

    # MSI Level -> Weight Mapping
    def get_msi_weight(met_idx):
        met_name = compounds[met_idx]
        pathway = met_pathways.get(met_name, 'Other')
        if pathway in ['Phenylpropanoid', 'Flavonoid']:
            return 1.0  # Level 2 equivalent
        else:
            return 0.4  # Level 3-4 equivalent

    edge_weights = torch.tensor(
        [get_msi_weight(dst) for dst in edges_dst],
        dtype=torch.float
    )
    data['Enzyme', 'catalyzes', 'Metabolite'].edge_weight = edge_weights
    print(f"Edge weights: Mean={edge_weights.mean():.3f}")

    # Reverse edge
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = torch.stack([
        torch.tensor(edges_dst), torch.tensor(edges_src)
    ])
    data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_weight = edge_weights

    print(f"Merged Graph: {data}")
    torch.save(data, output_path)
    print(f"Saved to {output_path}")

    return data


if __name__ == "__main__":
    import sys
    import argparse
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/processed/graph.pt')
    parser.add_argument('--output', type=str, default='data/processed/legacy_bipartite_graph.pt')
    args = parser.parse_args()

    build_bipartite_graph(args.graph, args.output)
