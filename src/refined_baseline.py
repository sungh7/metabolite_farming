import torch
import torch_geometric
import pandas as pd
import numpy as np
import random
from src.bipartite_builder import build_bipartite_graph
from torch_geometric.transforms import RandomLinkSplit

def run_baseline_eval():
    torch.manual_seed(42)
    np.random.seed(42)
    random.seed(42)
    print("Building Graph...")
    # data = build_bipartite_graph() # No need to rebuild if we load the same PT file as trainer
    data = torch.load('data/processed/bipartite_graph.pt')
    
    # --- Replicate Split Logic (from refined_trainer.py) ---
    num_enzymes = data['Enzyme'].num_nodes
    indices = torch.randperm(num_enzymes) # Standard torch generator (seeded 42)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool)
    train_enz_mask[indices[:split]] = True
    test_enz_mask = ~train_enz_mask
    
    print(f"Split: {train_enz_mask.sum()} Train, {test_enz_mask.sum()} Test")
    
    # Identify Test Edges
    # Edge type: ('Enzyme', 'catalyzes', 'Metabolite')
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    
    src = edge_index[0]
    mask = train_enz_mask[src]
    
    # Test edges are those connected to Test Enzymes
    test_edges = edge_index[:, ~mask]
    # In bipartite_builder, src=Enzyme.
    
    # Trainer used:
    # train_edges = edge_index[:, mask]
    # test_edges = edge_index[:, ~mask]
    
    # Evaluate rank of Metabolite -> Enzyme?
    # Task is Metabolite (Query) -> Enzyme (Candidate).
    # Since undirected, we can view it either way.
    # Trainer evaluate() did:
    # Given Metabolite, Rank Test Enzymes.
    
    # Test Edges: [Enzyme, Metabolite]
    pos_edges = test_edges
    
    # Create the candidate pool (All Test Enzymes)
    test_enzyme_indices = torch.arange(num_enzymes)[test_enz_mask]
    all_enzymes = np.arange(num_enzymes)
    
    # import networkx as nx # causing syntax error
    
    # --- Network Baseline Setup (Custom BFS) ---
    print("Building Adjacency Matrix from Training Edges...")
    
    # Nodes: 0..N_Enz-1 (Enzymes), N_Enz..N_Enz+N_Met-1 (Metabolites)
    # Build Adj List
    # adj[node] = [neighbors]
    
    adj = {}
    
    train_edge_index = edge_index[:, mask]
    
    src_list = train_edge_index[0].tolist()
    dst_list = (train_edge_index[1] + num_enzymes).tolist() # Offset Metabolites
    
    for u, v in zip(src_list, dst_list):
        if u not in adj: adj[u] = []
        adj[u].append(v)
        if v not in adj: adj[v] = []
        adj[v].append(u)
        
    def get_shortest_paths(start_node, cutoff=4):
        # BFS
        dists = {start_node: 0}
        queue = [(start_node, 0)]
        
        while queue:
            u, d = queue.pop(0)
            if d >= cutoff: continue
            
            if u in adj:
                for v in adj[u]:
                    if v not in dists:
                        dists[v] = d + 1
                        queue.append((v, d + 1))
        return dists

    # --- Evaluation ---
    hits_20_random = 0
    hits_20_network = 0
    total = pos_edges.shape[1]
    
    print(f"Evaluating {total} test edges (Custom BFS)...")
    
    for i in range(total):
        u, v = pos_edges[:, i]
        m_idx = u.item() # Local Met Index
        true_e_idx = v.item() # Local Enz Index
        
        # Node IDs
        bfs_m_id = num_enzymes + m_idx
        
        # 1. Random
        if np.random.rand() < 20 / num_enzymes:
            hits_20_random += 1
            
        # 2. Network Proximity (BFS)
        # Calculate distances from Metabolite M to all nodes
        dists = get_shortest_paths(bfs_m_id, cutoff=5)
        
        candidate_scores = []
        for cand_e in test_enzyme_indices:
            cand_id = cand_e.item()
            d = dists.get(cand_id, 999)
            # Score = 1/d (Large score = close)
            # If disconnected (999), score 0
            score = 1.0 / d if d > 0 else 0.0
            candidate_scores.append(score)
            
        candidate_scores = np.array(candidate_scores)
        
        # True is at idx?
        true_pos_idx = (test_enzyme_indices == true_e_idx).nonzero(as_tuple=True)[0]
        if len(true_pos_idx) == 0: continue
        true_pos_idx = true_pos_idx.item()
        
        # Rank
        argsort = np.argsort(-candidate_scores) # Descending
        rank = np.where(argsort == true_pos_idx)[0][0] + 1
        
        if rank <= 20:
            hits_20_network += 1
            
    print("-" * 30)
    print(f"Total Test Edges: {total}")
    print(f"Hits@20 (Random): {hits_20_random / total * 100:.2f}%")
    print(f"Hits@20 (Network Proximity): {hits_20_network / total * 100:.2f}%")
    
    with open('results/baseline_metrics.txt', 'w') as f:
        f.write(f"Hits@20_Random: {hits_20_random / total * 100:.2f}\n")
        f.write(f"Hits@20_NetworkProximity: {hits_20_network / total * 100:.2f}\n")
    print("Saved results/baseline_metrics.txt")

if __name__ == "__main__":
    run_baseline_eval()
