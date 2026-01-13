import torch
import torch_geometric.transforms as T
from src.dataloader import StringDBLoader
import pandas as pd
import os
import gzip
from collections import deque

def find_paths_bfs(adj, start_node, end_node, max_depth=2):
    """
    Simple BFS to find all paths from start to end with length <= max_depth
    adj: dict of set {u: {v1, v2, ...}}
    """
    print(f"Searching paths depth {max_depth} from {start_node} to {end_node}")
    queue = deque([(start_node, [start_node])])
    paths = []
    
    while queue:
        curr, path = queue.popleft()
        
        if len(path) > max_depth + 1:
            continue
            
        if curr == end_node:
            if len(path) > 1: # Don't allow self-loop path if start==end (unlikely)
                paths.append(path)
            continue
        
        if len(path) == max_depth + 1:
            continue

        for neighbor in adj.get(curr, []):
            if neighbor not in path: # Avoid cycles in simple path
                queue.append((neighbor, path + [neighbor]))
                
    return paths

def visualize_candidate_subgraph(rank=1):
    candidates = pd.read_csv('results/candidates.csv')
    candidate = candidates.iloc[rank-1]
    
    tf_name = candidate['TF_Name']
    enz_name = candidate['Enzyme_Name']
    tf_id = candidate['TF_ID']
    enz_id = candidate['Enzyme_ID']
    score = candidate['Score']
    
    print(f"Generating Mermaid for Rank {rank}: {tf_name} <--> {enz_name}")
    
    # 1. Load Mappings (Heavy lift)
    loader = StringDBLoader()
    loader.load_protein_info()
    
    # Load Graph
    data = torch.load('data/processed/graph.pt')
    
    # Reconstruct Metadata for Type Mapping
    # (Simplified: Just load names to use in Mermaid)
    # We still need to know WHICH node in 'graph.pt' corresponds to our TF_ID/Enzyme_ID
    # We used candidates.py logic.
    
    # To save time, we will assume we can reconstruct the ID mapping exactly.
    # Sorted String IDs -> Assign Types.
    
    sorted_sids = [None] * len(loader.node_to_idx)
    for sid, idx in loader.node_to_idx.items():
        sorted_sids[idx] = sid
        
    # We need annotation to identify types
    node_annotations = {}
    with gzip.open(loader.protein_info_path, 'rt') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                node_annotations[parts[0]] = parts[3]
            elif len(parts) >= 2:
                 node_annotations[parts[0]] = ""
                 
    from src.graph_builder import identify_node_type
    
    # Global ID -> (Type, LocalIndex)
    global_to_typeidx = {}
    # Also TypeIdx -> Global ID (for path decoding)
    typeidx_to_global = {'Signaling': {}, 'TF': {}, 'Enzyme': {}, 'Protein': {}}
    
    type_counts = {'Signaling': 0, 'TF': 0, 'Enzyme': 0, 'Protein': 0}
    
    for gid, sid in enumerate(sorted_sids):
        name = loader.protein_map.get(sid, "")
        ann = node_annotations.get(sid, "")
        ntype = identify_node_type(ann, name)
        
        lidx = type_counts[ntype]
        type_counts[ntype] += 1
        
        global_to_typeidx[gid] = (ntype, lidx)
        typeidx_to_global[ntype][lidx] = gid
        
    # Identify Start/End Node in Graph Terms
    # TF_ID is a string. We need its Global ID first.
    if tf_id not in loader.node_to_idx or enz_id not in loader.node_to_idx:
        print("Candidate IDs not found in loader map.")
        return
        
    start_gid = loader.node_to_idx[tf_id]
    end_gid = loader.node_to_idx[enz_id]
    
    start_node_key = f"{global_to_typeidx[start_gid][0]}_{global_to_typeidx[start_gid][1]}"
    end_node_key = f"{global_to_typeidx[end_gid][0]}_{global_to_typeidx[end_gid][1]}"
    
    print(f"Searching path: {start_node_key} -> {end_node_key}")
    
    # 2. Build Adjacency List (Homogeneous)
    adj = {}
    
    edge_types = data.edge_types
    for (st, rel, dt) in edge_types:
        edge_index = data[st, rel, dt].edge_index
        srcs = edge_index[0].numpy()
        dsts = edge_index[1].numpy()
        
        for s, d in zip(srcs, dsts):
            u = f"{st}_{s}"
            v = f"{dt}_{d}"
            
            if u not in adj: adj[u] = set()
            if v not in adj: adj[v] = set()
            
            adj[u].add(v)
            adj[v].add(u) # Undirected
            
    # 3. Find Paths
    paths = find_paths_bfs(adj, start_node_key, end_node_key, max_depth=3)
    print(f"Found {len(paths)} paths.")
    
    # 4. Generate Mermaid
    mmd = ["graph LR"]
    
    # Style
    mmd.append("    classDef tf fill:#a2d2ff,stroke:#333,stroke-width:2px;")
    mmd.append("    classDef enz fill:#b5e48c,stroke:#333,stroke-width:2px;")
    mmd.append("    classDef sig fill:#fec89a,stroke:#333,stroke-width:2px;")
    mmd.append("    classDef prot fill:#e5e5e5,stroke:#333,stroke-width:1px;")
    
    added_nodes = set()
    added_edges = set()
    
    # Helper to get name
    def get_info(node_key):
        t, idx = node_key.split('_')
        idx = int(idx)
        gid = typeidx_to_global[t][idx]
        sid = sorted_sids[gid]
        name = loader.protein_map.get(sid, sid)
        return sid, name, t
        
    # Limit number of paths for clarity
    display_paths = paths[:10]
    
    for path in display_paths:
        for i in range(len(path)-1):
            u = path[i]
            v = path[i+1]
            
            edge = tuple(sorted((u, v)))
            if edge in added_edges:
                continue
            added_edges.add(edge)
            
            # Add nodes
            for node in [u, v]:
                if node not in added_nodes:
                    sid, name, t = get_info(node)
                    # Clean name for mermaid
                    safe_name = name.replace("(", "").replace(")", "").replace(".", "_")
                    
                    style = "prot"
                    if t == 'TF': style = 'tf'
                    elif t == 'Enzyme': style = 'enz'
                    elif t == 'Signaling': style = 'sig'
                    
                    # Highlight targets
                    if node == start_node_key: style = "tf;stroke:red"
                    if node == end_node_key: style = "enz;stroke:red"
                    
                    mmd.append(f"    {node}({safe_name}):::{style}")
                    added_nodes.add(node)
            
            mmd.append(f"    {u} --- {v}")
            
    if not display_paths:
        mmd.append(f"    {start_node_key} -.- {end_node_key}")
        
    content = "\n".join(mmd)
    
    os.makedirs('results/figures', exist_ok=True)
    with open(f'results/figures/explanation_rank_{rank}.mmd', 'w') as f:
        f.write(content)
        
    print(f"Saved mermaid to results/figures/explanation_rank_{rank}.mmd")
    print(content)
    return paths

if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    
    # Try ranks 1 to 5 until we find a path
    for r in range(1, 6):
        print(f"--- Trying Rank {r} ---")
        paths = visualize_candidate_subgraph(rank=r)
        if paths and len(paths) > 0:
            print(f"Success! Found paths for Rank {r}")
            # We already saved the mermaid file inside the function
             # But wait, the function returns None implicitly if I didn't return paths.
            break

