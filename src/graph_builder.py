import torch
from torch_geometric.data import HeteroData
from src.dataloader import StringDBLoader
import os

def identify_node_type(annotation, preferred_name):
    """
    Heuristic to assign node type based on annotation/name.
    """
    ann = str(annotation).lower()
    name = str(preferred_name).lower()
    
    # Priority 1: Ethylene Signaling Pathway components
    if any(k in name for k in ['etr1', 'etr2', 'ers1', 'ein2', 'ein3', 'ctr1', 'ebf1', 'ebf2', 'ran1']) or \
       any(k in ann for k in ['ethylene receptor', 'ethylene response', 'ethylene-insensitive', 'constitutive triple response']):
        return 'Signaling'
        
    # Priority 2: Key Isoflavonoid Enzymes
    if any(k in name for k in ['pal', 'c4h', '4cl', 'chs', 'chi', 'ifs', 'hid']) or \
       any(k in ann for k in ['phenylalanine ammonia-lyase', 'chalcone synthase', 'chalcone isomerase', 'isoflavone synthase']):
        return 'Enzyme'

    # Priority 3: Transcription Factors
    if 'transcription factor' in ann or 'zinc finger' in ann or 'myb' in ann or 'wrky' in ann or 'bhlh' in ann or 'erf' in ann or 'dna-binding' in ann:
        return 'TF'
        
    # Priority 4: General Enzymes (Broad category)
    if 'synthase' in ann or 'kinase' in ann or 'transferase' in ann or 'reductase' in ann:
        return 'Enzyme'
        
    return 'Protein' # Default category

def build_graph(threshold=700, strict=False):
    loader = StringDBLoader()
    protein_map = loader.load_protein_info() # string_id -> name
    
    # Get all nodes first to determine types
    # Loader.node_to_idx has the mapping from string_id to 0..N-1
    # We need to reconstruct the list of string_ids in order
    sorted_string_ids = [None] * len(loader.node_to_idx)
    for sid, idx in loader.node_to_idx.items():
        sorted_string_ids[idx] = sid
        
    # 1. Determine Type for each node
    # Global Index -> (Type, Local Index)
    node_types = []
    type_counts = {'Signaling': 0, 'TF': 0, 'Enzyme': 0, 'Protein': 0}
    global_to_local = {} # global_idx -> (type_str, local_idx)
    
    # ... (Node classification logic remains same) ...
    
    # Re-reading annotation for classification
    # This is inefficient but one-time cost.
    import gzip
    node_annotations = {}
    with gzip.open(loader.protein_info_path, 'rt') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                sid = parts[0]
                ann = parts[3]
                node_annotations[sid] = ann
            elif len(parts) >= 2:
                 node_annotations[parts[0]] = ""

    # Assign types
    for global_idx, sid in enumerate(sorted_string_ids):
        # Only process if global_idx is filled (safety)
        if sid is None: continue
        
        name = protein_map.get(sid, "")
        ann = node_annotations.get(sid, "")
        
        ntype = identify_node_type(ann, name)
        type_counts[ntype] += 1
        
        # Current count is the local index
        local_idx = type_counts[ntype] - 1
        global_to_local[global_idx] = (ntype, local_idx)
        
    print(f"Node Counts: {type_counts}")
    
    # 2. Parse Edges and Remap (Pass Strict Flag)
    edges = loader.load_interactions(threshold=threshold, strict=strict)
    
    # Store edges: edges_dict[(src_type, dst_type)] = [[src_locals], [dst_locals]]
    edges_dict = {}
    
    print("Building Hetero Graph Edges...")
    for src_global, dst_global in edges:
        if src_global not in global_to_local or dst_global not in global_to_local: continue
        
        src_type, src_local = global_to_local[src_global]
        dst_type, dst_local = global_to_local[dst_global]
        
        edge_type = (src_type, 'interacts', dst_type)
        
        if edge_type not in edges_dict:
            edges_dict[edge_type] = [[], []]
        
        edges_dict[edge_type][0].append(src_local)
        edges_dict[edge_type][1].append(dst_local)
        
    # 3. Create HeteroData
    data = HeteroData()
    data.num_nodes_dict = type_counts  # Metadata
    
    # Add Nodes
    for ntype, count in type_counts.items():
        if count > 0:
            # Initialize with random embeddings (Trainable)
            # Embedding Dim = 64
            # Random initialization allows learning structural identity
            data[ntype].x = torch.randn(count, 64)
            data[ntype].num_nodes = count
            
    # Add Edges
    for (src_type, rel, dst_type), (src_list, dst_list) in edges_dict.items():
        edge_index = torch.tensor([src_list, dst_list], dtype=torch.long)
        data[src_type, rel, dst_type].edge_index = edge_index
    
    # 4. Save Enzyme Mapping (string_id -> local enzyme index)
    import pandas as pd
    enzyme_mapping = []
    for global_idx, (ntype, local_idx) in global_to_local.items():
        if ntype == 'Enzyme':
            sid = sorted_string_ids[global_idx]
            # Extract UniProt ID from STRING ID (e.g., "3847.A0A075W8S1" -> "A0A075W8S1")
            uniprot = sid.split('.')[-1] if '.' in sid else sid
            enzyme_mapping.append({
                'enzyme_idx': local_idx,
                'string_id': sid,
                'uniprot_id': uniprot
            })
    pd.DataFrame(enzyme_mapping).to_csv('data/processed/enzyme_string_mapping.csv', index=False)
    print(f"Saved enzyme mapping: {len(enzyme_mapping)} enzymes")
        
    return data

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--strict', action='store_true', help='Remove text-mining evidence')
    parser.add_argument('--threshold', type=int, default=700)
    args = parser.parse_args()
    
    data = build_graph(threshold=args.threshold, strict=args.strict)
    print("\nGraph Summary:")
    print(data)
    
    # Save
    os.makedirs('data/processed', exist_ok=True)
    filename = 'data/processed/strict_graph.pt' if args.strict else 'data/processed/graph.pt'
    torch.save(data, filename)
    print(f"Graph saved to {filename}")
