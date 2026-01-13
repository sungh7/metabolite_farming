import torch
import torch_geometric.transforms as T
from torch_geometric.data import HeteroData
from src.model import HGT, LinkPredictor
from src.dataloader import StringDBLoader
import pandas as pd
import os


def predict_candidates(top_k=20):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # 1. Load Data
    print("Loading graph...")
    data = torch.load('data/processed/graph.pt')
    
    # Need metadata for model init
    data = T.ToUndirected()(data)
    data = T.AddSelfLoops()(data)
    metadata = data.metadata()
    data = data.to(device)
    
    # 2. Load Model
    print("Loading model...")
    model = HGT(metadata, 64, 64, 64, 4, 2).to(device)
    model.load_state_dict(torch.load('data/models/hgt_model.pth'))
    model.eval()
    
    predictor = LinkPredictor(64).to(device)
    predictor.load_state_dict(torch.load('data/models/predictor.pth'))
    predictor.eval()
    
    # 3. Get Embeddings
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        
    # 4. Generate Candidates
    # Target: TF -> Enzyme
    src_type = 'TF'
    dst_type = 'Enzyme'
    
    print(f"Generating candidates for {src_type} -> {dst_type}...")
    
    # Get all possible pairs?
    # TF count: ~2800, Enzyme count: ~3400 -> Total 9.5M pairs. Feasible on GPU.
    num_src = data[src_type].num_nodes
    num_dst = data[dst_type].num_nodes
    
    # Create indices
    src_indices = torch.arange(num_src, device=device).repeat_interleave(num_dst)
    dst_indices = torch.arange(num_dst, device=device).repeat(num_src)
    edge_label_index = torch.stack([src_indices, dst_indices], dim=0)
    
    # Predict
    with torch.no_grad():
        # Batching might be needed if OOM
        # 9.5M pairs * 64 dim * 4 bytes ~ 2.4 GB. Should fit.
        scores = predictor(x_dict[src_type], x_dict[dst_type], edge_label_index)
        scores = torch.sigmoid(scores)
        
    # 5. Filter out existing edges
    print("Filtering existing interactions...")
    # Get existing edges
    if ('TF', 'interacts', 'Enzyme') in data.edge_types:
        existing_edges = data['TF', 'interacts', 'Enzyme'].edge_index
        # Convert to set of tuples for fast lookup? Or sparse mask.
        # Since we flattened the full matrix (src * dst), we can flatten existing indices
        # existing_flat_idx = src * num_dst + dst
        existing_flat_idx = existing_edges[0] * num_dst + existing_edges[1]
        
        # Create a mask
        mask = torch.ones(scores.size(0), dtype=torch.bool, device=device)
        mask[existing_flat_idx] = False
        
        # Apply mask
        scores = scores[mask]
        edge_label_index = edge_label_index[:, mask]
        
    # 6. Extract Top-K
    print(f"Extracting Top-{top_k}...")
    top_scores, top_indices = torch.topk(scores, top_k)
    
    top_src_idx = edge_label_index[0][top_indices].cpu().numpy()
    top_dst_idx = edge_label_index[1][top_indices].cpu().numpy()
    top_scores = top_scores.cpu().numpy()
    
    # 7. Map back to names
    # We need the StringDBLoader to get names, but GraphBuilder didn't save the id_to_name map efficiently.
    # However, GraphBuilder processed node_types in order of loader.node_to_idx?
    # Let's check graph_builder.py logic.
    # "sorted_string_ids = [None] * len(loader.node_to_idx)" ...
    # "for global_idx, sid in enumerate(sorted_string_ids): ... global_to_local[global_idx] = (ntype, local_idx)"
    # We need to reverse: (ntype, local_idx) -> global_idx -> string_id -> name
    
    # Re-run loader logic to reconstruct mapping? Or save it.
    # Saving it is better, but I didn't save it. I must reconstruction.
    loader = StringDBLoader()
    protein_map = loader.load_protein_info()
    
    sorted_string_ids = [None] * len(loader.node_to_idx)
    for sid, idx in loader.node_to_idx.items():
        sorted_string_ids[idx] = sid
        
    # Re-classify to get local-to-global mapping
    # This assumes deterministic classification (same as before)
    # We need the annotations again... this is heavy.
    # Alternative: GraphBuilder saved `data` which usually doesn't have names.
    # Let's verify `graph.pt`. It has no names.
    
    # Hack: We must replicate the exact logic of graph_builder.py to get the mapping right.
    # Let's import the logic or Copy-Paste the classification logic.
    
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
                 
    from src.graph_builder import identify_node_type
    
    local_to_global = {'TF': {}, 'Enzyme': {}, 'Signaling': {}, 'Protein': {}}
    type_counts_replay = {'Signaling': 0, 'TF': 0, 'Enzyme': 0, 'Protein': 0}

    for global_idx, sid in enumerate(sorted_string_ids):
        name = protein_map.get(sid, "")
        ann = node_annotations.get(sid, "")
        ntype = identify_node_type(ann, name)
        
        local_idx = type_counts_replay[ntype] # 0-based
        local_to_global[ntype][local_idx] = (sid, name, ann)
        type_counts_replay[ntype] += 1
        
    # Now map
    results = []
    for i in range(top_k):
        s_idx = top_src_idx[i]
        d_idx = top_dst_idx[i]
        score = top_scores[i]
        
        tf_info = local_to_global['TF'].get(s_idx, ("?", "?", "?"))
        enz_info = local_to_global['Enzyme'].get(d_idx, ("?", "?", "?"))
        
        results.append({
            'Rank': i+1,
            'TF_ID': tf_info[0],
            'TF_Name': tf_info[1],
            'Enzyme_ID': enz_info[0],
            'Enzyme_Name': enz_info[1],
            'Score': float(score),
            'TF_Desc': tf_info[2],
            'Enzyme_Desc': enz_info[2]
        })
        
    df = pd.DataFrame(results)
    print(df)
    os.makedirs('results', exist_ok=True)
    df.to_csv('results/candidates.csv', index=False)
    print("Saved to results/candidates.csv")

if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    predict_candidates()
