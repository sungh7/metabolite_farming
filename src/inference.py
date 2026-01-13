import torch
from src.model import HGT, LinkPredictor
from src.dataloader import StringDBLoader
from src.graph_builder import identify_node_type 
import argparse
import os
import gzip

def run_inference(graph_path, model_path, output_dir):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Loading Graph from {graph_path}...")
    data = torch.load(graph_path, map_location=device)
    
    # 1. Load Model
    print(f"Loading Model from {model_path}...")
    metadata = data.metadata()
    model = HGT(metadata, 64, 64, 64, num_heads=4, num_layers=2).to(device)
    predictor = LinkPredictor(64).to(device)
    
    # Load Weights
    checkpoint = torch.load(model_path, map_location=device)
    model.load_state_dict(checkpoint)
    # predictor has no params
    
    model.eval()
    
    # 2. Forward Pass
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        tf_emb = x_dict['TF'] # [Num_TFs, 64]
        enz_emb = x_dict['Enzyme'] # [Num_Enzymes, 64]

    # 3. Reconstruct ID Mapping & Filter TFs
    loader = StringDBLoader()
    protein_map = loader.load_protein_info() # sid -> name
    
    # Reconstruct sorted order matches graph_builder
    sorted_sids = [None] * len(loader.node_to_idx)
    for sid, idx in loader.node_to_idx.items():
        sorted_sids[idx] = sid
        
    # Load Annotations
    node_annotations = {}
    with gzip.open(loader.protein_info_path, 'rt') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                node_annotations[parts[0]] = parts[3] if len(parts) >= 4 else ""

    tf_metadata = [] # List of (local_idx, sid, name, annotation)
    enz_metadata = [] # List of (local_idx, sid, name, annotation)
    
    tf_count = 0
    enz_count = 0
    
    verified_tf_indices = [] # Local indices
    
    print("Reconstructing Mappings & Filtering TFs...")
    for sid in sorted_sids:
        if sid is None: continue
        name = protein_map.get(sid, "")
        ann = node_annotations.get(sid, "")
        ntype = identify_node_type(ann, name)
        
        if ntype == 'TF':
            # Check for Verified TFs
            n_lower = name.lower()
            is_verified = False
            if any(x in n_lower for x in ['wrky', 'myb', 'erf', 'bhlh', 'nac']):
                is_verified = True
                verified_tf_indices.append(tf_count)
            
            tf_metadata.append({
                'local_idx': tf_count,
                'sid': sid,
                'name': name,
                'ann': ann,
                'verified': is_verified
            })
            tf_count += 1
            
        elif ntype == 'Enzyme':
            enz_metadata.append({
                'local_idx': enz_count,
                'sid': sid,
                'name': name,
                'ann': ann
            })
            enz_count += 1
            
    print(f"Total TFs: {tf_count}, Verified TFs: {len(verified_tf_indices)}")
    
    if len(verified_tf_indices) == 0:
        print("No verified TFs found! Check criteria.")
        return

    # 4. Score Verified Pairs
    verified_tensor = torch.tensor(verified_tf_indices, device=device)
    subset_tf_emb = tf_emb[verified_tensor] # [Num_Verified, 64]
    
    # Score Matrix: [Num_Verified, Num_Enzymes]
    scores = (subset_tf_emb @ enz_emb.t())
    scores = torch.sigmoid(scores)
    
    # 5. Mask Existing Edges
    # We need to find which (tf_local, enz_local) edges exist
    existing_edges = set()
    if ('TF', 'interacts', 'Enzyme') in data.edge_index_dict:
        edge_index = data['TF', 'interacts', 'Enzyme'].edge_index.cpu()
        src, dst = edge_index[0].numpy(), edge_index[1].numpy()
        for s, d in zip(src, dst):
            existing_edges.add((s, d))
            
    # Masking loop
    # Create mask tensor
    mask = torch.ones_like(scores, dtype=torch.bool)
    
    # Mapping: Verified Index (0..V-1) -> Original TF Local Index
    # We iterate and mask
    for i, original_tf_idx in enumerate(verified_tf_indices):
        for j in range(enz_count):
            if (original_tf_idx, j) in existing_edges:
                mask[i, j] = False
                
    masked_scores = scores * mask.float()
    
    # 6. Rank Top 50
    # Flatten and TopK
    flat_scores = masked_scores.flatten()
    topk_vals, topk_indices = torch.topk(flat_scores, 50)
    
    results = []
    for rank, (score, idx) in enumerate(zip(topk_vals.tolist(), topk_indices.tolist())):
        # Decode indices
        verified_idx = idx // enz_count
        enz_idx = idx % enz_count
        
        tf_meta = tf_metadata[verified_tf_indices[verified_idx]]
        enz_meta = enz_metadata[enz_idx]
        
        results.append({
            'rank': rank + 1,
            'tf_id': tf_meta['sid'],
            'tf_name': tf_meta['name'],
            'enz_id': enz_meta['sid'],
            'enz_name': enz_meta['name'],
            'enz_ann': enz_meta['ann'],
            'score': score
        })
        
    # 7. Save
    os.makedirs(output_dir, exist_ok=True)
    out_tsv = os.path.join(output_dir, 'top_novel_pairs.tsv')
    with open(out_tsv, 'w') as f:
        f.write("Rank\tTF_ID\tTF_Name\tEnzyme_ID\tEnzyme_Name\tEnzyme_Ann\tScore\n")
        f.write("-" * 80 + "\n")
        for r in results:
            f.write(f"{r['rank']}\t{r['tf_id']}\t{r['tf_name']}\t{r['enz_id']}\t{r['enz_name']}\t{r['enz_ann']}\t{r['score']:.4f}\n")
            
    print(f"Saved Top-50 Verified Pairs to {out_tsv}")
    
    # Save Top-1 MD
    top1 = results[0]
    md_content = f"""# Case Study: Top-1 Novel Verified TF-Enzyme Pair

**Pair**: **{top1['tf_name']}** (TF) -- **{top1['enz_name']}** (Enzyme)
**Score**: {top1['score']:.4f} (Rank #1 among Verified TFs)
**TF ID**: `{top1['tf_id']}`
**Enzyme ID**: `{top1['enz_id']}`

## Context
- **TF**: {top1['tf_name']} (Verified Transcription Factor)
- **Enzyme**: {top1['enz_name']}
    - Annotation: {top1['enz_ann']}
- **Novelty**: Interaction absent from training graph (Strict Mode).

## Biological Hypothesis
This high-confidence link proposes a direct regulatory axis between the verified TF **{top1['tf_name']}** and the enzyme **{top1['enz_name']}**, potentially linking transcriptional control to metabolic output.
"""
    with open(os.path.join(output_dir, 'top1_novel_pair.md'), 'w') as f:
        f.write(md_content)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--graph', type=str, default='data/processed/strict_bipartite.pt')
    parser.add_argument('--model', type=str, default='data/models/refined_hgt_strict.pth')
    parser.add_argument('--output', type=str, default='results/case_study')
    args = parser.parse_args()
    
    run_inference(args.graph, args.model, args.output)
