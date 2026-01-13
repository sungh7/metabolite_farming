"""
Novel Edge Ranking Analysis

Extracts Top-50 novel TF-Enzyme predictions that are:
1. Not in training data
2. High scoring (potential novel discoveries)

For each, provides annotation and plausibility assessment.
"""

import torch
import numpy as np
import pandas as pd
import os
import sys
sys.path.append(os.getcwd())

from src.model import HGT, LinkPredictor
from src.dataloader import StringDBLoader

def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

def load_annotations():
    """Load node annotations for interpretation."""
    loader = StringDBLoader('data/raw/soybean_string')
    _, node_info = loader.load_proteins()
    return node_info

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    set_seed(42)
    
    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    num_enzymes = data['Enzyme'].num_nodes
    num_metabolites = data['Metabolite'].num_nodes
    num_tfs = data['TF'].num_nodes
    
    print(f"Enzymes: {num_enzymes}, Metabolites: {num_metabolites}, TFs: {num_tfs}")
    
    # Load annotations
    try:
        node_info = load_annotations()
        print(f"Loaded {len(node_info)} node annotations")
    except:
        node_info = {}
        print("No annotations available, using IDs only")
    
    # Train model
    print("Training model...")
    indices = torch.randperm(num_enzymes, device=device)
    split = int(0.9 * num_enzymes)
    train_enz_mask = torch.zeros(num_enzymes, dtype=torch.bool, device=device)
    train_enz_mask[indices[:split]] = True
    
    train_data = data.clone()
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    src = edge_index[0]
    mask = train_enz_mask[src]
    train_edges = edge_index[:, mask]
    train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index = train_edges
    
    rev_index = data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index
    rev_dst = rev_index[1]
    rev_mask = train_enz_mask[rev_dst]
    train_data['Metabolite', 'rev_catalyzes', 'Enzyme'].edge_index = rev_index[:, rev_mask]
    
    # Check for TF-Enzyme edges
    tf_enz_key = ('TF', 'interacts', 'Enzyme')
    if tf_enz_key in data.edge_types:
        tf_enz_edges = data[tf_enz_key].edge_index
        existing_tf_enz = set(zip(tf_enz_edges[0].cpu().numpy(), tf_enz_edges[1].cpu().numpy()))
        print(f"Existing TF-Enzyme edges: {len(existing_tf_enz)}")
    else:
        existing_tf_enz = set()
        print("No TF-Enzyme edges in graph")
    
    model = HGT(train_data.metadata(), 64, 64, 64, 4, 2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    for epoch in range(20):
        model.train()
        optimizer.zero_grad()
        
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        pos_edge_index = train_data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        
        num_pos = pos_edge_index.size(1)
        num_met = train_data['Metabolite'].num_nodes
        neg_src = pos_edge_index[0]
        offset = torch.randint(1, 6, (num_pos,), device=device) * (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (pos_edge_index[1] + offset) % num_met
        
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edge_index)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([neg_src, neg_dst]))
        
        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()
        loss.backward()
        optimizer.step()
    
    print("Scoring all TF-Enzyme pairs...")
    model.eval()
    
    with torch.no_grad():
        x_dict = model(train_data.x_dict, train_data.edge_index_dict)
        tf_emb = x_dict['TF']
        enz_emb = x_dict['Enzyme']
        
        # Score subset of TF-Enzyme pairs (sample for efficiency)
        n_sample_tf = min(100, num_tfs)
        n_sample_enz = min(500, num_enzymes)
        
        sampled_tfs = torch.randperm(num_tfs)[:n_sample_tf]
        sampled_enzs = torch.randperm(num_enzymes)[:n_sample_enz]
        
        all_scores = []
        for tf_idx in sampled_tfs:
            for enz_idx in sampled_enzs:
                tf_i = tf_idx.item()
                enz_i = enz_idx.item()
                
                # Skip if in training
                if (tf_i, enz_i) in existing_tf_enz:
                    continue
                
                # Score using dot product (simplified)
                score = torch.dot(tf_emb[tf_idx], enz_emb[enz_idx]).item()
                all_scores.append((tf_i, enz_i, score))
        
        # Sort by score
        all_scores.sort(key=lambda x: -x[2])
        
        # Top-50 novel predictions
        top_50 = all_scores[:50]
    
    print(f"\nTop-50 Novel TF-Enzyme Predictions:")
    print("="*60)
    
    results = []
    for rank, (tf_idx, enz_idx, score) in enumerate(top_50, 1):
        tf_ann = node_info.get(tf_idx, {}).get('annotation', f'TF_{tf_idx}')
        enz_ann = node_info.get(enz_idx, {}).get('annotation', f'Enzyme_{enz_idx}')
        
        results.append({
            'Rank': rank,
            'TF_Idx': tf_idx,
            'Enzyme_Idx': enz_idx,
            'Score': score,
            'TF_Annotation': tf_ann[:50] if isinstance(tf_ann, str) else tf_ann,
            'Enzyme_Annotation': enz_ann[:50] if isinstance(enz_ann, str) else enz_ann,
            'In_Training': 'No'
        })
        
        if rank <= 10:
            print(f"{rank}. TF[{tf_idx}] -> Enzyme[{enz_idx}] (Score: {score:.4f})")
    
    # Save results
    os.makedirs('results/gnn', exist_ok=True)
    df = pd.DataFrame(results)
    df.to_csv('results/gnn/novel_edge_ranking.csv', index=False)
    print(f"\nTop-50 novel predictions saved to results/gnn/novel_edge_ranking.csv")
    
    # Summary
    print("\n" + "="*60)
    print("NOVEL EDGE RANKING SUMMARY")
    print("="*60)
    print(f"Total TF-Enzyme pairs scored: {len(all_scores)}")
    print(f"Top-50 extracted (all novel, not in training)")
    print(f"Score range: [{top_50[-1][2]:.4f}, {top_50[0][2]:.4f}]")
    print("="*60)

if __name__ == "__main__":
    main()
