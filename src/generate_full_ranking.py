"""
Generate Full Enzyme Ranking for All Available Metabolites
Output ready for AlphaFold/Docking pipeline
"""

import torch
import pandas as pd
import os
from src.model import HGT, LinkPredictor


def generate_full_ranking():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Load graph and model
    data = torch.load('data/processed/enhanced_bipartite_graph.pt').to(device)
    
    model = HGT(data.metadata(), 64, 64, 64, 4, 2).to(device)
    model.load_state_dict(torch.load('data/models/refined_hgt.pth', map_location=device))
    predictor = LinkPredictor(64).to(device)
    
    model.eval()
    
    # Load enzyme mapping
    enzyme_df = pd.read_csv('data/processed/enzyme_string_mapping.csv')
    enzyme_idx_to_info = {
        row['enzyme_idx']: {
            'uniprot_id': row['uniprot_id'],
            'string_id': row['string_id']
        }
        for _, row in enzyme_df.iterrows()
    }
    
    # Get metabolite list
    met_list = data['Metabolite'].compound_ids
    print(f"Total metabolites: {len(met_list)}")
    
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        all_results = []
        
        # For each metabolite, get top 20 enzymes
        for met_idx, met_kegg in enumerate(met_list):
            # Score all enzymes
            num_enz = enz_emb.size(0)
            eval_src = torch.arange(num_enz, device=device)
            eval_dst = torch.full((num_enz,), met_idx, device=device)
            eval_edges = torch.stack([eval_src, eval_dst])
            
            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid().cpu().numpy()
            
            # Get top 20
            top_indices = scores.argsort()[::-1][:20]
            
            for rank, enz_idx in enumerate(top_indices, 1):
                info = enzyme_idx_to_info.get(enz_idx, {'uniprot_id': 'Unknown', 'string_id': 'Unknown'})
                all_results.append({
                    'metabolite_kegg': met_kegg,
                    'metabolite_idx': met_idx,
                    'rank': rank,
                    'enzyme_idx': enz_idx,
                    'uniprot_id': info['uniprot_id'],
                    'string_id': info['string_id'],
                    'score': float(scores[enz_idx])
                })
    
    # Save
    results_df = pd.DataFrame(all_results)
    results_df.to_csv('results/gnn/enhanced_enzyme_ranking.csv', index=False)
    print(f"Saved {len(results_df)} rankings ({len(met_list)} metabolites × 20 enzymes)")
    
    # Summary statistics
    print("\n=== Summary ===")
    print(f"Unique enzymes in top-20: {results_df['uniprot_id'].nunique()}")
    print(f"Score range: {results_df['score'].min():.3f} - {results_df['score'].max():.3f}")
    
    # Top enzymes overall (appearing most frequently in top-5)
    top5_df = results_df[results_df['rank'] <= 5]
    freq = top5_df['uniprot_id'].value_counts().head(10)
    print("\n=== Most Frequent Enzymes in Top-5 ===")
    for uniprot, count in freq.items():
        print(f"  {uniprot}: {count} metabolites")
    
    # Create docking-ready file
    docking_df = results_df[results_df['rank'] == 1].copy()
    docking_df = docking_df[['metabolite_kegg', 'uniprot_id', 'score']].drop_duplicates()
    docking_df.to_csv('results/gnn/top1_enzymes_for_docking.csv', index=False)
    print(f"\nSaved top-1 enzymes for docking: {len(docking_df)} pairs")


if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    generate_full_ranking()
