"""
Generate Top-20 Enzyme Ranking with Real Gene IDs
"""

import torch
import pandas as pd
import os
from src.model import HGT, LinkPredictor


def generate_proper_ranking():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Load graph and model
    data = torch.load('data/processed/enhanced_bipartite_graph.pt').to(device)
    
    model = HGT(data.metadata(), 64, 64, 64, 4, 2).to(device)
    model.load_state_dict(torch.load('data/models/refined_hgt.pth', map_location=device))
    predictor = LinkPredictor(64).to(device)
    
    model.eval()
    
    # Load enzyme mapping
    enzyme_df = pd.read_csv('data/processed/enzyme_string_mapping.csv')
    enzyme_idx_to_uniprot = dict(zip(enzyme_df['enzyme_idx'], enzyme_df['uniprot_id']))
    
    # Get metabolite list (from KEGG)
    met_list = data['Metabolite'].compound_ids
    
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        results = []
        
        # For each target metabolite (e.g., isoflavonoid pathway)
        target_mets = {
            'C00858': 'Formononetin',  
            'C00509': 'Daidzein',
            'C00786': 'Genistein',
            'C00389': 'Quercetin',
            'C00073': 'L-Methionine'  # Control
        }
        
        for met_kegg, met_name in target_mets.items():
            if met_kegg in met_list:
                met_idx = met_list.index(met_kegg)
            else:
                print(f"Metabolite {met_kegg} not in graph")
                continue
            
            # Score all enzymes for this metabolite
            num_enz = enz_emb.size(0)
            eval_src = torch.arange(num_enz, device=device)
            eval_dst = torch.full((num_enz,), met_idx, device=device)
            eval_edges = torch.stack([eval_src, eval_dst])
            
            scores = predictor(enz_emb, met_emb, eval_edges).sigmoid().cpu().numpy()
            
            # Get top 20
            top_indices = scores.argsort()[::-1][:20]
            
            for rank, enz_idx in enumerate(top_indices, 1):
                uniprot = enzyme_idx_to_uniprot.get(enz_idx, "Unknown")
                results.append({
                    'metabolite_kegg': met_kegg,
                    'metabolite_name': met_name,
                    'rank': rank,
                    'enzyme_idx': enz_idx,
                    'uniprot_id': uniprot,
                    'score': float(scores[enz_idx])
                })
    
    # Save
    results_df = pd.DataFrame(results)
    results_df.to_csv('results/gnn/proper_enzyme_ranking.csv', index=False)
    print(f"Saved {len(results_df)} rankings")
    
    # Print summary
    print("\n=== Top-5 Enzymes per Metabolite ===")
    for met_kegg in target_mets.keys():
        subset = results_df[results_df['metabolite_kegg'] == met_kegg].head(5)
        if not subset.empty:
            print(f"\n{met_kegg} ({target_mets[met_kegg]}):")
            for _, row in subset.iterrows():
                print(f"  {row['rank']}. {row['uniprot_id']} (score: {row['score']:.3f})")


if __name__ == "__main__":
    import sys
    sys.path.append(os.getcwd())
    generate_proper_ranking()
