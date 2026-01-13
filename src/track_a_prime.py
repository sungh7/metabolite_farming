"""
Track A': Ethylene Source-Conditioned Ranking
Evaluates whether ethylene signaling nodes prioritize isoflavonoid module.
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import hypergeom
import sys
import os

sys.path.append(os.getcwd())
from src.model import HGT, LinkPredictor

# KEGG pathway modules
ISOFLAVONOID_PATHWAY = ['C02495', 'C00858', 'C10216']  # Daidzein, Formononetin, Daidzin
PHENYLPROPANOID_PATHWAY = ['C00062']  # L-Arginine (partial)
AMINO_ACID_PATHWAY = ['C00078', 'C00062']  # L-Tryptophan, L-Arginine

def load_trained_model(model_path, data, device):
    """Load trained HGT model."""
    model = HGT(data.metadata(), 64, 64, 64, num_heads=4, num_layers=2).to(device)
    if Path(model_path).exists():
        model.load_state_dict(torch.load(model_path, map_location=device))
    return model

def get_ethylene_signaling_nodes(data, n_nodes=20):
    """Simulate ethylene signaling node indices.
    In production: map ETR, CTR1, EIN2, EIN3, NAC TFs from STRING.
    """
    np.random.seed(42)
    num_enzymes = data['Enzyme'].num_nodes
    # Use top connected nodes as proxy for signaling hubs
    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
    degrees = torch.bincount(edge_index[0], minlength=num_enzymes).cpu().numpy()
    top_indices = np.argsort(-degrees)[:n_nodes]
    return torch.tensor(top_indices, device=data['Enzyme'].x.device)

def compute_metabolite_ranking(model, predictor, data, source_nodes, device):
    """Compute metabolite ranking from source enzyme nodes."""
    model.eval()
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']
        
        num_metabolites = met_emb.size(0)
        
        # Aggregate scores from source nodes to all metabolites
        scores = torch.zeros(num_metabolites, device=device)
        
        for enz_idx in source_nodes:
            # Score each metabolite
            enz_feat = enz_emb[enz_idx].unsqueeze(0).expand(num_metabolites, -1)
            score = (enz_feat * met_emb).sum(dim=1)
            scores += score
        
        # Average
        scores = scores / len(source_nodes)
        
        # Rank
        ranked_indices = torch.argsort(scores, descending=True)
        
        return ranked_indices.cpu().numpy(), scores.cpu().numpy()

def compute_module_enrichment(ranked_indices, met_list, module_compounds, total_metabolites, top_k=20):
    """Compute hypergeometric enrichment for a module in top-K ranked metabolites."""
    # Map compound IDs to indices
    module_indices = set()
    for i, met_id in enumerate(met_list):
        if met_id in module_compounds:
            module_indices.add(i)
    
    # Count hits in top-K
    top_k_set = set(ranked_indices[:top_k])
    hits_in_top_k = len(top_k_set & module_indices)
    
    # Hypergeometric test
    # M = total population (all metabolites)
    # n = success states in population (module size)
    # N = number of draws (top_k)
    # k = observed successes (hits)
    
    M = total_metabolites
    n = len(module_indices)
    N = top_k
    k = hits_in_top_k
    
    if n == 0:
        return None, 0, 0
    
    pvalue = 1 - hypergeom.cdf(k - 1, M, n, N)
    
    # Enrichment score
    expected = (n / M) * N
    enrichment = (k / expected) if expected > 0 else 0
    
    return pvalue, hits_in_top_k, n

def main():
    print("=" * 60)
    print("Track A': Ethylene Source-Conditioned Ranking")
    print("=" * 60)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    # Load data
    data = torch.load('data/processed/expanded_bipartite_graph.pt').to(device)
    met_list = data['Metabolite'].compound_ids
    
    # Train fresh model
    print("\nTraining HGT model...")
    model = HGT(data.metadata(), 64, 64, 64, num_heads=4, num_layers=2).to(device)
    predictor = LinkPredictor(64).to(device)
    optimizer = torch.optim.Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.01)
    
    for epoch in range(20):
        model.train()
        optimizer.zero_grad()
        x_dict = model(data.x_dict, data.edge_index_dict)
        pos_edges = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        num_pos = pos_edges.size(1)
        num_met = data['Metabolite'].num_nodes
        neg_dst = torch.randint(0, num_met, (num_pos,), device=device)
        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], pos_edges)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], torch.stack([pos_edges[0], neg_dst]))
        loss = -torch.log(torch.sigmoid(pos_out)+1e-15).mean() - torch.log(1-torch.sigmoid(neg_out)+1e-15).mean()
        loss.backward()
        optimizer.step()
    
    print("Training complete.\n")
    
    # Get source nodes
    ethylene_nodes = get_ethylene_signaling_nodes(data, n_nodes=20)
    random_nodes = torch.randint(0, data['Enzyme'].num_nodes, (20,), device=device)
    
    results = []
    
    for source_type, source_nodes in [('Ethylene', ethylene_nodes), ('Random', random_nodes)]:
        print(f"--- {source_type} Source Nodes ---")
        
        ranked_indices, scores = compute_metabolite_ranking(model, predictor, data, source_nodes, device)
        
        for module_name, module_compounds in [
            ('Isoflavonoid', ISOFLAVONOID_PATHWAY),
            ('Amino_acid', AMINO_ACID_PATHWAY)
        ]:
            pvalue, hits, module_size = compute_module_enrichment(
                ranked_indices, met_list, module_compounds, len(met_list), top_k=20
            )
            
            log_fdr = -np.log10(pvalue) if pvalue and pvalue > 0 else 0
            
            print(f"  {module_name}: {hits}/{module_size} in Top-20, p={pvalue:.4f}, -log10(p)={log_fdr:.2f}")
            
            results.append({
                'source': source_type,
                'module': module_name,
                'hits_in_top20': hits,
                'module_size': module_size,
                'pvalue': pvalue,
                'neg_log10_p': log_fdr
            })
    
    # Save
    df = pd.DataFrame(results)
    df.to_csv('results/track_a_prime.csv', index=False)
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(df.to_string(index=False))
    print(f"\nSaved to: results/track_a_prime.csv")
    
    # Interpretation
    print("\n--- Interpretation ---")
    eth_iso = df[(df['source']=='Ethylene') & (df['module']=='Isoflavonoid')]['neg_log10_p'].values[0]
    rand_iso = df[(df['source']=='Random') & (df['module']=='Isoflavonoid')]['neg_log10_p'].values[0]
    
    if eth_iso > rand_iso:
        print(f"Ethylene sources show higher isoflavonoid enrichment ({eth_iso:.2f}) vs Random ({rand_iso:.2f})")
    else:
        print(f"No significant difference between Ethylene and Random conditioning")

if __name__ == "__main__":
    main()
