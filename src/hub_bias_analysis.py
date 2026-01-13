"""
Hub Bias Verification for Novel Edge Ranking

Analyzes whether top-ranked predictions are dominated by high-degree nodes.
Reports:
1. Enzyme degree distribution in Top-50
2. Degree percentile of top candidates
3. Degree-controlled statistics
"""

import torch
import pandas as pd
import numpy as np
import os
import sys
sys.path.append(os.getcwd())

def analyze_hub_bias():
    """Analyze hub bias in novel edge predictions."""
    
    # Load data
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    
    # Calculate node degrees
    num_enzymes = data['Enzyme'].num_nodes
    num_tfs = data['TF'].num_nodes
    
    # Enzyme degrees (from all edge types)
    enz_degree = torch.zeros(num_enzymes)
    
    # TF-Enzyme edges
    if ('TF', 'interacts', 'Enzyme') in data.edge_types:
        tf_enz_edges = data['TF', 'interacts', 'Enzyme'].edge_index
        for enz_idx in tf_enz_edges[1]:
            enz_degree[enz_idx] += 1
    
    # Enzyme-Enzyme edges
    if ('Enzyme', 'interacts', 'Enzyme') in data.edge_types:
        enz_enz_edges = data['Enzyme', 'interacts', 'Enzyme'].edge_index
        for enz_idx in enz_enz_edges[0]:
            enz_degree[enz_idx] += 1
        for enz_idx in enz_enz_edges[1]:
            enz_degree[enz_idx] += 1
    
    # Enzyme-Metabolite edges
    if ('Enzyme', 'catalyzes', 'Metabolite') in data.edge_types:
        enz_met_edges = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index
        for enz_idx in enz_met_edges[0]:
            enz_degree[enz_idx] += 1
    
    # Protein-Enzyme edges
    if ('Protein', 'interacts', 'Enzyme') in data.edge_types:
        prot_enz_edges = data['Protein', 'interacts', 'Enzyme'].edge_index
        for enz_idx in prot_enz_edges[1]:
            enz_degree[enz_idx] += 1
    
    enz_degree = enz_degree.numpy()
    
    # TF degrees
    tf_degree = torch.zeros(num_tfs)
    
    if ('TF', 'interacts', 'Enzyme') in data.edge_types:
        tf_enz_edges = data['TF', 'interacts', 'Enzyme'].edge_index
        for tf_idx in tf_enz_edges[0]:
            tf_degree[tf_idx] += 1
    
    if ('TF', 'interacts', 'TF') in data.edge_types:
        tf_tf_edges = data['TF', 'interacts', 'TF'].edge_index
        for tf_idx in tf_tf_edges[0]:
            tf_degree[tf_idx] += 1
        for tf_idx in tf_tf_edges[1]:
            tf_degree[tf_idx] += 1
    
    if ('TF', 'interacts', 'Protein') in data.edge_types:
        tf_prot_edges = data['TF', 'interacts', 'Protein'].edge_index
        for tf_idx in tf_prot_edges[0]:
            tf_degree[tf_idx] += 1
    
    tf_degree = tf_degree.numpy()
    
    print("="*60)
    print("NODE DEGREE STATISTICS")
    print("="*60)
    print(f"Enzyme degrees: mean={enz_degree.mean():.2f}, max={enz_degree.max():.0f}, median={np.median(enz_degree):.0f}")
    print(f"TF degrees: mean={tf_degree.mean():.2f}, max={tf_degree.max():.0f}, median={np.median(tf_degree):.0f}")
    
    # Load novel edge ranking
    novel_df = pd.read_csv('results/gnn/novel_edge_ranking.csv')
    
    # Analyze Top-50
    print("\n" + "="*60)
    print("TOP-50 HUB BIAS ANALYSIS")
    print("="*60)
    
    # Count occurrences
    enz_counts = novel_df['Enzyme_Idx'].value_counts()
    tf_counts = novel_df['TF_Idx'].value_counts()
    
    print(f"\nTop-50 Enzyme distribution:")
    print(f"  Unique enzymes: {len(enz_counts)} / 50")
    print(f"  Most frequent enzyme: Enzyme_{enz_counts.index[0]} ({enz_counts.iloc[0]} times)")
    
    print(f"\nTop-50 TF distribution:")
    print(f"  Unique TFs: {len(tf_counts)} / 50")
    print(f"  Most frequent TF: TF_{tf_counts.index[0]} ({tf_counts.iloc[0]} times)")
    
    # Degree percentiles for top candidates
    print("\n" + "="*60)
    print("DEGREE PERCENTILES OF TOP CANDIDATES")
    print("="*60)
    
    results = []
    
    for rank in range(min(10, len(novel_df))):
        row = novel_df.iloc[rank]
        enz_idx = int(row['Enzyme_Idx'])
        tf_idx = int(row['TF_Idx'])
        
        enz_deg = enz_degree[enz_idx]
        tf_deg = tf_degree[tf_idx]
        
        enz_pct = (enz_degree <= enz_deg).mean() * 100
        tf_pct = (tf_degree <= tf_deg).mean() * 100
        
        results.append({
            'Rank': rank + 1,
            'TF_Idx': tf_idx,
            'Enzyme_Idx': enz_idx,
            'TF_Degree': int(tf_deg),
            'TF_Percentile': enz_pct,
            'Enz_Degree': int(enz_deg),
            'Enz_Percentile': enz_pct
        })
        
        print(f"Rank {rank+1}: TF[{tf_idx}] (deg={int(tf_deg)}, pct={tf_pct:.1f}%) -> Enz[{enz_idx}] (deg={int(enz_deg)}, pct={enz_pct:.1f}%)")
    
    # Hub bias assessment
    top_enz_idx = enz_counts.index[0]
    top_enz_deg = enz_degree[top_enz_idx]
    top_enz_pct = (enz_degree <= top_enz_deg).mean() * 100
    
    top_tf_idx = tf_counts.index[0]
    top_tf_deg = tf_degree[top_tf_idx]
    top_tf_pct = (tf_degree <= top_tf_deg).mean() * 100
    
    print("\n" + "="*60)
    print("HUB BIAS SUMMARY")
    print("="*60)
    print(f"Most frequent Enzyme ({enz_counts.iloc[0]}/50 predictions): Enzyme_{top_enz_idx}")
    print(f"  - Degree: {int(top_enz_deg)}")
    print(f"  - Percentile: {top_enz_pct:.1f}% (higher = more connected)")
    
    print(f"\nMost frequent TF ({tf_counts.iloc[0]}/50 predictions): TF_{top_tf_idx}")
    print(f"  - Degree: {int(top_tf_deg)}")
    print(f"  - Percentile: {top_tf_pct:.1f}%")
    
    # Interpretation
    if top_enz_pct > 90:
        print("\n⚠️ WARNING: Top enzyme is in top 10% by degree - potential hub bias")
        interpretation = "hub_bias_possible"
    elif top_enz_pct > 75:
        print("\n⚠️ CAUTION: Top enzyme is moderately connected - mixed signal")
        interpretation = "moderate_degree"
    else:
        print("\n✓ Top enzyme is NOT a high-degree hub - predictions may be meaningful")
        interpretation = "low_degree_good"
    
    # Save results
    os.makedirs('results/gnn', exist_ok=True)
    with open('results/gnn/hub_bias_analysis.txt', 'w') as f:
        f.write("Hub Bias Analysis for Novel Edge Ranking\n")
        f.write("="*50 + "\n\n")
        f.write(f"Top-50 Enzyme distribution: {len(enz_counts)} unique / 50 total\n")
        f.write(f"Top-50 TF distribution: {len(tf_counts)} unique / 50 total\n\n")
        f.write(f"Most frequent Enzyme: Enzyme_{top_enz_idx}\n")
        f.write(f"  - Occurrences in Top-50: {enz_counts.iloc[0]}\n")
        f.write(f"  - Degree: {int(top_enz_deg)}\n")
        f.write(f"  - Degree Percentile: {top_enz_pct:.1f}%\n\n")
        f.write(f"Most frequent TF: TF_{top_tf_idx}\n")
        f.write(f"  - Occurrences in Top-50: {tf_counts.iloc[0]}\n")
        f.write(f"  - Degree: {int(top_tf_deg)}\n")
        f.write(f"  - Degree Percentile: {top_tf_pct:.1f}%\n\n")
        f.write(f"Interpretation: {interpretation}\n")
    
    print("\nResults saved to results/gnn/hub_bias_analysis.txt")
    
    return results

if __name__ == "__main__":
    analyze_hub_bias()
