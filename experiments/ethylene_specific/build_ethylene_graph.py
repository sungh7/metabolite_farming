#!/usr/bin/env python3
"""
Ethylene-Specific Graph Builder
- Metabolite 노드에 Log2FC, -log10(P-value) 피처 추가
- Enzyme 노드에 Proteomics FC 피처 추가 (있는 경우)
"""

import torch
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import sys
sys.path.append('/data/ethylene')

def build_ethylene_specific_graph():
    print("=" * 60)
    print("Building Ethylene-Specific Graph")
    print("=" * 60)
    
    # Paths
    base_graph_path = Path('/data/ethylene/data/processed/enhanced_bipartite_graph.pt')
    metabolomics_path = Path('/data/ethylene/results/Supplementary_Table_S3_All_Metabolites.csv')
    proteomics_path = Path('/data/ethylene/results/IFS_IFR_CHI_Evidence.csv')
    kegg_met_path = Path('/data/ethylene/data/kegg/metabolites.csv')
    output_dir = Path('/data/ethylene/experiments/ethylene_specific')
    
    # 1. Load base graph
    print("\n[1] Loading base graph...")
    data = torch.load(base_graph_path)
    
    # Get metabolite compound IDs
    met_list = data['Metabolite'].compound_ids
    n_metabolites = len(met_list)
    print(f"  Metabolites in graph: {n_metabolites}")
    
    # 2. Load metabolomics data
    print("\n[2] Loading metabolomics data...")
    met_df = pd.read_csv(metabolomics_path)
    print(f"  Total metabolites measured: {len(met_df)}")
    
    # KEGG mapping (compound_id -> log2fc, pvalue)
    kegg_met_df = pd.read_csv(kegg_met_path)
    kegg_to_fc = dict(zip(kegg_met_df['compound_id'], kegg_met_df['log2fc']))
    kegg_to_pval = dict(zip(kegg_met_df['compound_id'], kegg_met_df['pvalue']))
    print(f"  KEGG-mapped metabolites: {len(kegg_to_fc)}")
    
    # 3. Build metabolite features
    print("\n[3] Building metabolite features...")
    met_features = []
    matched = 0
    
    for met_id in met_list:
        if met_id in kegg_to_fc:
            fc = kegg_to_fc[met_id]
            pval = kegg_to_pval[met_id]
            # -log10(p-value), clipped to avoid inf
            neg_log_p = -np.log10(max(pval, 1e-50))
            # Significance flag
            is_sig = 1.0 if pval < 0.05 else 0.0
            met_features.append([fc, neg_log_p, is_sig])
            matched += 1
        else:
            # Non-experimental metabolites: 0 features (neutral)
            met_features.append([0.0, 0.0, 0.0])
    
    met_features = np.array(met_features, dtype=np.float32)
    print(f"  Matched experimental metabolites: {matched}/{n_metabolites}")
    print(f"  Feature dim: {met_features.shape[1]} (log2FC, -log10P, isSig)")
    
    # Normalize features (z-score for non-zero entries)
    for i in range(met_features.shape[1]):
        col = met_features[:, i]
        nonzero = col != 0
        if nonzero.sum() > 1:
            mean = col[nonzero].mean()
            std = col[nonzero].std() + 1e-8
            met_features[nonzero, i] = (col[nonzero] - mean) / std
    
    # 4. Combine with random embedding (64-dim)
    # Option A: Replace entirely with 3-dim experimental features
    # Option B: Concatenate experimental (3) + random (61) = 64 dim
    # We use Option B to maintain compatibility
    
    random_embed = torch.randn(n_metabolites, 61)
    exp_embed = torch.tensor(met_features, dtype=torch.float32)
    combined = torch.cat([exp_embed, random_embed], dim=1)  # 3 + 61 = 64
    
    print(f"\n[4] Final metabolite features: {combined.shape}")
    data['Metabolite'].x = combined
    
    # 5. Load proteomics for enzyme features (optional enhancement)
    print("\n[5] Loading proteomics data...")
    prot_df = pd.read_csv(proteomics_path)
    print(f"  Enzyme entries: {len(prot_df)}")
    
    # Gene ID -> Log2FC mapping
    gene_to_fc = {}
    for _, row in prot_df.iterrows():
        gene_id = row['Gene ID']
        fc = row['Log2 Fold Change']
        gene_to_fc[gene_id] = fc
    
    # For now, we don't have gene_id -> enzyme_idx mapping readily available
    # The enzyme features would require more complex mapping
    # Skip for this version (can be added later)
    print("  (Enzyme feature enhancement skipped - requires gene mapping)")
    
    # 6. Save
    output_path = output_dir / 'ethylene_specific_graph.pt'
    torch.save(data, output_path)
    print(f"\n[6] Saved to: {output_path}")
    
    # 7. Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Metabolites: {n_metabolites} (Experimental: {matched})")
    print(f"Metabolite Features: 3 experimental + 61 random = 64 dim")
    print(f"  - log2FC: 에틸렌 처리 시 변화량")
    print(f"  - -log10(P): 통계적 유의성")
    print(f"  - isSig: P<0.05 여부 (0/1)")
    
    return data

if __name__ == "__main__":
    build_ethylene_specific_graph()
