import pandas as pd
import scipy.stats as stats
import torch
import os
import gzip

# Load Data
MET_PATH = 'data/processed/mtbls531_differential.csv'
PROT_PATH = 'data/processed/pxd006989_mapped.csv'
GRAPH_PATH = 'data/processed/bipartite_graph.pt' # from Phase 2
STRING_MAP_PATH = 'data/processed/string_mapping.csv' # Assuming generated in Phase 1? Or we can use the graph indices.

# We need to link Metabolites (KEGG) -> Enzymes (STRING ID).
# The bipartite graph stores edges indices.
# AND we need the node mapping (Index -> ID).
# Phase 2 `bipartite_builder.py` saved `data/processed/bipartite_graph.pt`.
# But did it save the ID maps? 
# Usually `HeteroData` stores features, not IDs unless we added them.
# The graph builder loaded 'metabolite_nodes.csv' and 'enzyme_nodes.csv' likely.

def load_graph_mapping():
    # We need to reconstruct or load the map.
    # In Phase 4/5, we used 'src/candidates.py'.
    # Let's verify what files define the nodes.
    # 'data/processed/enzyme_nodes.csv' (from STRING graph?)
    # 'data/processed/metabolite_nodes.csv' (Synthentic 200?)
    
    # Actually, in Phase 6, we used MTBLS531 (Real Data). The synthetic nodes are irrelevant.
    # BUT the "Edges" in our system are defined by KEGG/PlantCyc logic.
    # We need "KEGG Compound -> KEGG Enzyme (EC) -> Uniprot/STRING ID".
    
    # Since we don't have a perfect "Graph File" that persists this logic in CSV,
    # we might need to query the Reaction Logic again or reuse `bipartite_builder` logic.
    
    # Simplification:
    # 1. Start with Significant Metabolites (KEGG IDs).
    # 2. Get their "Linked EC numbers" (from KEGG API / reaction map).
    # 3. Get Proteins annotated with those ECs (from InterPro or STRING annotations?).
    # 4. Check overlap with Significant Proteins.
    
    # Or rely on the "Bipartite Graph" if it used Real KEGG logic.
    # Let's Assume Phase 2 did "Metabolite (Simulated ID) -> Enzyme".
    # Since we replaced the Metabolites, the old graph "Metabolite Nodes" are mismatch.
    # BUT the "Enzyme Nodes" are fixed (STRING).
    
    # So we need: Sig_Metabolite -> EC -> Gene/Protein.
    pass

import urllib.request
import time

def get_linked_genes(kegg_compound_id):
    """
    Simulate linkage via KEGG: Compound -> Reaction -> Enzyme -> Gene (Soybean)
    This is hard via API.
    Easier: Compound -> Pathway -> Genes in Pathway?
    Or Compound -> Enzyme (EC) -> Glyma Genes?
    """
    # For validation, let's use Pathway Overlap.
    # Metabolites enriched in "Secondary Metabolism".
    # Are Proteins also enriched in "Secondary Metabolism"?
    # If yes, CONCORDANCE.
    return []

def run_pathway_concordance():
    print("Loading datasets...")
    met_df = pd.read_csv(MET_PATH)
    prot_df = pd.read_csv(PROT_PATH)
    
    sig_mets = met_df[met_df['P_Value'] < 0.05]['KEGG'].dropna().unique()
    sig_prots = prot_df[(prot_df['P_Value'] < 0.05) & (abs(prot_df['Log2FC']) > 1)]['STRING_ID'].dropna().unique()
    
    print(f"Sig Mets: {len(sig_mets)}")
    print(f"Sig Prots: {len(sig_prots)}")
    
    # Mapping Proteins -> KEGG Pathways
    # We need STRING annotations or KEGG mapping for the Proteins.
    # `3847.protein.info` had comments like "GST superfamily".
    # Better: Use KEGG API `link pathway <gene_id>`
    # But we have STRING IDs.
    # Use mapped alias (Glyma) -> KEGG ID (gmx:...) -> Pathway using KEGG API.
    
    # Let's map a sample of Sig Prots to Pathways.
    # Extract Glyma ID from `pxd006989_mapped.csv` (it was in the input `Protein IDs`).
    
    prot_df_sig = prot_df[(prot_df['P_Value'] < 0.05) & (abs(prot_df['Log2FC']) > 1)].copy()
    
    # Extract first Glyma ID
    # Column 'Protein IDs' has "Glyma.01G..."
    valid_genes = []
    for raw_id in prot_df_sig['Protein IDs']:
        # "Glyma.01G000700.1.p" -> "gmx:GLYMA_01G000700" (KEGG format for Soybean)
        # KEGG uses "gmx:100816..." or "gmx:Glyma..."?
        # Checking KEGG Soybean organism code: 'gmx'
        # Let's guess format or use API find.
        if 'Glyma.' in str(raw_id):
            parts = str(raw_id).split(';')
            for p in parts:
                if 'Glyma.' in p:
                    # Glyma.01G000700... -> GLYMA_01G000700
                    # KEGG often uses `gmx:GLYMA_01G000700` (older) or Entrez.
                    # Let's just track the stripped ID.
                    base = p.split('.')[0] + "." + p.split('.')[1] # Glyma.01G000700
                    valid_genes.append(base)
                    break
    
    print(f"Valid Glyma Genes: {len(valid_genes)}")
    
    # Load STRING Annotations
    INFO_PATH = 'data/raw/3847.protein.info.v12.0.txt.gz'
    print(f"Loading {INFO_PATH}...")
    annotations = {}
    with gzip.open(INFO_PATH, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            # string_protein_id \t preferred_name \t protein_size \t annotation
            if len(parts) >= 4:
                annotations[parts[0]] = parts[3]
    
    cols = []
    matches = []
    
    targets = ['PAL', 'CHS', '4CL', 'HIDH', 'CHI', 'IFS', 'I3H', 'Chalcone', 'Phenylalanine', 'Isoflavone'] 
    
    for idx, row in prot_df_sig.iterrows():
        sid = row['STRING_ID']
        # Get annotation
        annot = annotations.get(sid, "")
        
        # Concat all info
        desc = str(row['Gene names']) + " " + str(row['Protein IDs']) + " " + annot
        
        # Simple string match
        found = []
        for t in targets:
            if t.upper() in desc.upper():
                found.append(t)
        
        if found:
            matches.append((found, row['Protein IDs'], row['Log2FC'], annot))
                
    print(f"\n--- Candidate Enzyme Matches in Proteomics ({len(matches)}) ---")
    for m in matches:
        print(m)
        
    # Output result
    with open('results/omics_cross_validation.txt', 'w') as f:
        f.write(f"Sig Metabolites: {len(sig_mets)}\n")
        f.write(f"Sig Proteins: {len(sig_prots)}\n")
        f.write(f"\nTarget Enzyme Overlap (Phenylpropanoid/Isoflavonoid):\n")
        for m in matches:
            f.write(f"Targets: {m[0]}, ID: {m[1]}, FC: {m[2]:.2f}, Annot: {m[3]}\n")
        
    # Also check STRING annotation for 'secondary metabolism'
    # ...
    
    # Output result
    with open('results/omics_cross_validation.txt', 'w') as f:
        f.write(f"Sig Metabolites: {len(sig_mets)}\n")
        f.write(f"Sig Proteins: {len(sig_prots)}\n")
        f.write(f"\nTarget Enzyme Overlap (Phenylpropanoid/Isoflavonoid):\n")
        for m in matches:
            f.write(f"{m}\n")
            
    print("Saved results/omics_cross_validation.txt")

if __name__ == "__main__":
    run_pathway_concordance()
