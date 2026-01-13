"""
Corrected Coverage Analysis
Verifies compound→reaction existence via KEGG API directly.
"""

import requests
import time
import pandas as pd
from pathlib import Path

KEGG_DELAY = 0.35

def get_compound_reactions(compound_id: str) -> list:
    """Get reactions for a compound directly from KEGG."""
    url = f"https://rest.kegg.jp/link/reaction/{compound_id}"
    try:
        response = requests.get(url, timeout=10)
        time.sleep(KEGG_DELAY)
        if response.status_code == 200 and response.text.strip():
            reactions = []
            for line in response.text.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    rxn = parts[1].replace('rn:', '')
                    reactions.append(rxn)
            return reactions
        return []
    except:
        return []

def get_reaction_enzymes(reaction_id: str) -> list:
    """Get EC numbers for a reaction."""
    url = f"https://rest.kegg.jp/link/enzyme/{reaction_id}"
    try:
        response = requests.get(url, timeout=10)
        time.sleep(KEGG_DELAY)
        if response.status_code == 200 and response.text.strip():
            enzymes = []
            for line in response.text.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    ec = parts[1].replace('ec:', '')
                    enzymes.append(ec)
            return enzymes
        return []
    except:
        return []

def main():
    print("=" * 60)
    print("Corrected Coverage Analysis")
    print("=" * 60)
    
    kegg_dir = Path("data/kegg")
    
    # Load MTBLS531 metabolites
    mtbls = pd.read_csv(kegg_dir / "metabolites.csv")
    
    # Load gmx EC list
    gene_ec = pd.read_csv(kegg_dir / "gene_ec_mapping.tsv", sep='\t')
    gmx_ecs = set(gene_ec['ec'].unique())
    print(f"gmx ECs available: {len(gmx_ecs)}")
    
    # Load covered compounds from full fetch
    full_edges = pd.read_csv(kegg_dir / "full_enzyme_metabolite_edges.tsv", sep='\t')
    covered_in_graph = set(full_edges['metabolite_id'].unique())
    
    results = []
    
    for _, row in mtbls.iterrows():
        compound_id = row['compound_id']
        name = row['name']
        
        # Get reactions from KEGG
        reactions = get_compound_reactions(compound_id)
        
        # Check if any reaction EC is in gmx
        linked_ecs = []
        for rxn in reactions[:5]:  # Limit for speed
            ecs = get_reaction_enzymes(rxn)
            linked_ecs.extend(ecs)
        
        gmx_matched_ecs = set(linked_ecs) & gmx_ecs
        
        # Classification
        if compound_id in covered_in_graph:
            status = "COVERED"
        elif len(reactions) > 0 and len(gmx_matched_ecs) > 0:
            status = "REACTION_EXISTS_GMX_EC_FOUND"
        elif len(reactions) > 0:
            status = "REACTION_EXISTS_NO_GMX_EC"
        else:
            status = "NO_REACTION"
        
        results.append({
            'compound_id': compound_id,
            'name': name[:40],
            'n_reactions': len(reactions),
            'reactions': ';'.join(reactions[:3]),
            'n_ecs_total': len(linked_ecs),
            'n_ecs_in_gmx': len(gmx_matched_ecs),
            'gmx_ecs': ';'.join(list(gmx_matched_ecs)[:3]),
            'status': status
        })
        
        print(f"{compound_id}: {len(reactions)} rxns, {len(gmx_matched_ecs)} gmx ECs -> {status}")
    
    # Save
    df = pd.DataFrame(results)
    df.to_csv(kegg_dir / "corrected_coverage_analysis.csv", index=False)
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary by Status")
    print("=" * 60)
    for status in df['status'].unique():
        count = len(df[df['status'] == status])
        print(f"  {status}: {count}")
    
    # Detail on uncovered
    uncovered = df[df['status'] != 'COVERED']
    print(f"\nUncovered compounds: {len(uncovered)}")
    for _, row in uncovered.iterrows():
        print(f"  {row['compound_id']}: {row['status']} ({row['n_reactions']} rxns, {row['n_ecs_in_gmx']} gmx ECs)")

if __name__ == "__main__":
    main()
