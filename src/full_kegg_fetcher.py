#!/usr/bin/env python3
"""
Full KEGG Reaction Fetcher
Processes ALL EC numbers from Glycine max to maximize edge coverage.
Includes progress saving for resumability.
"""

import pandas as pd
import time
import requests
import re
import json
from pathlib import Path
from collections import defaultdict

KEGG_DELAY = 0.35

CURRENCY_METABOLITES = {
    'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006',
    'C00008', 'C00009', 'C00010', 'C00011', 'C00013', 'C00014',
    'C00020', 'C00027', 'C00044', 'C00080'
}

def fetch_ec_reactions(ec_id: str) -> list:
    """Fetch reaction IDs for a given EC number."""
    url = f"https://rest.kegg.jp/link/reaction/{ec_id}"
    try:
        response = requests.get(url, timeout=10)
        time.sleep(KEGG_DELAY)
        reactions = []
        if response.status_code == 200 and response.text.strip():
            for line in response.text.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    rxn_id = parts[1].replace('rn:', '')
                    reactions.append(rxn_id)
        return reactions
    except:
        return []

def fetch_reaction_compounds(reaction_id: str) -> dict:
    """Fetch substrate/product metabolites for a reaction."""
    url = f"https://rest.kegg.jp/get/rn:{reaction_id}"
    try:
        response = requests.get(url, timeout=10)
        time.sleep(KEGG_DELAY)
        if response.status_code != 200:
            return None
        info = {'substrates': [], 'products': []}
        for line in response.text.split('\n'):
            if line.startswith('EQUATION'):
                eq = line.replace('EQUATION', '').strip()
                if '<=>' in eq:
                    left, right = eq.split('<=>')
                elif '=>' in eq:
                    left, right = eq.split('=>')
                elif '<=' in eq:
                    right, left = eq.split('<=')
                else:
                    left, right = eq, ''
                info['substrates'] = re.findall(r'C\d{5}', left)
                info['products'] = re.findall(r'C\d{5}', right)
                break
        return info
    except:
        return None

def main():
    print("=" * 60)
    print("Full KEGG Reaction Fetcher - All ECs")
    print("=" * 60)
    
    kegg_dir = Path("data/kegg")
    
    # Load gene-EC mapping (already fetched)
    gene_ec_path = kegg_dir / "gene_ec_mapping.tsv"
    if not gene_ec_path.exists():
        print("Error: gene_ec_mapping.tsv not found. Run kegg_reaction_fetcher.py first.")
        return
    
    gene_ec_df = pd.read_csv(gene_ec_path, sep='\t')
    all_ecs = sorted(gene_ec_df['ec'].unique())
    print(f"Total unique ECs to process: {len(all_ecs)}")
    
    # Progress tracking
    progress_file = kegg_dir / "full_fetch_progress.json"
    processed_ecs = set()
    reaction_cache = {}
    
    if progress_file.exists():
        with open(progress_file, 'r') as f:
            progress = json.load(f)
            processed_ecs = set(progress.get('processed_ecs', []))
            reaction_cache = progress.get('reaction_cache', {})
        print(f"Resuming from checkpoint: {len(processed_ecs)} ECs already processed")
    
    # Collect all edges
    all_edges = []
    all_metabolites = set()
    
    remaining_ecs = [ec for ec in all_ecs if ec not in processed_ecs]
    print(f"ECs remaining: {len(remaining_ecs)}")
    
    for i, ec in enumerate(remaining_ecs):
        if i % 50 == 0 and i > 0:
            print(f"Progress: {i}/{len(remaining_ecs)} ({len(all_edges)} edges so far)")
            # Save checkpoint
            with open(progress_file, 'w') as f:
                json.dump({
                    'processed_ecs': list(processed_ecs),
                    'reaction_cache': reaction_cache
                }, f)
        
        reactions = fetch_ec_reactions(ec)
        
        for rxn_id in reactions[:10]:  # Limit per EC
            if rxn_id in reaction_cache:
                rxn_info = reaction_cache[rxn_id]
            else:
                rxn_info = fetch_reaction_compounds(rxn_id)
                if rxn_info:
                    reaction_cache[rxn_id] = rxn_info
            
            if rxn_info is None:
                continue
            
            substrates = [s for s in rxn_info['substrates'] if s not in CURRENCY_METABOLITES]
            products = [p for p in rxn_info['products'] if p not in CURRENCY_METABOLITES]
            all_mets = set(substrates + products)
            
            for met in all_mets:
                all_edges.append({
                    'enzyme_ec': ec,
                    'metabolite_id': met,
                    'reaction_id': rxn_id,
                    'is_substrate': met in substrates,
                    'is_product': met in products
                })
                all_metabolites.add(met)
        
        processed_ecs.add(ec)
    
    # Final save
    with open(progress_file, 'w') as f:
        json.dump({
            'processed_ecs': list(processed_ecs),
            'reaction_cache': reaction_cache
        }, f)
    
    # Deduplicate and save
    edges_df = pd.DataFrame(all_edges)
    edges_df = edges_df.drop_duplicates(subset=['enzyme_ec', 'metabolite_id'])
    edges_df.to_csv(kegg_dir / "full_enzyme_metabolite_edges.tsv", sep='\t', index=False)
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total ECs processed: {len(processed_ecs)}")
    print(f"Unique reactions cached: {len(reaction_cache)}")
    print(f"Unique metabolites found: {len(all_metabolites)}")
    print(f"Total edges (deduplicated): {len(edges_df)}")
    
    if len(edges_df) > 0:
        print(f"Unique EC numbers in edges: {edges_df['enzyme_ec'].nunique()}")
    
    # Check overlap with MTBLS531
    mtbls_path = kegg_dir / "metabolites.csv"
    if mtbls_path.exists():
        mtbls = pd.read_csv(mtbls_path)
        target_mets = set(mtbls['compound_id'])
        covered = all_metabolites & target_mets
        print(f"\nMTBLS531 coverage: {len(covered)}/{len(target_mets)}")
        print(f"Covered: {sorted(covered)}")
    
    print(f"\nSaved to: {kegg_dir / 'full_enzyme_metabolite_edges.tsv'}")

if __name__ == "__main__":
    main()
