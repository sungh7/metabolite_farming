import pandas as pd
import scipy.stats as stats
import urllib.request
import time
import os

def get_pathways(kegg_id):
    """
    Fetch pathways for a compound.
    Returns list of pathway IDs (e.g. map00940).
    """
    if pd.isna(kegg_id) or str(kegg_id).strip() == '':
        return []
        
    url = f"http://rest.kegg.jp/link/pathway/{kegg_id}"
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8').strip()
            if data:
                # Format: "cpd:C00123\tpath:map00940"
                paths = []
                for line in data.split('\n'):
                    p = line.split('\t')[1]
                    if 'map' in p:
                        paths.append(p.replace('path:', ''))
                return paths
    except:
        pass
    return []

def run_analysis(input_csv, output_csv):
    print(f"Loading {input_csv}...")
    df = pd.read_csv(input_csv)
    
    # Filter valid IDs
    valid_df = df.dropna(subset=['KEGG'])
    print(f"Valid KEGG IDs: {len(valid_df)} / {len(df)}")
    
    # Significant set
    sig_df = valid_df[valid_df['P_Value'] < 0.05]
    print(f"Significant: {len(sig_df)}")
    
    # Background set
    bg_df = valid_df
    
    # Map to pathways
    print("Fetching pathway info...")
    pathway_counts = {} # path -> {'sig': 0, 'bg': 0}
    pathway_names = {}
    
    # Helper to cache
    compound_paths = {}
    
    all_keggs = valid_df['KEGG'].unique()
    for i, kid in enumerate(all_keggs):
        paths = get_pathways(kid)
        compound_paths[kid] = paths
        if i % 5 == 0: time.sleep(0.2)
        
    # Count
    for _, row in valid_df.iterrows():
        kid = row['KEGG']
        paths = compound_paths.get(kid, [])
        is_sig = row['P_Value'] < 0.05
        
        for p in paths:
            if p not in pathway_counts:
                pathway_counts[p] = {'sig': 0, 'bg': 0}
            pathway_counts[p]['bg'] += 1
            if is_sig:
                pathway_counts[p]['sig'] += 1
                
    # Fisher Test
    results = []
    total_bg = len(valid_df)
    total_sig = len(sig_df)
    
    print("Calculating enrichment...")
    for pid, counts in pathway_counts.items():
        # Contingency Table
        #      Sig   NotSig
        # InPath   a     b
        # NotPath  c     d
        
        a = counts['sig']
        b = counts['bg'] - a
        c = total_sig - a
        d = (total_bg - total_sig) - b
        
        odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        
        results.append({
            'Pathway': pid,
            'Sig_Count': a,
            'Bg_Count': counts['bg'],
            'P_Value': pval,
            'Enrichment_Score': a / total_sig # Simple ratio
        })
        
    res_df = pd.DataFrame(results).sort_values('P_Value')
    
    # Add names (optional, hard to fetch efficiently)
    # Just save ID
    
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    res_df.to_csv(output_csv, index=False)
    print(f"Saved enrichment to {output_csv}")
    print(res_df.head(10))

if __name__ == "__main__":
    run_analysis('data/processed/mtbls531_differential.csv', 'results/table1_metabolomics_real.csv')
