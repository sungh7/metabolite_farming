import pandas as pd
import numpy as np
import scipy.stats as stats
import urllib.request
import time
import os

def map_chebi_to_kegg(chebi_id):
    """
    Maps ChEBI ID (e.g. CHEBI:138781) to KEGG Compound ID (e.g. C00123).
    """
    if pd.isna(chebi_id) or str(chebi_id).strip() == '':
        return None
    
    # Clean ID: CHEBI:138781 -> chebi:138781
    cid = str(chebi_id).lower().replace('chebi:', 'chebi:')
    if not cid.startswith('chebi:'):
        cid = 'chebi:' + cid.split(':')[-1]
        
    url = f"http://rest.kegg.jp/conv/compound/{cid}"
    
    try:
        # print(f"Querying {url}...")
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8').strip()
            if data:
                # Format: "chebi:138781\tcpd:C00000"
                parts = data.split('\t')
                if len(parts) > 1:
                    kegg_id = parts[1].replace('cpd:', '')
                    return kegg_id
    except Exception as e:
        # print(f"Failed to map {cid}: {e}")
        pass
        
    return None

def map_name_to_kegg(name):
    """
    Maps Metabolite Name to KEGG Compound ID using 'find' API.
    """
    if pd.isna(name) or str(name).strip() == '':
        return None
        
    # Clean name: remove "DL-", "(E)-", etc? KEGG might be sensitive.
    # Try exact first
    query = urllib.parse.quote(str(name).strip())
    url = f"http://rest.kegg.jp/find/compound/{query}"
    
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8').strip()
            if data:
                # Format: "cpd:C00123\tName"
                # Return first match
                first_line = data.split('\n')[0]
                kegg_id = first_line.split('\t')[0].replace('cpd:', '')
                return kegg_id
    except:
        pass
        
    return None

def process_data(maf_path, output_path):
    print(f"Loading {maf_path}...")
    df = pd.read_csv(maf_path, sep='\t')
    
    # Columns of interest
    control_cols = ['Control_1', 'Control_2', 'Control_3', 'Control_4']
    ethylene_cols = ['ethylene_1', 'ethylene_2', 'ethylene_3', 'ethylene_4']
    
    results = []
    
    print("Processing rows...")
    for idx, row in df.iterrows():
        chebi_id = row.get('database_identifier', None)
        name = row.get('metabolite_identification', 'Unknown')
        
        # Extract values and handle non-numeric
        try:
            con_vals = [float(row[c]) for c in control_cols if str(row[c]).replace('.','').isnumeric()]
            eth_vals = [float(row[c]) for c in ethylene_cols if str(row[c]).replace('.','').isnumeric()]
        except:
            continue
            
        if len(con_vals) < 2 or len(eth_vals) < 2:
            continue
            
        # Stats
        t_stat, p_val = stats.ttest_ind(eth_vals, con_vals, equal_var=False)
        mean_con = np.mean(con_vals)
        mean_eth = np.mean(eth_vals)
        
        # Log2FC
        if mean_con <= 0.001:
            log2fc = np.log2((mean_eth + 0.001) / 0.001)
        else:
            log2fc = np.log2((mean_eth + 0.001) / (mean_con + 0.001))
            
        results.append({
            'ChEBI': chebi_id,
            'Name': name,
            'Control_Mean': mean_con,
            'Ethylene_Mean': mean_eth,
            'Log2FC': log2fc,
            'P_Value': p_val if not np.isnan(p_val) else 1.0
        })
        
    res_df = pd.DataFrame(results)
    
    # Calculate False Discovery Rate (Bonferroni needed? or just map first)
    # Let's simple mapping first
    print(f"Mapping {len(res_df)} metabolites to KEGG IDs... (This may take a minute)")
    
    # Rate limit: KEGG doesn't like burst, but 80 is small.
    # We'll batch or just loop with sleep
    # Manual Map for difficult names
    MANUAL_MAP = {
        '13-OxoODE': 'C14765',
        'DL-Phenylalanine': 'C00503',
        "1,2-Di-(9Z,12Z,15Z-octadecatrienoyl)-3-(Galactosyl-alpha-1-6-Galactosyl-beta-1)-glycerol": "C06037", # DGDG
        "(Z)-4',6-Dihydroxyaurone 6-glucoside": "C10435", # Aureusidin 6-glucoside
        "1-(10Z-Heptadecenoyl)-sn-glycero-3-phospho-(1'-myo-inositol)": "C01177", # PI generic
        "(E)-3-Hexadecenoic acid": "C08320" # Palmitoleic acid (Proxy)
    }

    mapped_ids = []
    print("Mapping IDs...")
    for i, (cid, name) in enumerate(zip(res_df['ChEBI'], res_df['Name'])):
        kegg = None
        # Try manual
        if name in MANUAL_MAP:
            kegg = MANUAL_MAP[name]
            
        if not kegg:
            kegg = map_chebi_to_kegg(cid)
            
        if not kegg:
            # Fallback to Name
            kegg = map_name_to_kegg(name)
            
        mapped_ids.append(kegg)
        if i % 5 == 0: time.sleep(0.5)
        
    res_df['KEGG'] = mapped_ids
    
    # Filter significant
    # P < 0.05
    sig_df = res_df[res_df['P_Value'] < 0.05].copy()
    print(f"Found {len(sig_df)} significant metabolites (P<0.05)")
    
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    res_df.to_csv(output_path, index=False)
    print(f"Saved full results to {output_path}")

if __name__ == "__main__":
    process_data('data/experimental/maf.tsv', 'data/processed/mtbls531_differential.csv')
