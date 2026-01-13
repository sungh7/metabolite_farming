import pandas as pd
import gzip
import os

def load_aliases(alias_path):
    print(f"Loading {alias_path}...")
    # Format: string_protein_id alias source
    # We want: alias -> string_protein_id
    # Filter for GLYMA_...
    
    mapping = {}
    with gzip.open(alias_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            string_id = parts[0]
            alias = parts[1]
            
            # Optimization: Only keep GLYMA_ IDs or Uniprot
            if 'GLYMA_' in alias or 'Glyma' in alias:
                # Normalize alias if needed? 
                # The file has "GLYMA_01G000700"
                mapping[alias] = string_id
                
    print(f"Loaded {len(mapping)} mappings.")
    return mapping

def run_mapping( diff_path, alias_path, output_path):
    print(f"Loading {diff_path}...")
    df = pd.read_csv(diff_path)
    
    alias_map = load_aliases(alias_path)
    
    mapped_ids = []
    
    for ids in df['Protein IDs']:
        # ids might be "Glyma.01G000700.1.p;Glyma..."
        # We try to map any of them
        found_id = None
        
        candidates = str(ids).split(';')
        for cand in candidates:
            # Clean candidate: Glyma.01G000700.1.p -> GLYMA_01G000700
            if 'Glyma.' in cand:
                # Extract 01G... part
                # Glyma.01G000700.1.p
                parts = cand.split('.')
                # usually parts[0]="Glyma", parts[1]="01G000700"
                if len(parts) >= 2:
                    core = parts[1] 
                    # specific fix for "Glyma.01G..." vs "GLYMA_01G..."
                    query = f"GLYMA_{core}"
                    if query in alias_map:
                        found_id = alias_map[query]
                        break
                        
            # Try direct match
            if cand in alias_map:
                found_id = alias_map[cand]
                break
                
        mapped_ids.append(found_id)
        
    df['STRING_ID'] = mapped_ids
    
    # Filter mapped
    mapped_df = df.dropna(subset=['STRING_ID'])
    print(f"Mapped: {len(mapped_df)} / {len(df)}")
    
    # Significant
    sig_df = mapped_df[(mapped_df['P_Value'] < 0.05) & (abs(mapped_df['Log2FC']) > 1)]
    print(f"Significant & Mapped: {len(sig_df)}")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    mapped_df.to_csv(output_path, index=False)
    print(f"Saved to {output_path}")

if __name__ == "__main__":
    run_mapping('data/processed/pxd006989_differential.csv', 
                'data/raw/3847.protein.aliases.v12.0.txt.gz',
                'data/processed/pxd006989_mapped.csv')
