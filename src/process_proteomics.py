import pandas as pd
import numpy as np
import scipy.stats as stats
import os

def process_proteins(input_path, output_path):
    print(f"Loading {input_path}...")
    # Load specific columns to save memory
    # MaxQuant columns: 'Protein IDs', 'Gene names', 'LFQ intensity Control...', 'LFQ intensity Flooding...'
    # We need to inspect columns first usually.
    # But let's assume standard LFQ columns exist.
    
    try:
        df = pd.read_csv(input_path, sep='\t', nrows=5)
        print("Columns:", df.columns.tolist())
        
        # Identify LFQ columns
        lfq_cols = [c for c in df.columns if 'LFQ intensity' in c]
        print(f"LFQ Columns: {lfq_cols}")
        
    except Exception as e:
        print(f"Error reading header: {e}")
        return

    # Load full (or relevant)
    cols = ['Protein IDs', 'Gene names', 'Fasta headers'] + lfq_cols
    df = pd.read_csv(input_path, sep='\t', usecols=lambda c: c in cols)
    
    # 0 -> NaN
    df[lfq_cols] = df[lfq_cols].replace(0, np.nan)
    
    # Impute with global min / 2 (Simple approach)
    min_val = df[lfq_cols].min().min()
    impute_val = min_val / 2 if pd.notna(min_val) else 1e4
    df[lfq_cols] = df[lfq_cols].fillna(impute_val)
    
    # Identify groups: Control vs Flood at different Days?
    # Expected: LFQ intensity C2, C4, C6, F2, F4, F6?
    # We will print columns to debug first, but let's write generic logic.
    # Group by prefix?
    
    # Calculate Max Log2FC across all contrasts
    # Assume we have sets of C and F.
    # We'll calculate Log2(Mean_F / Mean_C) for each day if possible, 
    # OR just Log2(Mean_F_all / Mean_C_all) if columns are messy.
    
    # Let's try to detect pairs from column names.
    # "LFQ intensity C_2d", "LFQ intensity F_2d"
    
    results = []
    
    for idx, row in df.iterrows():
        # Heuristic: 
        # C columns: 'control'
        # F columns: 'ethylene' (and NOT 'ABA')
        
        c_vals = [row[c] for c in lfq_cols if 'control' in c.lower()]
        # Exclude 'ABA' to get pure ethylene
        f_vals = [row[c] for c in lfq_cols if 'ethylene' in c.lower() and 'aba' not in c.lower()]
        
        if not c_vals or not f_vals:
            continue
            
        mean_c = np.mean(c_vals)
        mean_f = np.mean(f_vals)
        
        log2fc = np.log2(mean_f / mean_c)
        
        # T-test
        try:
            _, pval = stats.ttest_ind(f_vals, c_vals, equal_var=False)
        except:
            pval = 1.0
            
        results.append({
            'Protein IDs': row['Protein IDs'],
            'Gene names': row['Gene names'],
            'Protein names': row.get('Protein names', ''),
            'Fasta headers': row.get('Fasta headers', ''),
            'Log2FC': log2fc,
            'P_Value': pval if pd.notna(pval) else 1.0,
            'Mean_Control': mean_c,
            'Mean_Ethylene': mean_f
        })
        
    res_df = pd.DataFrame(results)
    sig_df = res_df[(res_df['P_Value'] < 0.05) & (abs(res_df['Log2FC']) > 1)]
    
    print(f"Total Proteins: {len(df)}")
    print(f"Significant (P<0.05, FC>2): {len(sig_df)}")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # Ensure columns exist in output even if input didn't have them (get default '')
    res_df.to_csv(output_path, index=False)
    print(f"Saved to {output_path}")

if __name__ == "__main__":
    process_proteins('data/experimental/proteomics/proteinGroups.txt', 'data/processed/pxd006989_differential.csv')
