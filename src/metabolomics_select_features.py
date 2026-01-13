import pandas as pd
import numpy as np
import os
from statsmodels.stats.multitest import multipletests

def select_features(input_path, output_dir):
    print(f"Loading {input_path}...")
    df = pd.read_csv(input_path)
    
    # 1. QC Rules (Already processed mean/log2fc in previous step, checking detection logic)
    # Rules: Log2FC >= 1.0 (or <= -1.0), FDR <= 0.05
    
    # Calculate FDR if not present
    if 'FDR' not in df.columns:
        print("Calculating FDR...")
        # Drop NaNs in P_Value for calculation
        mask = df['P_Value'].notna()
        pvals = df.loc[mask, 'P_Value']
        reject, qvals, _, _ = multipletests(pvals, method='fdr_bh')
        df.loc[mask, 'FDR'] = qvals
    
    # 2. Filter
    # Up: FC >= 1.0, FDR <= 0.05
    up_mask = (df['Log2FC'] >= 1.0) & (df['FDR'] <= 0.05)
    up_features = df[up_mask].copy()
    up_features['Direction'] = 'Up'
    
    # Down: FC <= -1.0, FDR <= 0.05
    down_mask = (df['Log2FC'] <= -1.0) & (df['FDR'] <= 0.05)
    down_features = df[down_mask].copy()
    down_features['Direction'] = 'Down'
    
    print(f"Found {len(up_features)} Upregulated features.")
    print(f"Found {len(down_features)} Downregulated features.")
    
    # 3. Save
    os.makedirs(output_dir, exist_ok=True)
    
    up_path = os.path.join(output_dir, 'top_features_up.tsv')
    down_path = os.path.join(output_dir, 'top_features_down.tsv')
    
    # Select columns
    cols = ['ChEBI', 'Name', 'Log2FC', 'P_Value', 'FDR', 'KEGG']
    
    up_features[cols].to_csv(up_path, sep='\t', index=False)
    down_features[cols].to_csv(down_path, sep='\t', index=False)
    
    print(f"Saved to {up_path} and {down_path}")
    
    # Also save annotation level summary (Mock based on name)
    # If name starts with 'Compound' or similar -> Level 4. Else Level 2.
    annot_df = pd.concat([up_features, down_features])
    def get_level(name):
        if 'Unknown' in str(name): return 4
        return 2 # Assumed ID match
    
    annot_df['Level'] = annot_df['Name'].apply(get_level)
    annot_path = os.path.join(output_dir, 'annot_levels.tsv')
    annot_df[['ChEBI', 'Name', 'Level']].to_csv(annot_path, sep='\t', index=False)
    print("Saved annotation levels.")

if __name__ == "__main__":
    select_features('data/processed/mtbls531_differential.csv', 'results/metabolomics')
