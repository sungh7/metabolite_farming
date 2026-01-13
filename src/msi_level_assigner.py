"""
Assign MSI Levels to Metabolomics Features.

MSI Levels (Metabolomics Standards Initiative):
- Level 1: Confirmed (RT + MS/MS match with authentic standard) - NOT AVAILABLE
- Level 2: Putatively Annotated (Database match, e.g., ChEBI + KEGG)
- Level 3: Putatively Characterized (Class-level, e.g., ChEBI only)
- Level 4: Unknown (No annotation)

For this dataset (MTBLS531):
- Level 2: ChEBI ID + Name + KEGG ID (strong database evidence)
- Level 3: ChEBI ID + Name only (moderate database evidence)
- Level 4: No ChEBI ID

This is a conservative assignment as we lack MS/MS spectral library matching info.
"""

import pandas as pd
import os

def assign_msi_level(row):
    """Assign MSI Level based on available annotation evidence."""
    has_chebi = pd.notna(row.get('ChEBI')) and str(row.get('ChEBI')).strip() != ''
    has_name = pd.notna(row.get('Name')) and str(row.get('Name')).strip() != ''
    has_kegg = pd.notna(row.get('KEGG')) and str(row.get('KEGG')).strip() != ''
    
    if has_chebi and has_name and has_kegg:
        return 2  # Strong database evidence (multiple DBs agree)
    elif has_chebi and has_name:
        return 3  # Moderate database evidence (single DB)
    elif has_chebi or has_name:
        return 3  # Weak but still characterized
    else:
        return 4  # Unknown

def main():
    input_path = 'results/metabolomics/top_features_up.tsv'
    output_path = 'results/metabolomics/top_features_up_with_msi.tsv'
    
    df = pd.read_csv(input_path, sep='\t')
    
    # Assign MSI Levels
    df['MSI_Level'] = df.apply(assign_msi_level, axis=1)
    
    # Reorder columns
    cols = ['ChEBI', 'Name', 'MSI_Level', 'Log2FC', 'P_Value', 'FDR', 'KEGG']
    df = df[cols]
    
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved {len(df)} features with MSI Levels to {output_path}")
    print(df[['Name', 'MSI_Level', 'Log2FC']].to_string())

if __name__ == "__main__":
    main()
