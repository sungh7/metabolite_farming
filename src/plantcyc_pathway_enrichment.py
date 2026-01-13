"""
PlantCyc pathway enrichment analysis.

Performs Fisher's exact test to identify significantly enriched PlantCyc pathways
in ethylene-treated metabolomics data.
"""

import pandas as pd
import scipy.stats as stats
import os
from plantcyc_api import PlantCycClient
import time


def load_plantcyc_mappings(mapping_csv: str) -> pd.DataFrame:
    """
    Load PlantCyc metabolite-pathway mappings.

    Args:
        mapping_csv: Path to PlantCyc mapping CSV

    Returns:
        DataFrame with mappings
    """
    print(f"Loading PlantCyc mappings from {mapping_csv}...")
    df = pd.read_csv(mapping_csv)

    # Remove rows without pathway mappings
    df_valid = df.dropna(subset=['PlantCyc_Pathway_ID'])

    print(f"  Total rows: {len(df)}")
    print(f"  With pathway mappings: {len(df_valid)}")

    return df_valid


def merge_with_differential_data(
    plantcyc_df: pd.DataFrame,
    differential_csv: str
) -> pd.DataFrame:
    """
    Merge PlantCyc mappings with differential metabolomics data.

    Args:
        plantcyc_df: DataFrame with PlantCyc mappings
        differential_csv: Path to differential metabolomics CSV

    Returns:
        Merged DataFrame with statistical info
    """
    print(f"\nMerging with differential data from {differential_csv}...")
    diff_df = pd.read_csv(differential_csv)

    # Merge on metabolite name
    merged = plantcyc_df.merge(
        diff_df[['Name', 'Log2FC', 'P_Value', 'Control_Mean', 'Ethylene_Mean']],
        left_on='Metabolite_Name',
        right_on='Name',
        how='left'
    )

    print(f"  Merged rows: {len(merged)}")

    return merged


def get_pathway_names(pathway_ids: list, client: PlantCycClient) -> dict:
    """
    Retrieve pathway names from PlantCyc API.

    Args:
        pathway_ids: List of pathway IDs
        client: PlantCycClient instance

    Returns:
        Dictionary mapping pathway ID to name
    """
    print(f"\nRetrieving pathway names for {len(pathway_ids)} pathways...")
    pathway_names = {}

    for idx, pid in enumerate(pathway_ids):
        if idx % 10 == 0:
            print(f"  Progress: {idx}/{len(pathway_ids)}")

        info = client.get_pathway_info(pid)
        if info and 'name' in info:
            pathway_names[pid] = info['name']
        else:
            pathway_names[pid] = pid  # Use ID if name not available

        # Rate limiting
        time.sleep(1.1)

    return pathway_names


def run_enrichment_analysis(
    merged_df: pd.DataFrame,
    output_csv: str,
    p_threshold: float = 0.05,
    fetch_names: bool = False
) -> pd.DataFrame:
    """
    Run Fisher's exact test for pathway enrichment.

    Args:
        merged_df: DataFrame with pathway mappings and statistics
        output_csv: Path to output CSV
        p_threshold: P-value threshold for significance
        fetch_names: Whether to fetch pathway names from API

    Returns:
        DataFrame with enrichment results
    """
    print(f"\n{'='*60}")
    print("Running PlantCyc Pathway Enrichment Analysis")
    print(f"{'='*60}")

    # Define significant set
    sig_df = merged_df[merged_df['P_Value'] < p_threshold]
    total_sig = len(sig_df['Metabolite_Name'].unique())
    total_bg = len(merged_df['Metabolite_Name'].unique())

    print(f"\nSignificant metabolites (P < {p_threshold}): {total_sig}")
    print(f"Background metabolites: {total_bg}")

    # Count unique metabolites per pathway
    pathway_counts = {}

    for _, row in merged_df.iterrows():
        pid = row['PlantCyc_Pathway_ID']
        met_name = row['Metabolite_Name']
        is_sig = row['P_Value'] < p_threshold

        if pid not in pathway_counts:
            pathway_counts[pid] = {'sig': set(), 'bg': set()}

        pathway_counts[pid]['bg'].add(met_name)
        if is_sig:
            pathway_counts[pid]['sig'].add(met_name)

    # Convert sets to counts
    for pid in pathway_counts:
        pathway_counts[pid]['sig_count'] = len(pathway_counts[pid]['sig'])
        pathway_counts[pid]['bg_count'] = len(pathway_counts[pid]['bg'])

    # Fisher's exact test
    results = []

    print(f"\nCalculating enrichment for {len(pathway_counts)} pathways...")

    for pid, counts in pathway_counts.items():
        # Contingency table
        #        Sig   NotSig
        # InPath   a     b
        # NotPath  c     d

        a = counts['sig_count']
        b = counts['bg_count'] - a
        c = total_sig - a
        d = (total_bg - total_sig) - b

        # Skip if contingency table has negative values
        if a < 0 or b < 0 or c < 0 or d < 0:
            continue

        # Fisher's exact test (one-tailed, greater)
        try:
            odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        except ValueError:
            continue

        # Enrichment fold change
        if b > 0 and d > 0:
            fold_enrichment = (a / (a + b)) / (c / (c + d))
        else:
            fold_enrichment = float('inf') if a > 0 else 0

        results.append({
            'Pathway_ID': pid,
            'Sig_Count': a,
            'Bg_Count': counts['bg_count'],
            'P_Value': pval,
            'Fold_Enrichment': fold_enrichment,
            'Odds_Ratio': odds
        })

    # Convert to DataFrame and sort
    res_df = pd.DataFrame(results).sort_values('P_Value')

    # Fetch pathway names if requested
    if fetch_names:
        try:
            client = PlantCycClient(orgid='META')
            client.authenticate()

            unique_pathways = res_df['Pathway_ID'].unique().tolist()
            pathway_names = get_pathway_names(unique_pathways, client)

            res_df['Pathway_Name'] = res_df['Pathway_ID'].map(pathway_names)
        except Exception as e:
            print(f"\nWarning: Failed to fetch pathway names: {e}")
            res_df['Pathway_Name'] = res_df['Pathway_ID']
    else:
        res_df['Pathway_Name'] = res_df['Pathway_ID']

    # Save results
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    res_df.to_csv(output_csv, index=False)

    print(f"\n{'='*60}")
    print(f"Saved enrichment results to {output_csv}")
    print(f"{'='*60}")

    # Print top results
    print(f"\nTop 10 enriched pathways:")
    print(res_df[['Pathway_ID', 'Pathway_Name', 'Sig_Count', 'P_Value', 'Fold_Enrichment']].head(10).to_string(index=False))

    # Print summary
    sig_pathways = res_df[res_df['P_Value'] < 0.05]
    print(f"\nSummary:")
    print(f"  Total pathways tested: {len(res_df)}")
    print(f"  Significantly enriched (P < 0.05): {len(sig_pathways)}")
    print(f"  Highly significant (P < 0.01): {len(res_df[res_df['P_Value'] < 0.01])}")

    return res_df


def main():
    """Main function to run PlantCyc pathway enrichment analysis."""

    # File paths
    mapping_csv = 'data/processed/plantcyc_metabolite_pathways.csv'
    differential_csv = 'data/processed/mtbls531_differential.csv'
    output_csv = 'results/plantcyc_pathway_enrichment.csv'

    # Check if mapping file exists
    if not os.path.exists(mapping_csv):
        print(f"Error: Mapping file not found: {mapping_csv}")
        print("\nPlease run PlantCyc mapping first:")
        print("  python src/plantcyc_api.py map")
        print("\nSee docs/PLANTCYC_SETUP.md for instructions.")
        return

    # Load data
    plantcyc_df = load_plantcyc_mappings(mapping_csv)

    if len(plantcyc_df) == 0:
        print("\nError: No valid PlantCyc pathway mappings found.")
        print("This may indicate that:")
        print("  1. BioCyc authentication failed")
        print("  2. Metabolite names don't match PlantCyc database")
        print("  3. API requests are being blocked")
        print("\nPlease check:")
        print("  - BIOCYC_EMAIL and BIOCYC_PASSWORD environment variables")
        print("  - BioCyc account is active")
        print("  - Review mapping log: data/processed/plantcyc_mapping.log")
        return

    # Merge with differential data
    merged_df = merge_with_differential_data(plantcyc_df, differential_csv)

    # Run enrichment analysis
    # Note: Set fetch_names=True to retrieve pathway names (slower, requires API)
    enrichment_df = run_enrichment_analysis(
        merged_df,
        output_csv,
        p_threshold=0.05,
        fetch_names=False  # Set to True if you want to fetch pathway names
    )

    print("\n" + "="*60)
    print("PlantCyc Pathway Enrichment Analysis Complete!")
    print("="*60)


if __name__ == "__main__":
    main()
