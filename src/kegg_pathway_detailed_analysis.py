"""
Detailed KEGG pathway enrichment analysis with pathway name annotation.

This script enhances the existing KEGG pathway enrichment results by:
1. Fetching pathway names from KEGG API
2. Categorizing pathways by biological function
3. Generating detailed statistics for significant pathways
4. Creating publication-ready tables
"""

import pandas as pd
import urllib.request
import time
import os


def get_pathway_info(pathway_id: str) -> dict:
    """
    Fetch pathway information from KEGG API.

    Args:
        pathway_id: KEGG pathway ID (e.g., 'map01110')

    Returns:
        Dictionary with pathway info
    """
    info = {'id': pathway_id, 'name': None, 'class': None}

    try:
        url = f"http://rest.kegg.jp/get/{pathway_id}"
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8')

            for line in data.split('\n'):
                if line.startswith('NAME'):
                    info['name'] = line.replace('NAME', '').strip()
                elif line.startswith('CLASS'):
                    info['class'] = line.replace('CLASS', '').strip()

        time.sleep(0.3)  # Rate limiting
    except Exception as e:
        print(f"  Warning: Failed to fetch info for {pathway_id}: {e}")

    return info


def categorize_pathway(pathway_class: str, pathway_name: str) -> str:
    """
    Categorize pathway by biological function.

    Args:
        pathway_class: KEGG pathway class
        pathway_name: Pathway name

    Returns:
        Category name
    """
    if not pathway_class and not pathway_name:
        return "Unknown"

    text = f"{pathway_class or ''} {pathway_name or ''}".lower()

    if any(kw in text for kw in ['secondary metabolite', 'phenylpropanoid', 'flavonoid', 'isoflavonoid']):
        return "Secondary Metabolism"
    elif any(kw in text for kw in ['amino acid', 'aminoacyl', 'trna']):
        return "Amino Acid Metabolism"
    elif any(kw in text for kw in ['lipid', 'fatty acid', 'glycerophospholipid']):
        return "Lipid Metabolism"
    elif any(kw in text for kw in ['carbohydrate', 'glycolysis', 'tca']):
        return "Carbohydrate Metabolism"
    elif any(kw in text for kw in ['xenobiotics', 'degradation']):
        return "Xenobiotic Degradation"
    elif any(kw in text for kw in ['biosynthesis']):
        return "Biosynthesis"
    elif any(kw in text for kw in ['signal', 'transduction']):
        return "Signal Transduction"
    elif any(kw in text for kw in ['disease', 'infection', 'cancer']):
        return "Human Disease"
    else:
        return "Other Metabolism"


def analyze_kegg_pathways(
    enrichment_csv: str,
    metabolite_csv: str,
    output_csv: str,
    fetch_names: bool = True
) -> pd.DataFrame:
    """
    Detailed analysis of KEGG pathway enrichment results.

    Args:
        enrichment_csv: Path to pathway enrichment CSV
        metabolite_csv: Path to metabolite differential CSV
        output_csv: Path to output detailed analysis CSV
        fetch_names: Whether to fetch pathway names from KEGG

    Returns:
        DataFrame with detailed pathway analysis
    """
    print("=" * 70)
    print("KEGG Pathway Detailed Analysis")
    print("=" * 70)

    # Load enrichment results
    print(f"\nLoading enrichment results from {enrichment_csv}...")
    enrich_df = pd.read_csv(enrichment_csv)

    print(f"  Total pathways: {len(enrich_df)}")

    # Fetch pathway information
    if fetch_names:
        print(f"\nFetching pathway information from KEGG API...")
        pathway_info = []

        for idx, row in enrich_df.iterrows():
            pid = row['Pathway']
            print(f"  [{idx+1}/{len(enrich_df)}] {pid}")

            info = get_pathway_info(pid)
            pathway_info.append(info)

        # Merge with enrichment data
        info_df = pd.DataFrame(pathway_info)
        enrich_df = enrich_df.merge(info_df, left_on='Pathway', right_on='id', how='left')

        # Categorize pathways
        enrich_df['Category'] = enrich_df.apply(
            lambda row: categorize_pathway(row.get('class'), row.get('name')),
            axis=1
        )
    else:
        enrich_df['name'] = enrich_df['Pathway']
        enrich_df['class'] = None
        enrich_df['Category'] = "Unknown"

    # Load metabolite data to find which metabolites belong to each pathway
    print(f"\nLoading metabolite data from {metabolite_csv}...")
    met_df = pd.read_csv(metabolite_csv)

    # Get pathway members (requires additional KEGG API calls)
    print(f"\nIdentifying pathway members...")
    # This is a simplified version - full implementation would require
    # querying KEGG for each pathway's compounds

    # Calculate additional statistics
    enrich_df['Significant'] = enrich_df['P_Value'] < 0.05
    enrich_df['Highly_Significant'] = enrich_df['P_Value'] < 0.01
    enrich_df['-log10(P)'] = -np.log10(enrich_df['P_Value'] + 1e-300)  # Avoid log(0)

    # Sort by P-value
    enrich_df = enrich_df.sort_values('P_Value')

    # Save results
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    enrich_df.to_csv(output_csv, index=False)

    print(f"\n{'='*70}")
    print(f"Saved detailed analysis to {output_csv}")
    print(f"{'='*70}")

    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"  Total pathways: {len(enrich_df)}")
    print(f"  Significant (P < 0.05): {enrich_df['Significant'].sum()}")
    print(f"  Highly significant (P < 0.01): {enrich_df['Highly_Significant'].sum()}")

    # Category breakdown
    print(f"\nPathway Category Breakdown:")
    category_counts = enrich_df['Category'].value_counts()
    for cat, count in category_counts.items():
        print(f"  {cat}: {count}")

    # Top enriched pathways
    print(f"\nTop 10 Enriched Pathways:")
    top_cols = ['Pathway', 'name', 'Category', 'Sig_Count', 'P_Value', 'Enrichment_Score']
    available_cols = [col for col in top_cols if col in enrich_df.columns]
    print(enrich_df[available_cols].head(10).to_string(index=False))

    # Significant pathways detail
    sig_pathways = enrich_df[enrich_df['Significant']]
    if len(sig_pathways) > 0:
        print(f"\nSignificant Pathways (P < 0.05):")
        for idx, row in sig_pathways.iterrows():
            print(f"\n  {row['Pathway']}: {row.get('name', 'N/A')}")
            print(f"    P-value: {row['P_Value']:.4e}")
            print(f"    Significant metabolites: {row['Sig_Count']}")
            print(f"    Background metabolites: {row['Bg_Count']}")
            print(f"    Enrichment score: {row['Enrichment_Score']:.3f}")
            print(f"    Category: {row.get('Category', 'Unknown')}")

    return enrich_df


def generate_publication_table(detailed_csv: str, output_csv: str):
    """
    Generate publication-ready table from detailed analysis.

    Args:
        detailed_csv: Path to detailed analysis CSV
        output_csv: Path to output publication table CSV
    """
    print(f"\nGenerating publication-ready table...")

    df = pd.read_csv(detailed_csv)

    # Select and rename columns for publication
    pub_df = df[[
        'Pathway',
        'name',
        'Category',
        'Sig_Count',
        'Bg_Count',
        'P_Value',
        'Enrichment_Score'
    ]].copy()

    pub_df.columns = [
        'Pathway ID',
        'Pathway Name',
        'Category',
        'Significant Metabolites',
        'Total Metabolites',
        'P-value',
        'Enrichment Score'
    ]

    # Format P-values
    pub_df['P-value'] = pub_df['P-value'].apply(lambda x: f"{x:.4e}" if x < 0.001 else f"{x:.4f}")

    # Round enrichment scores
    pub_df['Enrichment Score'] = pub_df['Enrichment Score'].round(3)

    # Filter significant pathways
    df_sig = df[df['Significant']]
    if len(df_sig) > 0:
        pub_df_sig = pub_df.iloc[:len(df_sig)]
    else:
        pub_df_sig = pub_df.head(10)  # Top 10 if none significant

    # Save
    pub_df_sig.to_csv(output_csv, index=False)

    print(f"  Saved publication table to {output_csv}")
    print(f"  Included {len(pub_df_sig)} pathways")


if __name__ == "__main__":
    import numpy as np

    # File paths
    enrichment_csv = 'results/table1_metabolomics_real.csv'
    metabolite_csv = 'data/processed/mtbls531_differential.csv'
    output_detailed_csv = 'results/kegg_pathway_detailed.csv'
    output_publication_csv = 'results/kegg_pathway_publication_table.csv'

    # Check if enrichment file exists
    if not os.path.exists(enrichment_csv):
        print(f"Error: Enrichment file not found: {enrichment_csv}")
        print("\nPlease run KEGG pathway analysis first:")
        print("  python src/pathway_analysis.py")
        exit(1)

    # Run detailed analysis
    detailed_df = analyze_kegg_pathways(
        enrichment_csv,
        metabolite_csv,
        output_detailed_csv,
        fetch_names=True
    )

    # Generate publication table
    generate_publication_table(output_detailed_csv, output_publication_csv)

    print("\n" + "="*70)
    print("KEGG Pathway Detailed Analysis Complete!")
    print("="*70)
    print("\nOutput files:")
    print(f"  1. Detailed analysis: {output_detailed_csv}")
    print(f"  2. Publication table: {output_publication_csv}")
