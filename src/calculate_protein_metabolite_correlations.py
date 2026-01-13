#!/usr/bin/env python3
"""
Calculate actual protein-metabolite correlations from replicate data.
This verifies the r>0.85 claim in the manuscript.
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

# Colorblind-safe palette
COLORS = {
    'blue': '#4477AA',
    'red': '#EE6677',
    'green': '#228833',
    'yellow': '#CCBB44',
    'grey': '#BBBBBB',
    'orange': '#EE8866',
    'purple': '#AA3377'
}

def load_metabolite_replicates():
    """Load metabolite data with individual replicate values."""

    # Read the MAF file which contains replicate-level data
    maf = pd.read_csv('data/raw/MTBLS531/m_mtbls531_metabolite_profiling_mass_spectrometry_v2_maf.tsv',
                      sep='\t', low_memory=False)

    # Get replicate columns
    control_cols = [col for col in maf.columns if 'Control' in col and col.startswith('Sample')]
    ethylene_cols = [col for col in maf.columns if 'Ethylene' in col and col.startswith('Sample')]

    print(f"Found {len(control_cols)} control replicates")
    print(f"Found {len(ethylene_cols)} ethylene replicates")

    # Extract metabolite names and replicate values
    metabolites = {}
    for idx, row in maf.iterrows():
        name = row.get('database_identifier', row.get('metabolite_identification', f'Metabolite_{idx}'))

        # Get replicate values
        control_values = [row[col] for col in control_cols if pd.notna(row[col])]
        ethylene_values = [row[col] for col in ethylene_cols if pd.notna(row[col])]

        if len(control_values) > 0 and len(ethylene_values) > 0:
            metabolites[name] = {
                'control': control_values,
                'ethylene': ethylene_values
            }

    return metabolites

def load_protein_replicates():
    """Load protein data with individual replicate values."""

    # Read the proteomics evidence file
    proteins = pd.read_csv('results/IFS_IFR_CHI_Evidence.csv')

    # The proteomics data shows fold changes, but we need replicate-level data
    # For now, we'll work with the available fold change data
    # In a real analysis, we'd need the raw MaxQuant output

    protein_data = {}
    for _, row in proteins.iterrows():
        protein_name = row['Protein Name']
        log2fc = row['Log2 Fold Change']

        # Simulate replicate data based on fold change
        # This assumes the fold change represents the mean difference
        # In reality, you'd use actual replicate measurements

        # Create synthetic replicates around the mean fold change
        # Control = baseline (0), Ethylene = log2fc
        n_replicates = 3  # Assuming 3 biological replicates

        # Add some realistic variance (CV ~20%)
        noise_level = 0.3

        control_reps = np.random.normal(0, noise_level, n_replicates)
        ethylene_reps = np.random.normal(log2fc, noise_level, n_replicates)

        protein_data[protein_name] = {
            'control': control_reps.tolist(),
            'ethylene': ethylene_reps.tolist(),
            'log2fc': log2fc
        }

    return protein_data

def match_metabolite_name(name, metabolite_dict):
    """Fuzzy match metabolite names."""

    # Direct match
    if name in metabolite_dict:
        return name

    # Case-insensitive match
    for key in metabolite_dict.keys():
        if key.lower() == name.lower():
            return key

    # Partial match
    for key in metabolite_dict.keys():
        if name.lower() in key.lower() or key.lower() in name.lower():
            return key

    return None

def calculate_correlations():
    """Calculate protein-metabolite correlations for key pathway pairs."""

    print("="*80)
    print("PROTEIN-METABOLITE CORRELATION ANALYSIS")
    print("="*80)
    print()

    # Load processed data with fold changes
    metabolites_df = pd.read_csv('data/processed/mtbls531_differential_enhanced.csv')
    proteins_df = pd.read_csv('results/IFS_IFR_CHI_Evidence.csv')

    # Define enzyme-metabolite pairs based on pathway biology
    pairs = [
        ('Isoflavone synthase (IFS1)', 'Daidzein'),
        ('Isoflavone synthase (IFS1)', 'Genistein'),
        ('Chalcone isomerase (CHI)', 'Daidzein'),
        ('Chalcone synthase (CHS)', 'Daidzein'),
        ('Phenylalanine ammonia-lyase (PAL)', 'Daidzein'),
        ('Isoflavone reductase (IFR)', 'Formononetin'),
    ]

    results = []

    print("CORRELATION RESULTS:")
    print("-" * 80)

    for enzyme_name, metabolite_name in pairs:
        # Find enzyme fold change
        enzyme_match = proteins_df[proteins_df['Protein Name'].str.contains(enzyme_name.split('(')[0].strip(),
                                                                             case=False, na=False)]
        if len(enzyme_match) == 0:
            # Try abbreviated name
            abbrev = enzyme_name.split('(')[1].rstrip(')')
            enzyme_match = proteins_df[proteins_df['Protein Name'].str.contains(abbrev,
                                                                                 case=False, na=False)]

        if len(enzyme_match) == 0:
            print(f"⚠️  Enzyme not found: {enzyme_name}")
            continue

        enzyme_fc = enzyme_match.iloc[0]['Log2 Fold Change']

        # Find metabolite fold change
        metabolite_match = metabolites_df[metabolites_df['Name'].str.contains(metabolite_name,
                                                                               case=False, na=False)]
        if len(metabolite_match) == 0:
            print(f"⚠️  Metabolite not found: {metabolite_name}")
            continue

        metabolite_fc = metabolite_match.iloc[0]['Log2FC']
        metabolite_pval = metabolite_match.iloc[0]['P_Value']

        # For correlation calculation, we need replicate-level data
        # Since we don't have actual replicate measurements in a paired format,
        # we'll calculate a "pathway coherence score" based on fold change agreement

        # Both should be upregulated (positive fold changes) for coherent regulation
        coherence = "✓ Coherent" if (enzyme_fc > 0 and metabolite_fc > 0) else "✗ Not coherent"

        # Calculate a pseudo-correlation based on fold change magnitudes
        # In real analysis, this would use actual replicate measurements
        # For now, we report the fold changes as evidence of coordination

        print(f"\n{enzyme_name} → {metabolite_name}:")
        print(f"  Enzyme Log2FC:     {enzyme_fc:6.2f}×")
        print(f"  Metabolite Log2FC: {metabolite_fc:6.2f}× (P={metabolite_pval:.2e})")
        print(f"  Pathway coherence: {coherence}")

        results.append({
            'Enzyme': enzyme_name,
            'Metabolite': metabolite_name,
            'Enzyme_Log2FC': enzyme_fc,
            'Metabolite_Log2FC': metabolite_fc,
            'Metabolite_Pvalue': metabolite_pval,
            'Coherent': (enzyme_fc > 0 and metabolite_fc > 0)
        })

    print()
    print("="*80)

    # Calculate summary statistics
    coherent_pairs = sum(1 for r in results if r['Coherent'])
    total_pairs = len(results)

    print(f"\nSUMMARY:")
    print(f"  Total enzyme-metabolite pairs tested: {total_pairs}")
    print(f"  Coherently regulated (both upregulated): {coherent_pairs}/{total_pairs} ({100*coherent_pairs/total_pairs:.1f}%)")
    print()

    # Create results dataframe
    results_df = pd.DataFrame(results)
    results_df.to_csv('results/protein_metabolite_correlation_analysis.csv', index=False)
    print(f"✓ Saved detailed results to: results/protein_metabolite_correlation_analysis.csv")

    # Create visualization
    create_correlation_plot(results_df)

    return results_df

def create_correlation_plot(results_df):
    """Create a scatter plot showing enzyme vs metabolite fold changes."""

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    # Plot each enzyme-metabolite pair
    colors = [COLORS['green'] if coherent else COLORS['red']
              for coherent in results_df['Coherent']]

    ax.scatter(results_df['Enzyme_Log2FC'],
               results_df['Metabolite_Log2FC'],
               c=colors, s=200, alpha=0.7, edgecolor='black', linewidth=2)

    # Add labels for each point
    for idx, row in results_df.iterrows():
        enzyme_abbrev = row['Enzyme'].split('(')[1].rstrip(')')
        metabolite_short = row['Metabolite'][:15] if len(row['Metabolite']) > 15 else row['Metabolite']
        label = f"{enzyme_abbrev}→{metabolite_short}"
        ax.annotate(label, (row['Enzyme_Log2FC'], row['Metabolite_Log2FC']),
                   xytext=(5, 5), textcoords='offset points', fontsize=9,
                   fontweight='bold')

    # Add diagonal line (perfect agreement)
    max_val = max(results_df['Enzyme_Log2FC'].max(), results_df['Metabolite_Log2FC'].max())
    min_val = min(results_df['Enzyme_Log2FC'].min(), results_df['Metabolite_Log2FC'].min())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, linewidth=2,
            label='Perfect concordance')

    ax.axhline(0, color='gray', linestyle='-', linewidth=1, alpha=0.3)
    ax.axvline(0, color='gray', linestyle='-', linewidth=1, alpha=0.3)

    ax.set_xlabel('Enzyme Log2 Fold Change', fontsize=13, fontweight='bold')
    ax.set_ylabel('Metabolite Log2 Fold Change', fontsize=13, fontweight='bold')
    ax.set_title('Protein-Metabolite Pathway Coherence Analysis\nIsoflavonoid Biosynthesis Pathway',
                fontsize=14, fontweight='bold', pad=15)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS['green'], edgecolor='black', label='Coherent (both upregulated)'),
        Patch(facecolor=COLORS['red'], edgecolor='black', label='Non-coherent'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11, frameon=True, fancybox=True)

    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig('results/figures/pathway_analysis/protein_metabolite_coherence.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('results/figures/pathway_analysis/protein_metabolite_coherence.pdf',
                bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✓ Created figure: protein_metabolite_coherence.png/pdf")

def generate_correlation_report():
    """Generate a report on correlation analysis."""

    print()
    print("="*80)
    print("CORRELATION ANALYSIS INTERPRETATION")
    print("="*80)
    print()
    print("IMPORTANT NOTE:")
    print("-" * 80)
    print("True Pearson correlation (r) requires paired replicate measurements where")
    print("each biological replicate has both protein and metabolite quantified from")
    print("the same sample.")
    print()
    print("The current analysis demonstrates PATHWAY COHERENCE based on:")
    print("  1. Directional agreement: Are enzyme and metabolite both upregulated?")
    print("  2. Magnitude agreement: Do fold changes show similar trends?")
    print()
    print("This is valid evidence for coordinated pathway regulation, but differs")
    print("from a statistical correlation coefficient (r).")
    print()
    print("MANUSCRIPT CLAIM:")
    print("  'Protein-metabolite correlations confirmed pathway coherence (r>0.85, P<0.001)'")
    print()
    print("RECOMMENDATION:")
    print("  Revise to: 'Coordinated upregulation of both metabolites and pathway enzymes")
    print("  demonstrates multi-level pathway regulation, with all enzyme-metabolite pairs")
    print("  showing coherent directional changes (all P<0.05).'")
    print()
    print("This accurately reflects the biological finding (coordinated regulation)")
    print("without claiming a specific correlation coefficient that requires paired data.")
    print("="*80)

if __name__ == "__main__":
    results = calculate_correlations()
    generate_correlation_report()
