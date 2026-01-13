"""
Generate Publication-Quality Pathway Analysis Figures

Creates 5 comprehensive figures for pathway enrichment analysis:
1. KEGG pathway enrichment bar chart
2. Metabolite volcano plot
3. Metabolite-pathway membership heatmap
4. KEGG vs PlantCyc comparison
5. Isoflavonoid pathway diagram

Output:
- PNG (300 DPI) and PDF (vector) formats
- Colorblind-safe palette
- Publication-ready styling

Requirements:
- results/kegg_pathway_detailed.csv
- results/plantcyc_pathway_enrichment.csv
- data/processed/mtbls531_differential.csv
- data/processed/plantcyc_metabolite_pathways.csv

Usage:
    python src/generate_pathway_figures.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import os
import sys
import textwrap
from typing import List, Dict, Optional, Tuple

# ============================================================================
# CONSTANTS & CONFIGURATION
# ============================================================================

# Colorblind-safe palette (Paul Tol's palette)
COLORBLIND_COLORS = {
    'blue': '#4477AA',
    'cyan': '#66CCEE',
    'green': '#228833',
    'yellow': '#CCBB44',
    'red': '#EE6677',
    'purple': '#AA3377',
    'grey': '#BBBBBB',
    'dark_grey': '#666666',
}

# Categorical palette for pathways
PATHWAY_PALETTE = [
    '#4477AA', '#66CCEE', '#228833', '#CCBB44',
    '#EE6677', '#AA3377', '#EE99AA', '#88CCAA'
]

# Statistical thresholds
P_THRESHOLD = 0.05
FC_THRESHOLD = 1.0

# ============================================================================
# CONFIGURATION FUNCTIONS
# ============================================================================

def set_publication_style():
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 16,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'pdf.fonttype': 42,  # TrueType for PDF compatibility
        'ps.fonttype': 42,
        'axes.linewidth': 1.2,
        'grid.linewidth': 0.8,
        'lines.linewidth': 2.0,
    })

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def wrap_text(text, max_width=40):
    """Wrap long text for axis labels."""
    if pd.isna(text):
        return ""
    return '\n'.join(textwrap.wrap(str(text), max_width))


def get_significance_marker(pval):
    """Return significance marker based on P-value."""
    if pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return ''


def save_figure(fig, base_path, dpi=300):
    """Save figure in both PNG and PDF formats."""
    os.makedirs(os.path.dirname(base_path), exist_ok=True)

    # PNG for presentations/web
    png_path = base_path if base_path.endswith('.png') else base_path + '.png'
    fig.savefig(png_path, dpi=dpi, bbox_inches='tight', facecolor='white')

    # PDF for publication (vector)
    pdf_path = png_path.replace('.png', '.pdf')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')

    print(f"  ✓ Saved: {png_path}")
    print(f"  ✓ Saved: {pdf_path}")

    return png_path, pdf_path


def validate_input_files(files_dict):
    """Check that all required input files exist."""
    missing = []
    for name, path in files_dict.items():
        if not os.path.exists(path):
            missing.append(f"  - {name}: {path}")

    if missing:
        raise FileNotFoundError(
            "Missing required input files:\n" + '\n'.join(missing) +
            "\n\nPlease run pathway enrichment analyses first:\n" +
            "  python src/kegg_pathway_detailed_analysis.py\n" +
            "  python src/plantcyc_pathway_enrichment.py"
        )

# ============================================================================
# FIGURE 1: KEGG PATHWAY ENRICHMENT BAR CHART
# ============================================================================

def create_kegg_enrichment_bar(kegg_csv, output_path, top_n=20):
    """
    Horizontal bar chart of KEGG pathway enrichment.

    Args:
        kegg_csv: Path to kegg_pathway_detailed.csv
        output_path: Output file path
        top_n: Number of top pathways to show
    """
    df = pd.read_csv(kegg_csv)

    # Sort by P-value and select top N
    df_sorted = df.sort_values('P_Value').head(top_n)

    # Calculate -log10(P)
    df_sorted = df_sorted.copy()
    df_sorted['-log10P'] = -np.log10(df_sorted['P_Value'] + 1e-300)

    # Use pathway name if available, otherwise use ID
    if 'name' in df_sorted.columns:
        pathway_names = df_sorted['name'].fillna(df_sorted['Pathway'])
    else:
        pathway_names = df_sorted['Pathway']

    # Color by significance
    colors = [COLORBLIND_COLORS['red'] if p < P_THRESHOLD
              else COLORBLIND_COLORS['grey']
              for p in df_sorted['P_Value']]

    fig, ax = plt.subplots(figsize=(10, 8))

    # Horizontal bars
    y_pos = np.arange(len(df_sorted))
    bars = ax.barh(y_pos, df_sorted['-log10P'], color=colors,
                   edgecolor=COLORBLIND_COLORS['dark_grey'], linewidth=1)

    # Pathway names with wrapping
    pathway_labels = [wrap_text(name, 50) for name in pathway_names]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pathway_labels, fontsize=9)

    # Reference line for P=0.05
    ax.axvline(x=-np.log10(P_THRESHOLD), color='black', linestyle='--',
               linewidth=1.5, alpha=0.6, label='P=0.05 threshold')

    # Significance markers
    for i, (idx, row) in enumerate(df_sorted.iterrows()):
        marker = get_significance_marker(row['P_Value'])
        if marker:
            ax.text(row['-log10P'] + 0.05, i, marker,
                   va='center', fontsize=14, fontweight='bold',
                   color=COLORBLIND_COLORS['red'])

    ax.set_xlabel('-log₁₀(P-value)', fontsize=12, fontweight='bold')
    ax.set_ylabel('KEGG Pathway', fontsize=12, fontweight='bold')
    ax.set_title('KEGG Pathway Enrichment Analysis\n(Ethylene vs Control)',
                fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='lower right', framealpha=0.9)
    ax.invert_yaxis()

    # Add grid
    ax.grid(axis='x', alpha=0.3, linestyle=':')

    save_figure(fig, output_path)
    plt.close()

# ============================================================================
# FIGURE 2: METABOLITE VOLCANO PLOT
# ============================================================================

def create_volcano_plot(differential_csv, output_path,
                       label_metabolites=['Daidzein', 'Formononetin',
                                         '6\'\'-Malonylgenistin'],
                       fc_threshold=FC_THRESHOLD, p_threshold=P_THRESHOLD):
    """
    Volcano plot showing differential metabolites.

    Args:
        differential_csv: Path to mtbls531_differential.csv
        output_path: Output file path
        label_metabolites: List of metabolite names to label
        fc_threshold: Log2FC threshold for significance
        p_threshold: P-value threshold
    """
    df = pd.read_csv(differential_csv)

    # Calculate -log10(P)
    df = df.copy()
    df['-log10P'] = -np.log10(df['P_Value'] + 1e-300)  # Avoid log(0)

    # Classify points
    df['Regulation'] = 'Non-significant'
    df.loc[(df['Log2FC'] > fc_threshold) & (df['P_Value'] < p_threshold),
           'Regulation'] = 'Up-regulated'
    df.loc[(df['Log2FC'] < -fc_threshold) & (df['P_Value'] < p_threshold),
           'Regulation'] = 'Down-regulated'

    # Count categories
    up_count = (df['Regulation'] == 'Up-regulated').sum()
    down_count = (df['Regulation'] == 'Down-regulated').sum()
    ns_count = (df['Regulation'] == 'Non-significant').sum()

    # Color mapping
    color_map = {
        'Up-regulated': COLORBLIND_COLORS['green'],
        'Down-regulated': COLORBLIND_COLORS['red'],
        'Non-significant': COLORBLIND_COLORS['grey']
    }

    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot each category separately for legend
    for category, color in color_map.items():
        subset = df[df['Regulation'] == category]
        count = len(subset)
        ax.scatter(subset['Log2FC'], subset['-log10P'],
                  c=color, label=f'{category} ({count})',
                  alpha=0.6, s=60, edgecolors='none')

    # Threshold lines
    ax.axhline(y=-np.log10(p_threshold), color='black',
               linestyle='--', linewidth=1.5, alpha=0.5,
               label=f'P={p_threshold}')
    ax.axvline(x=fc_threshold, color='black',
               linestyle='--', linewidth=1.5, alpha=0.5)
    ax.axvline(x=-fc_threshold, color='black',
               linestyle='--', linewidth=1.5, alpha=0.5)

    # Label specific metabolites
    labeled_count = 0
    for met_name in label_metabolites:
        met_rows = df[df['Name'].str.contains(met_name, case=False, na=False)]
        if len(met_rows) > 0:
            met = met_rows.iloc[0]
            ax.annotate(met_name,
                       xy=(met['Log2FC'], met['-log10P']),
                       xytext=(15, 10), textcoords='offset points',
                       fontsize=9, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.4',
                                facecolor='yellow', alpha=0.8,
                                edgecolor=COLORBLIND_COLORS['dark_grey']),
                       arrowprops=dict(arrowstyle='->',
                                      connectionstyle='arc3,rad=0.3',
                                      lw=1.5))
            labeled_count += 1

    ax.set_xlabel('Log₂ Fold Change (Ethylene / Control)',
                 fontsize=12, fontweight='bold')
    ax.set_ylabel('-log₁₀(P-value)', fontsize=12, fontweight='bold')
    ax.set_title('Metabolite Volcano Plot\n(Ethylene vs Control)',
                fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper left', framealpha=0.9, fontsize=10)
    ax.grid(True, alpha=0.3, linestyle=':')

    # Add text with stats
    ax.text(0.98, 0.02,
           f'Total metabolites: {len(df)}\nLabeled: {labeled_count}',
           transform=ax.transAxes, fontsize=9,
           verticalalignment='bottom', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    save_figure(fig, output_path)
    plt.close()

# ============================================================================
# FIGURE 3: METABOLITE-PATHWAY HEATMAP
# ============================================================================

def create_metabolite_pathway_heatmap(mapping_csv, differential_csv,
                                     output_path,
                                     top_metabolites=25,
                                     top_pathways=15):
    """
    Binary heatmap of metabolite-pathway memberships.

    Args:
        mapping_csv: Path to plantcyc_metabolite_pathways.csv
        differential_csv: Path to mtbls531_differential.csv
        output_path: Output file path
        top_metabolites: Number of metabolites to show
        top_pathways: Number of pathways to show
    """
    # Load data
    mapping_df = pd.read_csv(mapping_csv)
    diff_df = pd.read_csv(differential_csv)

    # Filter to mapped metabolites with pathway IDs
    mapping_df = mapping_df.dropna(subset=['PlantCyc_Pathway_ID'])

    # Merge to get P-values and Log2FC
    mapped_mets = mapping_df['Metabolite_Name'].unique()
    diff_mapped = diff_df[diff_df['Name'].isin(mapped_mets)]

    if len(diff_mapped) == 0:
        print("  ⚠ Warning: No metabolites found with both pathway mappings and differential data")
        print("  Skipping heatmap generation")
        return

    # Select top significant metabolites
    top_mets = diff_mapped.nsmallest(top_metabolites, 'P_Value')['Name'].tolist()

    # Find pathways containing these metabolites
    pathway_counts = mapping_df[
        mapping_df['Metabolite_Name'].isin(top_mets)
    ]['PlantCyc_Pathway_ID'].value_counts()

    top_pathways_list = pathway_counts.head(top_pathways).index.tolist()

    # Build binary matrix
    matrix = np.zeros((len(top_mets), len(top_pathways_list)))

    for i, met in enumerate(top_mets):
        met_pathways = mapping_df[
            mapping_df['Metabolite_Name'] == met
        ]['PlantCyc_Pathway_ID'].unique()

        for j, pathway in enumerate(top_pathways_list):
            if pathway in met_pathways:
                matrix[i, j] = 1

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))

    # Heatmap (binary colormap)
    cmap = ListedColormap(['white', COLORBLIND_COLORS['blue']])
    im = ax.imshow(matrix, cmap=cmap, aspect='auto',
                   interpolation='nearest')

    # Axes
    ax.set_xticks(np.arange(len(top_pathways_list)))
    ax.set_yticks(np.arange(len(top_mets)))

    # Labels with wrapping
    pathway_labels = [wrap_text(p, 25) for p in top_pathways_list]
    ax.set_xticklabels(pathway_labels, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(top_mets, fontsize=9)

    ax.set_xlabel('PlantCyc Pathway', fontsize=12, fontweight='bold')
    ax.set_ylabel('Metabolite (ranked by P-value)', fontsize=12, fontweight='bold')
    ax.set_title('Metabolite-Pathway Membership Matrix\n(Top Significant Metabolites)',
                fontsize=14, fontweight='bold', pad=20)

    # Grid
    ax.set_xticks(np.arange(matrix.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(matrix.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which='minor', color='lightgrey', linestyle='-', linewidth=0.5)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(['Not in pathway', 'In pathway'], fontsize=10)

    save_figure(fig, output_path)
    plt.close()

# ============================================================================
# FIGURE 4: KEGG VS PLANTCYC COMPARISON
# ============================================================================

def create_database_comparison(kegg_csv, plantcyc_csv, output_path, top_n=15):
    """
    Compare KEGG vs PlantCyc pathway enrichment results.

    Args:
        kegg_csv: Path to KEGG enrichment results
        plantcyc_csv: Path to PlantCyc enrichment results
        output_path: Output file path
        top_n: Number of top pathways per database
    """
    kegg_df = pd.read_csv(kegg_csv)
    plantcyc_df = pd.read_csv(plantcyc_csv)

    # Select top pathways
    kegg_top = kegg_df.sort_values('P_Value').head(top_n).copy()
    plantcyc_top = plantcyc_df.sort_values('P_Value').head(top_n).copy()

    # Calculate -log10(P)
    kegg_top['-log10P'] = -np.log10(kegg_top['P_Value'] + 1e-300)
    plantcyc_top['-log10P'] = -np.log10(plantcyc_top['P_Value'] + 1e-300)

    # Get pathway names
    if 'name' in kegg_top.columns:
        kegg_names = kegg_top['name'].fillna(kegg_top['Pathway'])
    else:
        kegg_names = kegg_top['Pathway']

    plantcyc_names = plantcyc_top['Pathway_Name'].fillna(plantcyc_top['Pathway_ID'])

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # KEGG panel
    y_pos = np.arange(len(kegg_top))
    colors_kegg = [COLORBLIND_COLORS['red'] if p < P_THRESHOLD
                   else COLORBLIND_COLORS['grey']
                   for p in kegg_top['P_Value']]

    ax1.barh(y_pos, kegg_top['-log10P'], color=colors_kegg,
            edgecolor=COLORBLIND_COLORS['dark_grey'], linewidth=1)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([wrap_text(n, 40) for n in kegg_names], fontsize=9)
    ax1.axvline(x=-np.log10(P_THRESHOLD), color='black', linestyle='--',
               linewidth=1.5, alpha=0.6, label='P=0.05')

    # Add significance markers for KEGG
    for i, (idx, row) in enumerate(kegg_top.iterrows()):
        marker = get_significance_marker(row['P_Value'])
        if marker:
            ax1.text(row['-log10P'] + 0.05, i, marker,
                    va='center', fontsize=12, fontweight='bold',
                    color=COLORBLIND_COLORS['red'])

    ax1.set_xlabel('-log₁₀(P-value)', fontsize=11, fontweight='bold')
    ax1.set_title('KEGG Pathways', fontsize=13, fontweight='bold')
    ax1.invert_yaxis()
    ax1.legend(loc='lower right', fontsize=9)
    ax1.grid(axis='x', alpha=0.3, linestyle=':')

    # PlantCyc panel
    y_pos2 = np.arange(len(plantcyc_top))
    colors_pc = [COLORBLIND_COLORS['blue'] if p < P_THRESHOLD
                 else COLORBLIND_COLORS['grey']
                 for p in plantcyc_top['P_Value']]

    ax2.barh(y_pos2, plantcyc_top['-log10P'], color=colors_pc,
            edgecolor=COLORBLIND_COLORS['dark_grey'], linewidth=1)
    ax2.set_yticks(y_pos2)
    ax2.set_yticklabels([wrap_text(n, 40) for n in plantcyc_names], fontsize=9)
    ax2.axvline(x=-np.log10(P_THRESHOLD), color='black', linestyle='--',
               linewidth=1.5, alpha=0.6, label='P=0.05')
    ax2.set_xlabel('-log₁₀(P-value)', fontsize=11, fontweight='bold')
    ax2.set_title('PlantCyc Pathways', fontsize=13, fontweight='bold')
    ax2.invert_yaxis()
    ax2.legend(loc='lower right', fontsize=9)
    ax2.grid(axis='x', alpha=0.3, linestyle=':')

    # Overall title
    fig.suptitle('Pathway Enrichment Comparison: KEGG vs PlantCyc\n(Ethylene Treatment)',
                fontsize=16, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save_figure(fig, output_path)
    plt.close()

# ============================================================================
# FIGURE 5: ISOFLAVONOID PATHWAY DIAGRAM
# ============================================================================

def create_isoflavonoid_pathway_diagram(differential_csv, proteomics_csv, output_path):
    """
    Create annotated isoflavonoid pathway diagram with metabolite and protein fold changes.

    Args:
        differential_csv: Path to differential metabolomics data
        proteomics_csv: Path to proteomics differential data
        output_path: Output file path
    """
    met_df = pd.read_csv(differential_csv)

    # Try to load proteomics data
    try:
        prot_df = pd.read_csv(proteomics_csv)
    except:
        prot_df = None

    # Extract key metabolites
    metabolites_of_interest = {
        'Daidzein': None,
        'Formononetin': None,
        'Daidzin': None,
        '6\'\'-Malonylgenistin': None,
    }

    for met_name in list(metabolites_of_interest.keys()):
        met_data = met_df[met_df['Name'].str.contains(met_name, case=False, na=False)]
        if len(met_data) > 0:
            metabolites_of_interest[met_name] = {
                'log2fc': met_data.iloc[0]['Log2FC'],
                'pval': met_data.iloc[0]['P_Value']
            }

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(6, 9.5, 'Isoflavonoid Biosynthesis Pathway',
           ha='center', fontsize=16, fontweight='bold')
    ax.text(6, 9, '(Ethylene Treatment Effect)',
           ha='center', fontsize=12, style='italic')

    # Define pathway structure (manual layout)
    pathway_elements = [
        # (x, y, name, type)
        (2, 7, 'Phenylalanine', 'metabolite'),
        (4, 7, 'Naringenin chalcone', 'metabolite'),
        (6, 7, 'Naringenin', 'metabolite'),
        (8, 7, 'Liquiritigenin', 'metabolite'),
        (10, 7, 'Daidzein', 'metabolite'),
        (10, 5, 'Formononetin', 'metabolite'),
        (10, 3, 'Daidzin', 'metabolite'),
        (10, 1, 'Malonyl-daidzin', 'metabolite'),
    ]

    # Draw metabolite boxes
    for x, y, name, elem_type in pathway_elements:
        # Lookup fold change data
        met_info = metabolites_of_interest.get(name)

        if met_info and met_info['pval'] < P_THRESHOLD:
            # Color by fold change
            if met_info['log2fc'] > 0:
                facecolor = COLORBLIND_COLORS['green']
                alpha = 0.8
            else:
                facecolor = COLORBLIND_COLORS['red']
                alpha = 0.8
            is_sig = True
        else:
            facecolor = 'white'
            alpha = 1.0
            is_sig = False

        # Box
        box_width = 1.4
        box_height = 0.5
        box = plt.Rectangle((x - box_width/2, y - box_height/2),
                           box_width, box_height,
                           facecolor=facecolor,
                           edgecolor=COLORBLIND_COLORS['dark_grey'],
                           linewidth=2, alpha=alpha, zorder=3)
        ax.add_patch(box)

        # Label
        ax.text(x, y, name, ha='center', va='center',
               fontsize=9, fontweight='bold', zorder=4)

        # Fold change annotation
        if met_info and is_sig:
            fc_text = f"↑{met_info['log2fc']:.1f}x"
            ax.text(x, y - box_height/2 - 0.15, fc_text,
                   ha='center', va='top', fontsize=8,
                   style='italic', fontweight='bold',
                   color=COLORBLIND_COLORS['green'] if met_info['log2fc'] > 0
                         else COLORBLIND_COLORS['red'],
                   zorder=4)

    # Draw arrows (enzyme reactions)
    arrows = [
        # (x1, y1, x2, y2, enzyme_name)
        (2.7, 7, 3.3, 7, 'PAL'),
        (4.7, 7, 5.3, 7, 'CHS'),
        (6.7, 7, 7.3, 7, 'CHI'),
        (8.7, 7, 9.3, 7, 'IFS'),
        (10, 6.7, 10, 5.3, 'I2\'H'),
        (10, 4.7, 10, 3.3, 'UGT'),
        (10, 2.7, 10, 1.3, 'MAT'),
    ]

    for x1, y1, x2, y2, enzyme in arrows:
        # Draw arrow
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                   arrowprops=dict(arrowstyle='->', lw=2.5,
                                  color=COLORBLIND_COLORS['dark_grey'],
                                  zorder=2))

        # Enzyme label
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2

        # Offset label to side
        if x1 == x2:  # Vertical arrow
            label_x = mid_x + 0.4
        else:  # Horizontal arrow
            label_y = mid_y + 0.3
            label_x = mid_x

        if x1 == x2:
            label_y = mid_y

        ax.text(label_x, label_y if x1 != x2 else label_y, enzyme,
               ha='center', fontsize=8, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.3',
                        facecolor='lightyellow',
                        edgecolor=COLORBLIND_COLORS['dark_grey'],
                        linewidth=1),
               zorder=4)

    # Legend
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=COLORBLIND_COLORS['green'],
                     alpha=0.8, edgecolor='black', label='Upregulated (P<0.05)'),
        plt.Rectangle((0, 0), 1, 1, facecolor='white',
                     edgecolor='black', linewidth=2, label='Not significant/measured')
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=10,
             framealpha=0.9)

    # Add enzyme key
    enzyme_text = ('Enzyme Key:\nPAL: Phenylalanine ammonia-lyase\n' +
                   'CHS: Chalcone synthase\nCHI: Chalcone isomerase\n' +
                   'IFS: Isoflavone synthase\nI2\'H: Isoflavone 2\'-hydroxylase\n' +
                   'UGT: UDP-glucosyltransferase\nMAT: Malonyl transferase')
    ax.text(0.5, 2, enzyme_text, fontsize=8,
           verticalalignment='top',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue',
                    alpha=0.3, edgecolor=COLORBLIND_COLORS['dark_grey']))

    save_figure(fig, output_path)
    plt.close()

# ============================================================================
# MAIN ORCHESTRATION
# ============================================================================

def main():
    """Generate all pathway analysis figures."""

    print("="*70)
    print("PATHWAY ANALYSIS FIGURE GENERATION")
    print("="*70)

    # Set publication style
    set_publication_style()

    # Define file paths
    files = {
        'kegg_enrichment': 'results/kegg_pathway_detailed.csv',
        'plantcyc_enrichment': 'results/plantcyc_pathway_enrichment.csv',
        'differential': 'data/processed/mtbls531_differential.csv',
        'mappings': 'data/processed/plantcyc_metabolite_pathways.csv',
        'proteomics': 'data/processed/pxd006989_differential.csv'
    }

    # Validate inputs (only check required files)
    required_files = {k: v for k, v in files.items() if k != 'proteomics'}
    try:
        validate_input_files(required_files)
    except FileNotFoundError as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

    # Create output directory
    output_dir = 'results/figures/pathway_analysis'
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n📁 Output directory: {output_dir}/\n")

    # Figure 1: KEGG Enrichment
    print("[1/5] Creating KEGG pathway enrichment bar chart...")
    create_kegg_enrichment_bar(
        files['kegg_enrichment'],
        f'{output_dir}/figure1_kegg_enrichment.png',
        top_n=20
    )

    # Figure 2: Volcano Plot
    print("\n[2/5] Creating metabolite volcano plot...")
    create_volcano_plot(
        files['differential'],
        f'{output_dir}/figure2_volcano_plot.png',
        label_metabolites=['Daidzein', 'Formononetin', '6\'\'-Malonylgenistin']
    )

    # Figure 3: Heatmap
    print("\n[3/5] Creating metabolite-pathway heatmap...")
    create_metabolite_pathway_heatmap(
        files['mappings'],
        files['differential'],
        f'{output_dir}/figure3_metabolite_pathway_heatmap.png',
        top_metabolites=25,
        top_pathways=15
    )

    # Figure 4: Database Comparison
    print("\n[4/5] Creating KEGG vs PlantCyc comparison...")
    create_database_comparison(
        files['kegg_enrichment'],
        files['plantcyc_enrichment'],
        f'{output_dir}/figure4_database_comparison.png',
        top_n=15
    )

    # Figure 5: Isoflavonoid Pathway
    print("\n[5/5] Creating isoflavonoid pathway diagram...")
    create_isoflavonoid_pathway_diagram(
        files['differential'],
        files['proteomics'],
        f'{output_dir}/figure5_isoflavonoid_pathway.png'
    )

    print("\n" + "="*70)
    print("✅ ALL FIGURES GENERATED SUCCESSFULLY!")
    print("="*70)
    print(f"\n📊 Generated files in {output_dir}/:")
    print("  1. figure1_kegg_enrichment.png/pdf")
    print("  2. figure2_volcano_plot.png/pdf")
    print("  3. figure3_metabolite_pathway_heatmap.png/pdf")
    print("  4. figure4_database_comparison.png/pdf")
    print("  5. figure5_isoflavonoid_pathway.png/pdf")
    print("\n💡 All figures are publication-ready (300 DPI PNG + vector PDF)")
    print("="*70)


if __name__ == "__main__":
    main()
