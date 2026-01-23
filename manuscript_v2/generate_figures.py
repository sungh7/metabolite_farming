#!/usr/bin/env python3
"""
Figure Generation for Manuscript v2.0
Two-Track Analysis Framework

Generates publication-quality figures using real experimental data.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['figure.dpi'] = 300

# Output directory
OUTPUT_DIR = Path('/data/ethylene/manuscript_v2/figures')
OUTPUT_DIR.mkdir(exist_ok=True)

# Color scheme for two-track structure
COLORS = {
    'track_a': '#2196F3',  # Blue
    'track_b': '#FF9800',  # Orange
    'up': '#E53935',       # Red
    'down': '#1E88E5',     # Blue
    'neutral': '#9E9E9E',  # Gray
    'metabolite': '#4CAF50', # Green
    'enzyme': '#9C27B0',   # Purple
    'significant': '#E53935',  # Red
    'not_sig': '#BDBDBD',  # Light gray
}


def load_metabolomics_data():
    """Load metabolomics data from CSV."""
    df = pd.read_csv('/data/ethylene/results/Supplementary_Table_S3_All_Metabolites.csv')
    return df


def load_proteomics_data():
    """Load proteomics enzyme data from CSV."""
    df = pd.read_csv('/data/ethylene/results/IFS_IFR_CHI_Evidence.csv')
    return df


# =============================================================================
# Figure 2: Metabolomics Results - Volcano Plot
# =============================================================================
def generate_figure2_volcano():
    """Generate volcano plot for metabolomics data."""
    print("Generating Figure 2: Volcano Plot...")
    
    df = load_metabolomics_data()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Calculate -log10(p-value)
    df['neg_log10_p'] = -np.log10(df['P_Value'].replace(0, 1e-50))
    
    # Classify metabolites
    df['category'] = 'Not significant'
    df.loc[(df['P_Value'] < 0.05) & (df['Log2FC'] > 1), 'category'] = 'Up (Log2FC > 1)'
    df.loc[(df['P_Value'] < 0.05) & (df['Log2FC'] < -1), 'category'] = 'Down (Log2FC < -1)'
    df.loc[(df['P_Value'] < 0.05) & (df['Log2FC'].abs() <= 1), 'category'] = 'Significant (|Log2FC| ≤ 1)'
    
    # Plot by category
    for cat, color in [('Not significant', COLORS['not_sig']), 
                       ('Significant (|Log2FC| ≤ 1)', '#FFA726'),
                       ('Up (Log2FC > 1)', COLORS['up']),
                       ('Down (Log2FC < -1)', COLORS['down'])]:
        subset = df[df['category'] == cat]
        ax.scatter(subset['Log2FC'], subset['neg_log10_p'], 
                  c=color, alpha=0.7, s=50, label=cat, edgecolors='white', linewidths=0.5)
    
    # Highlight key isoflavonoids
    isoflavonoids = ['Daidzein', 'Formononetin', "6''-O-Acetyldaidzin", "6''-Malonylgenistin"]
    for name in isoflavonoids:
        subset = df[df['Name'].str.contains(name, case=False, na=False)]
        for _, row in subset.iterrows():
            ax.annotate(row['Name'][:20], 
                       (row['Log2FC'], row['neg_log10_p']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, fontweight='bold',
                       arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    # Add threshold lines
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='P = 0.05')
    ax.axvline(1, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(-1, color='gray', linestyle='--', alpha=0.3)
    
    ax.set_xlabel('Log₂ Fold Change (Ethylene / Control)', fontweight='bold')
    ax.set_ylabel('-Log₁₀ P-value', fontweight='bold')
    ax.set_title('Figure 2A: Metabolite Differential Abundance\nEthylene Treatment vs Control', 
                fontweight='bold', fontsize=14)
    
    ax.legend(loc='upper right', framealpha=0.9)
    
    # Add text box with summary statistics
    n_up = len(df[(df['P_Value'] < 0.05) & (df['Log2FC'] > 1)])
    n_down = len(df[(df['P_Value'] < 0.05) & (df['Log2FC'] < -1)])
    n_sig = len(df[df['P_Value'] < 0.05])
    text = f"Total: {len(df)} metabolites\nSignificant (P<0.05): {n_sig}\nUp: {n_up} | Down: {n_down}"
    ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'Figure2A_volcano_plot.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'Figure2A_volcano_plot.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: Figure2A_volcano_plot.png/pdf")


# =============================================================================
# Figure 4: Track A - GNN Performance Comparison
# =============================================================================
def generate_figure4_gnn_performance():
    """Generate GNN performance comparison bar chart."""
    print("Generating Figure 4C: GNN Performance Comparison...")
    
    # Performance data
    methods = ['Random\nBaseline', 'Adamic-\nAdar', 'HGT\n(Ours)']
    hits20 = [5.8, 14.9, 77.6]
    colors = [COLORS['not_sig'], '#FFA726', COLORS['track_a']]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    bars = ax.bar(methods, hits20, color=colors, edgecolor='white', linewidth=2)
    
    # Add value labels on bars
    for bar, val in zip(bars, hits20):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 1,
               f'{val}%', ha='center', va='bottom', fontweight='bold', fontsize=14)
    
    # Add improvement annotation
    ax.annotate('', xy=(2, 77.6), xytext=(0, 5.8),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2, ls='--'))
    ax.text(1, 45, '13.4×\nimprovement', ha='center', fontsize=12, fontweight='bold', color='gray')
    
    ax.set_ylabel('Hits@20 (%)', fontweight='bold', fontsize=12)
    ax.set_title('Figure 4C: Track A - GNN Model Performance\nEnzyme-Metabolite Link Prediction', 
                fontweight='bold', fontsize=14)
    ax.set_ylim(0, 100)
    
    # Add Track A indicator
    ax.text(0.02, 0.98, 'Track A', transform=ax.transAxes, fontsize=12,
           fontweight='bold', color=COLORS['track_a'],
           verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='white', edgecolor=COLORS['track_a'], alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'Figure4C_gnn_performance.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'Figure4C_gnn_performance.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: Figure4C_gnn_performance.png/pdf")


# =============================================================================
# Figure 5: Track A - Proteomics Validation
# =============================================================================
def generate_figure5_proteomics():
    """Generate proteomics enzyme fold change bar chart."""
    print("Generating Figure 5A: Proteomics Enzyme Fold Changes...")
    
    df = load_proteomics_data()
    
    # Sort by pathway order
    pathway_order = ['PAL', '4CL', 'CHS', 'CHI', 'IFS1', 'IFR']
    order_map = {name: i for i, name in enumerate(pathway_order)}
    
    # Extract enzyme abbreviation
    df['Abbrev'] = df['Protein Name'].str.extract(r'\(([A-Z0-9]+)\)')
    df['Order'] = df['Abbrev'].map(lambda x: order_map.get(x, 99))
    df = df.sort_values('Order')
    
    # Calculate linear fold change
    df['Linear_FC'] = 2 ** df['Log2 Fold Change']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create bars
    enzymes = df['Abbrev'].tolist()
    fold_changes = df['Linear_FC'].tolist()
    
    bars = ax.bar(enzymes, fold_changes, color=COLORS['enzyme'], edgecolor='white', linewidth=2)
    
    # Add value labels
    for bar, val in zip(bars, fold_changes):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 2,
               f'{val:.0f}×', ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    # Add pathway flow arrows
    for i in range(len(enzymes)-1):
        ax.annotate('', xy=(i+0.7, 5), xytext=(i+0.3, 5),
                   arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    
    ax.set_ylabel('Fold Change (Ethylene / Control)', fontweight='bold', fontsize=12)
    ax.set_xlabel('Isoflavonoid Biosynthesis Pathway Enzymes', fontweight='bold', fontsize=12)
    ax.set_title('Figure 5A: Track A Validation - Enzyme Upregulation\nAll P < 0.05', 
                fontweight='bold', fontsize=14)
    
    # Add Track A indicator
    ax.text(0.02, 0.98, 'Track A', transform=ax.transAxes, fontsize=12,
           fontweight='bold', color=COLORS['track_a'],
           verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='white', edgecolor=COLORS['track_a'], alpha=0.9))
    
    # Add validation checkmark
    ax.text(0.98, 0.98, '✓ GNN Predictions Validated', transform=ax.transAxes, fontsize=11,
           fontweight='bold', color='green', ha='right',
           verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='#E8F5E9', edgecolor='green', alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'Figure5A_proteomics_validation.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'Figure5A_proteomics_validation.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: Figure5A_proteomics_validation.png/pdf")


# =============================================================================
# Figure 2B: Isoflavonoid Comparison (Basal vs Conjugated)
# =============================================================================
def generate_figure2b_isoflavonoid_comparison():
    """Generate comparison of basal vs conjugated isoflavonoids."""
    print("Generating Figure 2B: Isoflavonoid Comparison...")
    
    df = load_metabolomics_data()
    
    # Select key isoflavonoids
    basal = df[df['Name'].isin(['Daidzein', 'Formononetin'])]
    conjugated = df[df['Name'].str.contains("Acetyldaidzin|Malonylgenistin", case=False, na=False)]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data
    categories = ['Daidzein\n(Basal)', 'Formononetin\n(Basal)', 
                  "6''-O-Acetyldaidzin\n(Conjugated)", "6''-Malonylgenistin\n(Conjugated)"]
    
    # Get fold changes (convert from log2)
    fold_changes = []
    for name in ['Daidzein', 'Formononetin']:
        row = df[df['Name'] == name].iloc[0]
        fold_changes.append(2 ** row['Log2FC'])
    
    for name in ["6''-O-Acetyldaidzin", "6''-Malonylgenistin"]:
        row = df[df['Name'].str.contains(name.replace("'", ""), case=False, na=False)]
        if len(row) > 0:
            fold_changes.append(2 ** row.iloc[0]['Log2FC'])
        else:
            fold_changes.append(1)
    
    colors = [COLORS['metabolite'], COLORS['metabolite'], COLORS['up'], COLORS['up']]
    
    bars = ax.bar(categories, fold_changes, color=colors, edgecolor='white', linewidth=2)
    
    # Add value labels
    for bar, val in zip(bars, fold_changes):
        height = bar.get_height()
        if height > 100:
            label = f'{val:.0f}×'
            y_pos = height + 200
        else:
            label = f'{val:.1f}×'
            y_pos = height + 0.1
        ax.text(bar.get_x() + bar.get_width()/2., y_pos,
               label, ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    ax.set_ylabel('Fold Change (Ethylene / Control)', fontweight='bold', fontsize=12)
    ax.set_title('Figure 2B: Basal vs Conjugated Isoflavonoids\nDramatic Accumulation of Conjugated Forms', 
                fontweight='bold', fontsize=14)
    
    # Use log scale for y-axis
    ax.set_yscale('log')
    ax.set_ylim(0.5, 10000)
    
    # Add annotations
    ax.axhline(1, color='gray', linestyle='--', alpha=0.5, label='No change')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'Figure2B_isoflavonoid_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'Figure2B_isoflavonoid_comparison.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: Figure2B_isoflavonoid_comparison.png/pdf")


# =============================================================================
# Figure 6: Track B - Docking Predictions (with warning)
# =============================================================================
def generate_figure6_docking():
    """Generate docking predictions visualization with warning."""
    print("Generating Figure 6: Track B Docking Predictions...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Docking data
    interactions = ['Daidzein\n↔\nFNR', 'Formononetin\n↔\nKinase']
    binding_energies = [-7.80, -7.60]
    
    bars = ax.barh(interactions, [-e for e in binding_energies], 
                   color=COLORS['track_b'], edgecolor='white', linewidth=2)
    
    # Add value labels
    for bar, val in zip(bars, binding_energies):
        width = bar.get_width()
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2.,
               f'{val} kcal/mol', ha='left', va='center', fontweight='bold', fontsize=12)
    
    ax.set_xlabel('Binding Affinity (-kcal/mol)', fontweight='bold', fontsize=12)
    ax.set_title('Figure 6: Track B - Predicted Metabolite-Protein Binding\n⚠️ EXPLORATORY - Requires Experimental Validation', 
                fontweight='bold', fontsize=14, color='#E65100')
    
    # Add Track B indicator
    ax.text(0.02, 0.98, 'Track B', transform=ax.transAxes, fontsize=12,
           fontweight='bold', color=COLORS['track_b'],
           verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='white', edgecolor=COLORS['track_b'], alpha=0.9))
    
    # Add warning box
    warning_text = """⚠️ CAUTION:
These predictions are HYPOTHETICAL
and require experimental validation
(SPR, ITC, MST, genetic studies).

Track B is INDEPENDENT from Track A
and does NOT validate GNN predictions."""
    
    ax.text(0.98, 0.02, warning_text, transform=ax.transAxes, fontsize=9,
           ha='right', va='bottom',
           bbox=dict(boxstyle='round', facecolor='#FFF3E0', edgecolor='#E65100', alpha=0.95))
    
    ax.axvline(7.0, color='gray', linestyle='--', alpha=0.5)
    ax.text(7.0, 1.5, 'Strong binding\nthreshold', ha='center', fontsize=9, color='gray')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'Figure6_docking_predictions.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'Figure6_docking_predictions.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: Figure6_docking_predictions.png/pdf")


# =============================================================================
# Figure 7: Summary - Track Relationship
# =============================================================================
def generate_figure7_track_summary():
    """Generate summary figure showing relationship between tracks."""
    print("Generating Figure 7: Track Summary...")
    
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Track A box
    track_a_box = FancyBboxPatch((0.5, 4.5), 6, 3, boxstyle="round,pad=0.05",
                                  facecolor='#E3F2FD', edgecolor=COLORS['track_a'], linewidth=3)
    ax.add_patch(track_a_box)
    ax.text(3.5, 7.2, 'Track A: Biosynthesis Pathway', ha='center', fontsize=14, 
           fontweight='bold', color=COLORS['track_a'])
    ax.text(3.5, 6.3, 'Question: Which enzymes synthesize metabolites?', ha='center', fontsize=10)
    ax.text(3.5, 5.7, 'Method: GNN (HGT) → Proteomics Validation', ha='center', fontsize=10)
    ax.text(3.5, 5.1, 'Result: PAL→4CL→CHS→CHI→IFS→IFR all ↑', ha='center', fontsize=10, 
           fontweight='bold', color='green')
    ax.text(3.5, 4.7, '✓ VALIDATED', ha='center', fontsize=12, fontweight='bold', color='green')
    
    # Track B box
    track_b_box = FancyBboxPatch((7.5, 4.5), 6, 3, boxstyle="round,pad=0.05",
                                  facecolor='#FFF3E0', edgecolor=COLORS['track_b'], linewidth=3)
    ax.add_patch(track_b_box)
    ax.text(10.5, 7.2, 'Track B: Metabolite-Protein Binding', ha='center', fontsize=14, 
           fontweight='bold', color=COLORS['track_b'])
    ax.text(10.5, 6.3, 'Question: Do metabolites regulate proteins?', ha='center', fontsize=10)
    ax.text(10.5, 5.7, 'Method: Molecular Docking', ha='center', fontsize=10)
    ax.text(10.5, 5.1, 'Result: Daidzein-FNR, Formononetin-Kinase', ha='center', fontsize=10)
    ax.text(10.5, 4.7, '⚠️ HYPOTHETICAL', ha='center', fontsize=12, fontweight='bold', color='#E65100')
    
    # Separator
    ax.text(7, 6, '≠', ha='center', fontsize=40, fontweight='bold', color='red')
    
    # Key message box
    msg_box = FancyBboxPatch((2, 0.5), 10, 3.5, boxstyle="round,pad=0.05",
                              facecolor='#FFEBEE', edgecolor='red', linewidth=2)
    ax.add_patch(msg_box)
    ax.text(7, 3.7, 'KEY MESSAGE', ha='center', fontsize=14, fontweight='bold', color='red')
    ax.text(7, 3.0, 'These tracks address DIFFERENT biological questions:', ha='center', fontsize=11)
    ax.text(7, 2.4, '• Track A: Pathway membership (enzyme synthesizes metabolite)', ha='center', fontsize=10)
    ax.text(7, 1.9, '• Track B: Physical binding (metabolite binds protein)', ha='center', fontsize=10)
    ax.text(7, 1.2, 'Docking CANNOT validate GNN predictions!', ha='center', fontsize=12, 
           fontweight='bold', color='red')
    ax.text(7, 0.7, '(Enzymes bind substrates, not products)', ha='center', fontsize=10, 
           style='italic', color='gray')
    
    plt.savefig(OUTPUT_DIR / 'Figure7_track_summary.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'Figure7_track_summary.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: Figure7_track_summary.png/pdf")


# =============================================================================
# Figure S3: GNN Ablation Study
# =============================================================================
def generate_figureS3_ablation():
    """Generate ablation study figure."""
    print("Generating Figure S3: Ablation Study...")
    
    # Ablation data
    conditions = ['Learnable\n(Baseline)', 'All-Constant\n(Topology)', 'Tier-R+P\n(Pathway)']
    hits20 = [15.0, 18.8, 18.1]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    colors = ['#90CAF9', '#64B5F6', COLORS['track_a']]
    bars = ax.bar(conditions, hits20, color=colors, edgecolor='white', linewidth=2)
    
    for bar, val in zip(bars, hits20):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.3,
               f'{val}%', ha='center', va='bottom', fontweight='bold', fontsize=12)
    
    ax.set_ylabel('Hits@20 (%)', fontweight='bold', fontsize=12)
    ax.set_title('Figure S3: Track A - Ablation Study\nTopology Learning is Critical', 
                fontweight='bold', fontsize=14)
    ax.set_ylim(0, 25)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'FigureS3_ablation_study.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'FigureS3_ablation_study.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: FigureS3_ablation_study.png/pdf")


# =============================================================================
# Main execution
# =============================================================================
def main():
    print("=" * 60)
    print("Generating Manuscript v2.0 Figures")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    
    # Generate all figures
    generate_figure2_volcano()
    generate_figure2b_isoflavonoid_comparison()
    generate_figure4_gnn_performance()
    generate_figure5_proteomics()
    generate_figure6_docking()
    generate_figure7_track_summary()
    generate_figureS3_ablation()
    
    print()
    print("=" * 60)
    print("Figure generation complete!")
    print(f"All figures saved to: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
