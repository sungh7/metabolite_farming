#!/usr/bin/env python3
"""
Create supplementary figure showing database coverage and metabolite categories.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

def create_coverage_figure():
    df = pd.read_csv('data/processed/mtbls531_differential_enhanced.csv')

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=300)

    # Panel A: Database mapping pie chart
    ax = axes[0, 0]
    mapped = df['KEGG'].notna().sum()
    unmapped = df['KEGG'].isna().sum()

    ax.pie([mapped, unmapped], labels=['Mapped to KEGG\n(n=29, 36.7%)',
                                       'Specialized metabolites\n(n=50, 63.3%)'],
           colors=[COLORS['blue'], COLORS['grey']], autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax.set_title('A. KEGG Database Coverage', fontsize=13, fontweight='bold', pad=15)

    # Panel B: Significance by mapping status
    ax = axes[0, 1]

    mapped_sig = df[(df['KEGG'].notna()) & (df['P_Value'] < 0.05)].shape[0]
    mapped_ns = df[(df['KEGG'].notna()) & (df['P_Value'] >= 0.05)].shape[0]
    unmapped_sig = df[(df['KEGG'].isna()) & (df['P_Value'] < 0.05)].shape[0]
    unmapped_ns = df[(df['KEGG'].isna()) & (df['P_Value'] >= 0.05)].shape[0]

    x = np.arange(2)
    width = 0.35

    bars1 = ax.bar(x - width/2, [mapped_sig, unmapped_sig], width,
                   label='Significant (P<0.05)', color=COLORS['green'],
                   edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, [mapped_ns, unmapped_ns], width,
                   label='Not significant', color=COLORS['grey'],
                   edgecolor='black', linewidth=1.5)

    ax.set_ylabel('Number of metabolites', fontsize=11, fontweight='bold')
    ax.set_title('B. Significance by Mapping Status', fontsize=13, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(['Mapped to KEGG', 'Specialized'], fontweight='bold')
    ax.legend(fontsize=10, frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3, axis='y')

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Panel C: Top unmapped metabolites
    ax = axes[1, 0]

    unmapped_sig_df = df[(df['KEGG'].isna()) & (df['P_Value'] < 0.05)].sort_values('P_Value').head(10)

    y_pos = np.arange(len(unmapped_sig_df))
    colors_c = [COLORS['red'] if fc > 10 else COLORS['orange'] if fc > 5 else COLORS['yellow']
                for fc in unmapped_sig_df['Log2FC']]

    ax.barh(y_pos, -np.log10(unmapped_sig_df['P_Value']), color=colors_c,
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    labels = [name[:35] + '...' if len(name) > 35 else name
              for name in unmapped_sig_df['Name']]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('-log10(P-value)', fontsize=11, fontweight='bold')
    ax.set_title('C. Top Unmapped Significant Metabolites', fontsize=13, fontweight='bold', pad=15)
    ax.axvline(-np.log10(0.05), color='black', linestyle='--', linewidth=2, alpha=0.5)
    ax.grid(True, alpha=0.3, axis='x')

    # Panel D: Fold change distribution
    ax = axes[1, 1]

    mapped_fc = df[df['KEGG'].notna()]['Log2FC']
    unmapped_fc = df[df['KEGG'].isna()]['Log2FC']

    bins = np.linspace(-15, 15, 30)
    ax.hist(mapped_fc, bins=bins, alpha=0.7, label='Mapped to KEGG',
            color=COLORS['blue'], edgecolor='black', linewidth=0.5)
    ax.hist(unmapped_fc, bins=bins, alpha=0.7, label='Specialized metabolites',
            color=COLORS['orange'], edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Log2 Fold Change', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('D. Fold Change Distribution by Mapping', fontsize=13, fontweight='bold', pad=15)
    ax.legend(fontsize=10, frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3, axis='y')
    ax.axvline(0, color='black', linestyle='-', linewidth=1)

    plt.suptitle('Database Coverage Analysis: KEGG Mapping and Metabolite Categories',
                fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()

    plt.savefig('results/figures/pathway_analysis/figureS1_database_coverage.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('results/figures/pathway_analysis/figureS1_database_coverage.pdf',
                bbox_inches='tight', facecolor='white')
    plt.close()

    print('✓ Created Figure S1: Database coverage analysis')
    print('  Files: figureS1_database_coverage.png/pdf')

if __name__ == "__main__":
    create_coverage_figure()
