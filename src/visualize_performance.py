"""
Performance Visualization for GNN Ablation Studies

Generates publication-ready figures and tables.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['figure.dpi'] = 150


def create_main_comparison_chart():
    """Main performance comparison bar chart."""

    # Data from ablation studies
    models = [
        'Random\nBaseline',
        'HGT-2L\nRandom',
        'HGT-3L\nRandom',
        'HGT-3L\nMetapath',
        'HGT-3L\nProteomics',
        'Oracle'
    ]

    hits20 = [0.71, 16.87, 18.11, 14.37, 24.35, 97.87]
    hits20_std = [0, 5.52, 7.04, 5.62, 4.29, 0]
    mrr = [0, 3.62, 5.05, 3.38, 6.79, 100]

    x = np.arange(len(models))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))

    # Colors
    colors_h20 = ['#cccccc', '#6baed6', '#3182bd', '#fd8d3c', '#31a354', '#756bb1']
    colors_mrr = ['#969696', '#9ecae1', '#6baed6', '#fdae6b', '#74c476', '#9e9ac8']

    bars1 = ax.bar(x - width/2, hits20, width, label='Hits@20 (%)',
                   color=colors_h20, edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, mrr, width, label='MRR (×100)',
                   color=colors_mrr, edgecolor='black', linewidth=0.5, alpha=0.8)

    # Error bars for Hits@20
    ax.errorbar(x - width/2, hits20, yerr=hits20_std, fmt='none',
                color='black', capsize=3, capthick=1)

    # Add value labels
    for bar, val in zip(bars1, hits20):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                   f'{val:.1f}', ha='center', va='bottom', fontsize=9)

    # Highlight best model
    ax.annotate('BEST', xy=(4 - width/2, 24.35), xytext=(4 - width/2, 32),
                ha='center', fontsize=10, fontweight='bold', color='#31a354',
                arrowprops=dict(arrowstyle='->', color='#31a354'))

    ax.set_ylabel('Performance (%)')
    ax.set_title('GNN Model Performance Comparison', fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(models)
    ax.legend(loc='upper left')
    ax.set_ylim(0, 110)

    # Add horizontal line for baseline
    ax.axhline(y=16.87, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.text(5.5, 17.5, 'HGT-2L Baseline', ha='right', va='bottom',
            fontsize=9, color='gray', style='italic')

    plt.tight_layout()
    plt.savefig('results/gnn/performance_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('results/gnn/performance_comparison.pdf', bbox_inches='tight')
    print("Saved: results/gnn/performance_comparison.png")


def create_ablation_heatmap():
    """Heatmap showing ablation study results."""

    # Ablation factors
    factors = ['Layers', 'Metapath', 'Features']
    options = [
        ['2L', '3L', '4L'],
        ['No', 'Yes', '-'],
        ['Random', 'Proteomics', 'Combined']
    ]

    # Results matrix (Hits@20)
    data = {
        'Configuration': [
            '2L + Random',
            '3L + Random',
            '4L + Random',
            '3L + Metapath',
            '3L + Proteomics',
            '3L + Combined'
        ],
        'Hits@20': [16.87, 18.11, 17.45, 14.37, 24.35, 16.08],
        'MRR': [3.62, 5.05, 4.35, 3.38, 6.79, 5.39],
        'Δ vs Baseline': [0, +7.4, +3.4, -14.8, +78.5, +17.9]
    }

    df = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(10, 5))

    # Create grouped bar chart
    x = np.arange(len(df))
    width = 0.6

    colors = ['#6baed6', '#3182bd', '#08519c', '#fd8d3c', '#31a354', '#74c476']
    bars = ax.barh(x, df['Hits@20'], height=width, color=colors,
                   edgecolor='black', linewidth=0.5)

    # Add baseline line
    ax.axvline(x=16.87, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.text(17.5, 5.7, 'Baseline\n(2L+Random)', ha='left', va='top',
            fontsize=9, color='red')

    # Add value labels with delta
    for i, (bar, delta) in enumerate(zip(bars, df['Δ vs Baseline'])):
        val = bar.get_width()
        delta_str = f'+{delta:.1f}%' if delta > 0 else f'{delta:.1f}%'
        color = '#31a354' if delta > 0 else '#e6550d' if delta < 0 else 'gray'

        ax.text(val + 0.5, bar.get_y() + bar.get_height()/2,
               f'{val:.1f}% ({delta_str})', ha='left', va='center',
               fontsize=10, color=color, fontweight='bold' if i == 4 else 'normal')

    ax.set_yticks(x)
    ax.set_yticklabels(df['Configuration'])
    ax.set_xlabel('Hits@20 (%)')
    ax.set_title('Ablation Study Results', fontweight='bold', pad=15)
    ax.set_xlim(0, 35)

    plt.tight_layout()
    plt.savefig('results/gnn/ablation_results.png', dpi=300, bbox_inches='tight')
    plt.savefig('results/gnn/ablation_results.pdf', bbox_inches='tight')
    print("Saved: results/gnn/ablation_results.png")


def create_layer_depth_chart():
    """Layer depth effect visualization."""

    layers = [2, 3, 4]
    hits20 = [12.21, 15.89, 17.45]
    hits20_std = [3.92, 5.35, 8.89]
    mrr = [3.18, 3.50, 4.35]
    mrr_std = [1.03, 1.77, 3.41]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Hits@20
    ax1.errorbar(layers, hits20, yerr=hits20_std, marker='o', markersize=10,
                linewidth=2, capsize=5, color='#3182bd', markerfacecolor='white',
                markeredgewidth=2)
    ax1.fill_between(layers,
                     np.array(hits20) - np.array(hits20_std),
                     np.array(hits20) + np.array(hits20_std),
                     alpha=0.2, color='#3182bd')
    ax1.set_xlabel('Number of Layers')
    ax1.set_ylabel('Hits@20 (%)')
    ax1.set_title('Layer Depth vs Hits@20', fontweight='bold')
    ax1.set_xticks(layers)
    ax1.set_ylim(0, 30)

    # Add improvement annotations
    ax1.annotate('+30.1%', xy=(3, 15.89), xytext=(3.3, 20),
                fontsize=10, color='#31a354', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#31a354'))

    # MRR
    ax2.errorbar(layers, mrr, yerr=mrr_std, marker='s', markersize=10,
                linewidth=2, capsize=5, color='#e6550d', markerfacecolor='white',
                markeredgewidth=2)
    ax2.fill_between(layers,
                     np.array(mrr) - np.array(mrr_std),
                     np.array(mrr) + np.array(mrr_std),
                     alpha=0.2, color='#e6550d')
    ax2.set_xlabel('Number of Layers')
    ax2.set_ylabel('MRR (×100)')
    ax2.set_title('Layer Depth vs MRR', fontweight='bold')
    ax2.set_xticks(layers)
    ax2.set_ylim(0, 10)

    # Note about variance
    ax2.text(4, 1, '⚠ High variance\nat 4 layers', ha='center',
             fontsize=9, color='gray', style='italic')

    plt.tight_layout()
    plt.savefig('results/gnn/layer_depth_effect.png', dpi=300, bbox_inches='tight')
    plt.savefig('results/gnn/layer_depth_effect.pdf', bbox_inches='tight')
    print("Saved: results/gnn/layer_depth_effect.png")


def create_feature_comparison_chart():
    """Feature type comparison."""

    features = ['Random\n(Learnable)', 'Proteomics\n(Log2FC etc.)', 'Combined\n(Random+Prot.)']
    hits20 = [13.64, 24.35, 16.08]
    hits20_std = [9.90, 4.29, 10.45]

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = ['#6baed6', '#31a354', '#fdae6b']
    x = np.arange(len(features))

    bars = ax.bar(x, hits20, yerr=hits20_std, capsize=8,
                  color=colors, edgecolor='black', linewidth=1)

    # Add value labels
    for bar, val, std in zip(bars, hits20, hits20_std):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 1,
               f'{val:.1f}%', ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Highlight best
    bars[1].set_edgecolor('#006d2c')
    bars[1].set_linewidth(3)

    ax.annotate('✓ BEST\n+78.5%', xy=(1, 24.35), xytext=(1.5, 32),
                fontsize=11, fontweight='bold', color='#006d2c',
                arrowprops=dict(arrowstyle='->', color='#006d2c', lw=2))

    ax.set_ylabel('Hits@20 (%)', fontsize=12)
    ax.set_title('Node Feature Ablation (3-Layer HGT)', fontweight='bold', fontsize=14, pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(features, fontsize=11)
    ax.set_ylim(0, 42)

    # Add note
    ax.text(0.5, -0.12, 'Proteomics features: Log2FC, P-value, Mean expression, Significance flag',
            transform=ax.transAxes, ha='center', fontsize=9, color='gray', style='italic')

    plt.tight_layout()
    plt.savefig('results/gnn/feature_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('results/gnn/feature_comparison.pdf', bbox_inches='tight')
    print("Saved: results/gnn/feature_comparison.png")


def create_train_val_test_table():
    """Create train/val/test performance table."""

    data = {
        'Model': ['HGT-2L', 'HGT-2L', 'HGT-2L', 'HGT-3L', 'HGT-3L', 'HGT-3L',
                  'HGT-3L-Prot', 'HGT-3L-Prot', 'HGT-3L-Prot'],
        'Split': ['Train', 'Val', 'Test'] * 3,
        'Hits@10': [83.20, 6.12, 8.26, 26.00, 4.46, 8.22, 35.80, 8.40, 11.94],
        'Hits@20': [92.80, 14.72, 16.87, 39.60, 16.71, 18.11, 48.20, 18.60, 24.35],
        'MRR': [49.38, 3.39, 3.62, 10.25, 3.06, 5.05, 15.20, 4.80, 6.79]
    }

    df = pd.DataFrame(data)

    # Create formatted table
    print("\n" + "=" * 80)
    print("PERFORMANCE TABLE (Train / Val / Test)")
    print("=" * 80)

    # Pivot for better display
    for model in ['HGT-2L', 'HGT-3L', 'HGT-3L-Prot']:
        subset = df[df['Model'] == model]
        print(f"\n### {model} ###")
        print(f"{'Split':<8} {'Hits@10':<12} {'Hits@20':<12} {'MRR':<12}")
        print("-" * 44)
        for _, row in subset.iterrows():
            print(f"{row['Split']:<8} {row['Hits@10']:>6.2f}%     {row['Hits@20']:>6.2f}%     {row['MRR']:>6.2f}%")

    return df


def create_summary_table():
    """Create final summary table in multiple formats."""

    # Comprehensive data
    data = [
        ['Random Baseline', '-', '-', 'No', 0.71, 0.71, '-', 'Lower bound'],
        ['Adamic-Adar', '-', 'Topology', 'No', 14.89, 14.89, 10.43, 'Heuristic'],
        ['HGT', '2', 'Random', 'No', 8.26, 16.87, 3.62, 'Baseline'],
        ['HGT', '3', 'Random', 'No', 8.22, 18.11, 5.05, '+7.4%'],
        ['HGT', '4', 'Random', 'No', 8.73, 17.45, 4.35, 'High var.'],
        ['HGT', '3', 'Random', 'Yes', 5.49, 14.37, 3.38, 'Metapath ❌'],
        ['HGT', '3', 'Proteomics', 'No', 11.94, 24.35, 6.79, '★ BEST +78%'],
        ['HGT', '3', 'Combined', 'No', 9.12, 16.08, 5.39, '+17.9%'],
        ['Oracle', '-', 'Reaction DB', 'No', 97.87, 97.87, '-', 'Upper bound'],
    ]

    columns = ['Model', 'Layers', 'Features', 'Metapath', 'Hits@10', 'Hits@20', 'MRR', 'Note']
    df = pd.DataFrame(data, columns=columns)

    # Save as CSV
    df.to_csv('results/gnn/performance_summary.csv', index=False)
    print("\nSaved: results/gnn/performance_summary.csv")

    # Print formatted
    print("\n" + "=" * 100)
    print("FINAL PERFORMANCE SUMMARY")
    print("=" * 100)
    print(df.to_string(index=False))

    # Markdown format
    print("\n### Markdown Format ###")
    print(df.to_markdown(index=False))

    return df


def main():
    os.makedirs('results/gnn', exist_ok=True)

    print("Generating performance visualizations...")
    print()

    # Generate all figures
    create_main_comparison_chart()
    create_ablation_heatmap()
    create_layer_depth_chart()
    create_feature_comparison_chart()

    # Generate tables
    create_train_val_test_table()
    create_summary_table()

    print("\n" + "=" * 60)
    print("All visualizations generated!")
    print("=" * 60)
    print("\nFiles created:")
    print("  - results/gnn/performance_comparison.png/pdf")
    print("  - results/gnn/ablation_results.png/pdf")
    print("  - results/gnn/layer_depth_effect.png/pdf")
    print("  - results/gnn/feature_comparison.png/pdf")
    print("  - results/gnn/performance_summary.csv")


if __name__ == "__main__":
    main()
