"""
Generate Figures for Independent Validation

Creates:
- Figure 1: Module Enrichment Summary (bar plot)
- Figure 2: Rank Agreement (heatmap)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import sys
sys.path.append(os.getcwd())

def load_validation_metrics():
    """Load validation metrics from TSV."""
    df = pd.read_csv('results/independent_validation/validation_metrics.tsv', sep='\t')
    return df

def create_enrichment_figure(df, output_path):
    """
    Figure 1: Module Enrichment Summary
    Bar plot showing enrichment scores by dataset and module.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    datasets = df['Dataset'].unique()
    modules = df['Module'].unique()
    
    x = np.arange(len(datasets))
    width = 0.25
    
    colors = {'phenylpropanoid': '#E63946', 'flavonoid': '#457B9D', 'amino_acid': '#2A9D8F'}
    
    for i, module in enumerate(modules):
        module_data = df[df['Module'] == module]
        scores = [module_data[module_data['Dataset'] == d]['Enrichment_Score'].values[0] 
                  for d in datasets]
        bars = ax.bar(x + i * width, scores, width, label=module.replace('_', ' ').title(),
                      color=colors.get(module, '#666666'))
    
    ax.set_xlabel('Dataset', fontsize=12)
    ax.set_ylabel('Enrichment Score', fontsize=12)
    ax.set_title('Module-Level Enrichment in Independent Datasets', fontsize=14, fontweight='bold')
    ax.set_xticks(x + width)
    ax.set_xticklabels(datasets, fontsize=10)
    ax.legend(title='Module', loc='upper right')
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Expected')
    ax.set_ylim(0, max(df['Enrichment_Score'].max() * 1.2, 2))
    
    # Add significance markers
    for i, module in enumerate(modules):
        module_data = df[df['Module'] == module]
        for j, dataset in enumerate(datasets):
            row = module_data[module_data['Dataset'] == dataset]
            if len(row) > 0 and row['N_Significant'].values[0] > 0:
                score = row['Enrichment_Score'].values[0]
                ax.annotate('*', (x[j] + i * width, score + 0.05), 
                           ha='center', fontsize=12, color='red')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

def create_activity_heatmap(df, output_path):
    """
    Figure 2: Module Activity Heatmap
    Shows which modules are active in each dataset.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    
    datasets = df['Dataset'].unique()
    modules = df['Module'].unique()
    
    # Create matrix
    matrix = np.zeros((len(modules), len(datasets)))
    for i, module in enumerate(modules):
        for j, dataset in enumerate(datasets):
            row = df[(df['Module'] == module) & (df['Dataset'] == dataset)]
            if len(row) > 0:
                matrix[i, j] = row['N_Significant'].values[0]
    
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
    
    ax.set_xticks(np.arange(len(datasets)))
    ax.set_yticks(np.arange(len(modules)))
    ax.set_xticklabels(datasets, fontsize=10)
    ax.set_yticklabels([m.replace('_', ' ').title() for m in modules], fontsize=10)
    
    # Add text annotations
    for i in range(len(modules)):
        for j in range(len(datasets)):
            text = ax.text(j, i, int(matrix[i, j]),
                          ha='center', va='center', color='black' if matrix[i, j] < 5 else 'white',
                          fontsize=11, fontweight='bold')
    
    ax.set_title('Number of Significant Features per Module', fontsize=14, fontweight='bold')
    
    # Colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('N Significant', rotation=-90, va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

def create_performance_comparison(output_path):
    """
    Figure 3: Model Performance Comparison
    Bar plot comparing different model variants.
    """
    # Performance data
    models = ['Random', 'AA/RA', 'HGT\n(Hard Neg)', 'HGT\n(Frozen)', 'HGT\n(Attribute)', 
              'HGT\n(Inductive)', 'Oracle']
    hits20 = [0.71, 14.89, 15.03, 15.03, 13.55, 20.77, 97.87]
    stds = [0, 0, 3.87, 3.87, 4.24, 9.43, 0]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = ['#E0E0E0', '#90BE6D', '#F9844A', '#F9844A', '#F9844A', '#43AA8B', '#577590']
    
    bars = ax.bar(models, hits20, yerr=stds, capsize=5, color=colors, edgecolor='black')
    
    ax.set_ylabel('Hits@20 (%)', fontsize=12)
    ax.set_title('Link Prediction Performance Comparison', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 110)
    
    # Add value labels
    for bar, val, std in zip(bars, hits20, stds):
        height = bar.get_height()
        label = f'{val:.1f}%' if std == 0 else f'{val:.1f}±{std:.1f}%'
        ax.annotate(label,
                   xy=(bar.get_x() + bar.get_width() / 2, height + std + 1),
                   ha='center', va='bottom', fontsize=9)
    
    # Add reference lines
    ax.axhline(y=0.71, color='gray', linestyle='--', alpha=0.3)
    ax.axhline(y=15.03, color='orange', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

def main():
    os.makedirs('results/independent_validation', exist_ok=True)
    
    print("Generating Figures...")
    
    # Load data
    df = load_validation_metrics()
    
    # Figure 1: Enrichment
    create_enrichment_figure(df, 'results/independent_validation/figure1_enrichment.png')
    
    # Figure 2: Activity Heatmap
    create_activity_heatmap(df, 'results/independent_validation/figure2_activity_heatmap.png')
    
    # Figure 3: Performance Comparison
    create_performance_comparison('results/independent_validation/figure3_performance.png')
    
    print("\nAll figures generated!")

if __name__ == "__main__":
    main()
