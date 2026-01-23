"""
Pathway Report Generator

Generates comprehensive pathway-centered final outputs for Isoflavonoid biosynthesis:
1. Enzyme ranking per pathway step (CSV)
2. Associated transcription factors (CSV)
3. Proteomics validation results (CSV)
4. Pathway summary (TXT)
5. Pathway diagram (PNG)

Usage:
    python src/generate_pathway_report.py --target Daidzein
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pathway_predictor import (
    PATHWAY_TEMPLATES,
    ENZYME_KEYWORDS,
    load_enzyme_database,
    find_enzymes_by_class,
    load_proteomics_features,
    add_proteomics_to_graph,
    train_gnn_model,
    rank_enzymes_for_metabolite,
    HGTWithFeatures,
    LinkPredictor
)

import torch

# Optional matplotlib for visualization
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available. Diagram generation will be skipped.")


def generate_pathway_enzymes_csv(pathway_result, proteomics_features, output_path):
    """Generate CSV with step-wise enzyme rankings."""
    rows = []
    
    for step_info in pathway_result:
        step = step_info['step']
        enzyme_class = step_info['enzyme_class']
        reaction = step_info.get('reaction', '')
        
        if 'top_candidates' in step_info and step_info['top_candidates']:
            for rank, cand in enumerate(step_info['top_candidates'], 1):
                prot_info = proteomics_features.get(cand['uniprot'], None)
                log2fc = prot_info[0] if prot_info else None
                pvalue = prot_info[1] if prot_info else None
                significant = pvalue is not None and pvalue < 0.05
                
                rows.append({
                    'step': step,
                    'reaction': reaction,
                    'enzyme_class': enzyme_class,
                    'rank': rank,
                    'uniprot': cand['uniprot'],
                    'name': cand['name'],
                    'gnn_score': round(cand.get('gnn_score', 0), 4),
                    'combined_score': round(cand.get('combined_score', 0), 4),
                    'log2fc': round(log2fc, 4) if log2fc is not None else '',
                    'pvalue': round(pvalue, 6) if pvalue is not None else '',
                    'significant': significant
                })
        else:
            # No candidates found
            rows.append({
                'step': step,
                'reaction': reaction,
                'enzyme_class': enzyme_class,
                'rank': '',
                'uniprot': '',
                'name': f'No {enzyme_class} candidates found',
                'gnn_score': '',
                'combined_score': '',
                'log2fc': '',
                'pvalue': '',
                'significant': ''
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"  -> Saved: {output_path}")
    return df


def generate_associated_tfs_csv(associated_tfs, proteomics_features, output_path):
    """Generate CSV with associated transcription factors."""
    rows = []
    
    for tf in associated_tfs:
        prot_info = proteomics_features.get(tf['uniprot'], None)
        log2fc = prot_info[0] if prot_info else None
        pvalue = prot_info[1] if prot_info else None
        significant = tf.get('significant', False)
        
        rows.append({
            'tf_uniprot': tf['uniprot'],
            'tf_name': tf['name'],
            'connected_enzymes': len(tf.get('connected_enzymes', [])),
            'connection_count': tf.get('connection_count', 0),
            'log2fc': round(log2fc, 4) if log2fc is not None else '',
            'pvalue': round(pvalue, 6) if pvalue is not None else '',
            'significant': significant
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"  -> Saved: {output_path}")
    return df


def generate_proteomics_validation_csv(pathway_result, proteomics_features, output_path):
    """Generate CSV validating predicted enzymes against proteomics data."""
    rows = []
    
    for step_info in pathway_result:
        if 'top_candidates' not in step_info:
            continue
            
        for rank, cand in enumerate(step_info['top_candidates'], 1):
            prot_info = proteomics_features.get(cand['uniprot'], None)
            
            if prot_info:
                log2fc = prot_info[0]
                pvalue = prot_info[1]
                fold_change = 2 ** log2fc if log2fc else None
                validated = pvalue < 0.05 if pvalue else False
                
                rows.append({
                    'enzyme_class': step_info['enzyme_class'],
                    'enzyme': cand['uniprot'],
                    'name': cand['name'],
                    'predicted_rank': rank,
                    'log2fc': round(log2fc, 4),
                    'fold_change': round(fold_change, 2) if fold_change else '',
                    'pvalue': round(pvalue, 6),
                    'validated': validated,
                    'direction': 'Up' if log2fc > 0 else 'Down'
                })
            else:
                rows.append({
                    'enzyme_class': step_info['enzyme_class'],
                    'enzyme': cand['uniprot'],
                    'name': cand['name'],
                    'predicted_rank': rank,
                    'log2fc': '',
                    'fold_change': '',
                    'pvalue': '',
                    'validated': 'No data',
                    'direction': ''
                })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"  -> Saved: {output_path}")
    return df


def generate_pathway_summary_txt(target_metabolite, pathway_result, associated_tfs, 
                                  proteomics_features, template, output_path):
    """Generate text summary of the pathway."""
    lines = []
    
    lines.append("=" * 70)
    lines.append(f"ISOFLAVONOID BIOSYNTHESIS PATHWAY REPORT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Target Metabolite: {target_metabolite}")
    lines.append(f"Pathway: {template['name']}")
    lines.append(f"Description: {template['description']}")
    lines.append("")
    
    # Biosynthetic flow
    lines.append("=" * 70)
    lines.append("BIOSYNTHETIC FLOW")
    lines.append("=" * 70)
    lines.append("")
    
    for step_info in pathway_result:
        step = step_info['step']
        substrate = template['steps'][step - 1].get('substrate', '')
        product = template['steps'][step - 1].get('product', '')
        enzyme_class = step_info['enzyme_class']
        
        lines.append(f"{substrate}")
        
        if 'top_candidates' in step_info and step_info['top_candidates']:
            best = step_info['top_candidates'][0]
            prot_info = proteomics_features.get(best['uniprot'], None)
            
            stars = ""
            log2fc_str = ""
            if prot_info:
                log2fc = prot_info[0]
                pvalue = prot_info[1]
                if pvalue < 0.05:
                    stars = " ★★" if abs(log2fc) > 1 else " ★"
                    log2fc_str = f" (Log2FC={log2fc:+.2f}, p<0.05)"
                elif pvalue < 0.1:
                    log2fc_str = f" (Log2FC={log2fc:+.2f})"
            
            lines.append(f"  ↓ [{enzyme_class}: {best['name']} ({best['uniprot']})]{stars}{log2fc_str}")
        else:
            lines.append(f"  ↓ [{enzyme_class}: ?]")
    
    # Final product
    lines.append(f"{target_metabolite}")
    lines.append("")
    
    # Regulatory TFs
    lines.append("=" * 70)
    lines.append("REGULATORY TRANSCRIPTION FACTORS")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Total TFs connected to pathway enzymes: {len(associated_tfs)}")
    lines.append("")
    
    sig_tfs = [tf for tf in associated_tfs if tf.get('significant', False)]
    if sig_tfs:
        lines.append(f"Ethylene-responsive TFs (p<0.05): {len(sig_tfs)}")
        for tf in sig_tfs[:10]:
            prot_info = proteomics_features.get(tf['uniprot'], None)
            if prot_info:
                log2fc = prot_info[0]
                direction = "↑" if log2fc > 0 else "↓"
                lines.append(f"  - {tf['name']}: {direction} {abs(log2fc):.1f}x (connects to {tf['connection_count']} enzymes)")
    else:
        lines.append("No significantly changed TFs with proteomics data")
    
    lines.append("")
    lines.append("Top connected TFs:")
    for i, tf in enumerate(associated_tfs[:10], 1):
        lines.append(f"  {i}. {tf['name']} ({tf['uniprot']}) - {tf['connection_count']} enzyme connections")
    
    # Key findings
    lines.append("")
    lines.append("=" * 70)
    lines.append("KEY FINDINGS")
    lines.append("=" * 70)
    lines.append("")
    
    # Count significant enzymes
    total_steps = 0
    sig_steps = 0
    for step_info in pathway_result:
        if 'top_candidates' in step_info and step_info['top_candidates']:
            total_steps += 1
            best = step_info['top_candidates'][0]
            prot_info = proteomics_features.get(best['uniprot'], None)
            if prot_info and prot_info[1] < 0.05:
                sig_steps += 1
    
    lines.append(f"• {sig_steps}/{total_steps} pathway enzymes significantly regulated by ethylene")
    
    # Find strongest induction
    strongest = None
    strongest_fc = 0
    for step_info in pathway_result:
        if 'top_candidates' in step_info and step_info['top_candidates']:
            best = step_info['top_candidates'][0]
            prot_info = proteomics_features.get(best['uniprot'], None)
            if prot_info and prot_info[1] < 0.05:
                fc = abs(prot_info[0])
                if fc > strongest_fc:
                    strongest_fc = fc
                    strongest = (step_info['enzyme_class'], best['name'], prot_info[0])
    
    if strongest:
        direction = "upregulated" if strongest[2] > 0 else "downregulated"
        fold = 2 ** abs(strongest[2])
        lines.append(f"• {strongest[0]} ({strongest[1]}) shows strongest induction: {fold:.0f}× fold change")
    
    lines.append("")
    lines.append("=" * 70)
    lines.append("END OF REPORT")
    lines.append("=" * 70)
    
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"  -> Saved: {output_path}")
    return '\n'.join(lines)


def generate_pathway_diagram(pathway_result, proteomics_features, template, 
                              target_metabolite, output_path):
    """Generate clean vertical flowchart-style pathway diagram."""
    if not HAS_MATPLOTLIB:
        print("  -> Skipped diagram (matplotlib not available)")
        return
    
    n_steps = len(pathway_result)
    fig_height = max(12, 2 + n_steps * 1.8)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, fig_height))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, fig_height)
    ax.axis('off')
    ax.set_aspect('equal')
    
    # Colors
    COLORS = {
        'substrate': '#E3F2FD',
        'substrate_edge': '#1976D2',
        'product': '#C8E6C9',
        'product_edge': '#388E3C',
        'enzyme_nodata': '#ECEFF1',
        'enzyme_edge': '#455A64',
        'up_strong': '#EF5350',      # Red
        'up_weak': '#FFCDD2',        # Light red
        'down_strong': '#42A5F5',    # Blue  
        'down_weak': '#BBDEFB',      # Light blue
        'not_sig': '#FFF9C4',        # Light yellow
        'arrow': '#546E7A'
    }
    
    # Title
    title_y = fig_height - 0.5
    ax.text(5, title_y, f'Isoflavonoid Biosynthesis Pathway',
            ha='center', va='center', fontsize=14, fontweight='bold')
    ax.text(5, title_y - 0.4, f'Target: {target_metabolite} | Ethylene Response',
            ha='center', va='center', fontsize=10, color='#616161')
    
    # Draw pathway as vertical flowchart
    center_x = 5
    box_width = 3.5
    box_height = 0.55
    enzyme_width = 2.8
    enzyme_height = 0.5
    y_spacing = 1.5
    
    start_y = fig_height - 2.0
    
    for i, step_info in enumerate(pathway_result):
        y = start_y - i * y_spacing
        step = step_info['step']
        enzyme_class = step_info['enzyme_class']
        
        substrate = template['steps'][step - 1].get('substrate', f'Step {step}')
        
        # ===== SUBSTRATE BOX (centered) =====
        sub_box = FancyBboxPatch(
            (center_x - box_width/2, y - box_height/2), 
            box_width, box_height,
            boxstyle="round,pad=0.02,rounding_size=0.15",
            facecolor=COLORS['substrate'], 
            edgecolor=COLORS['substrate_edge'], 
            linewidth=2
        )
        ax.add_patch(sub_box)
        ax.text(center_x, y, substrate[:30], 
                ha='center', va='center', fontsize=9, fontweight='bold')
        
        # ===== DOWN ARROW =====
        arrow_y_start = y - box_height/2 - 0.05
        arrow_y_end = y - y_spacing/2 + enzyme_height/2 + 0.1
        ax.annotate('', 
                   xy=(center_x, arrow_y_end), 
                   xytext=(center_x, arrow_y_start),
                   arrowprops=dict(arrowstyle='-|>', color=COLORS['arrow'], lw=2))
        
        # ===== ENZYME BOX (centered, between arrows) =====
        enzyme_y = y - y_spacing/2
        
        if 'top_candidates' in step_info and step_info['top_candidates']:
            best = step_info['top_candidates'][0]
            prot_info = proteomics_features.get(best['uniprot'], None)
            
            # Determine color based on expression
            if prot_info and prot_info[1] < 0.05:
                log2fc = prot_info[0]
                if log2fc > 1:
                    enz_color = COLORS['up_strong']
                    edge_color = '#C62828'
                elif log2fc > 0:
                    enz_color = COLORS['up_weak']
                    edge_color = '#E57373'
                elif log2fc < -1:
                    enz_color = COLORS['down_strong']
                    edge_color = '#1565C0'
                else:
                    enz_color = COLORS['down_weak']
                    edge_color = '#64B5F6'
                fc_text = f"({log2fc:+.1f}*)"
            elif prot_info:
                enz_color = COLORS['not_sig']
                edge_color = '#FBC02D'
                fc_text = f"({prot_info[0]:+.1f})"
            else:
                enz_color = COLORS['enzyme_nodata']
                edge_color = COLORS['enzyme_edge']
                fc_text = ""
            
            enz_label = f"{enzyme_class}: {best['name'][:12]}"
        else:
            enz_color = COLORS['enzyme_nodata']
            edge_color = COLORS['enzyme_edge']
            enz_label = f"{enzyme_class}: ?"
            fc_text = ""
        
        enz_box = FancyBboxPatch(
            (center_x - enzyme_width/2, enzyme_y - enzyme_height/2), 
            enzyme_width, enzyme_height,
            boxstyle="round,pad=0.02,rounding_size=0.1",
            facecolor=enz_color, 
            edgecolor=edge_color, 
            linewidth=1.5
        )
        ax.add_patch(enz_box)
        
        display_text = f"{enz_label} {fc_text}".strip()
        ax.text(center_x, enzyme_y, display_text,
                ha='center', va='center', fontsize=8, fontweight='bold')
        
        # ===== DOWN ARROW to next substrate =====
        if i < n_steps - 1:
            arrow_y_start2 = enzyme_y - enzyme_height/2 - 0.05
            arrow_y_end2 = y - y_spacing + box_height/2 + 0.05
            ax.annotate('', 
                       xy=(center_x, arrow_y_end2), 
                       xytext=(center_x, arrow_y_start2),
                       arrowprops=dict(arrowstyle='-|>', color=COLORS['arrow'], lw=2))
    
    # ===== FINAL PRODUCT BOX =====
    final_y = start_y - n_steps * y_spacing + y_spacing/2
    final_box = FancyBboxPatch(
        (center_x - box_width/2, final_y - box_height/2), 
        box_width, box_height,
        boxstyle="round,pad=0.02,rounding_size=0.15",
        facecolor=COLORS['product'], 
        edgecolor=COLORS['product_edge'], 
        linewidth=2.5
    )
    ax.add_patch(final_box)
    ax.text(center_x, final_y, f"★ {target_metabolite} ★", 
            ha='center', va='center', fontsize=10, fontweight='bold', color='#1B5E20')
    
    # ===== LEGEND =====
    legend_y = 0.8
    legend_items = [
        (COLORS['up_strong'], 'Up (|Log2FC|>1, p<0.05)'),
        (COLORS['up_weak'], 'Up (Log2FC>0, p<0.05)'),
        (COLORS['down_strong'], 'Down (|Log2FC|>1, p<0.05)'),
        (COLORS['down_weak'], 'Down (Log2FC<0, p<0.05)'),
        (COLORS['not_sig'], 'Not Significant'),
        (COLORS['enzyme_nodata'], 'No Data'),
    ]
    
    legend_x_start = 0.5
    for i, (color, label) in enumerate(legend_items):
        x = legend_x_start + (i % 3) * 3.2
        y = legend_y - (i // 3) * 0.35
        rect = plt.Rectangle((x, y - 0.1), 0.25, 0.2, facecolor=color, edgecolor='#424242', lw=0.5)
        ax.add_patch(rect)
        ax.text(x + 0.35, y, label, fontsize=7, va='center')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  -> Saved: {output_path}")


def run_pathway_analysis(target_metabolite='Daidzein', pathway_type='isoflavonoid', 
                         top_k=3, output_dir='results/pathway'):
    """Run complete pathway analysis and generate all outputs."""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\n{'='*70}")
    print(f"PATHWAY REPORT GENERATOR")
    print(f"Target: {target_metabolite} | Pathway: {pathway_type}")
    print(f"{'='*70}\n")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print("[1/7] Loading graph data...")
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)
    
    print("[2/7] Loading mappings...")
    enzyme_mapping_df = pd.read_csv('data/processed/enzyme_string_mapping.csv')
    tf_mapping_df = pd.read_csv('data/processed/tf_string_mapping.csv')
    
    uniprot_to_idx = dict(zip(enzyme_mapping_df['uniprot_id'], enzyme_mapping_df['enzyme_idx']))
    idx_to_uniprot = dict(zip(enzyme_mapping_df['enzyme_idx'], enzyme_mapping_df['uniprot_id']))
    idx_to_tf_uniprot = dict(zip(tf_mapping_df['tf_idx'], tf_mapping_df['uniprot_id']))
    tf_idx_to_name = dict(zip(tf_mapping_df['tf_idx'], tf_mapping_df['uniprot_id']))
    
    enzyme_db = load_enzyme_database()
    
    print("[3/7] Loading proteomics features...")
    proteomics_features = load_proteomics_features('data/processed/pxd006989_mapped.csv')
    data = add_proteomics_to_graph(data, proteomics_features, enzyme_mapping_df)
    
    print("[4/7] Training GNN model...")
    model, predictor = train_gnn_model(data, device)
    
    # Get pathway template
    template = PATHWAY_TEMPLATES.get(pathway_type)
    if not template:
        print(f"ERROR: Unknown pathway type: {pathway_type}")
        return
    
    print(f"[5/7] Predicting pathway enzymes for {template['name']}...")
    
    # Predict each step
    pathway_result = []
    pathway_enzyme_indices = []
    
    for step_info in template['steps']:
        step_num = step_info['step']
        enzyme_class = step_info['enzyme_class']
        enzyme_name = step_info['enzyme_name']
        reaction = step_info['reaction']
        
        # Find candidate enzymes
        keywords = ENZYME_KEYWORDS.get(enzyme_class, [enzyme_name])
        candidates = find_enzymes_by_class(enzyme_db, enzyme_class, keywords)
        
        if not candidates:
            pathway_result.append({
                'step': step_num,
                'reaction': reaction,
                'enzyme_class': enzyme_class,
                'status': 'no_candidates'
            })
            continue
        
        # Get indices for candidates in graph
        candidate_data = []
        for uniprot, info, match_type in candidates:
            if uniprot in uniprot_to_idx:
                idx = uniprot_to_idx[uniprot]
                prot_info = proteomics_features.get(uniprot, None)
                
                candidate_data.append({
                    'uniprot': uniprot,
                    'idx': idx,
                    'name': info['name'],
                    'log2fc': prot_info[0] if prot_info else None,
                    'pvalue': prot_info[1] if prot_info else None
                })
        
        if not candidate_data:
            pathway_result.append({
                'step': step_num,
                'reaction': reaction,
                'enzyme_class': enzyme_class,
                'status': 'not_in_graph'
            })
            continue
        
        # Rank by GNN
        candidate_indices = [c['idx'] for c in candidate_data]
        gnn_scores = rank_enzymes_for_metabolite(model, predictor, data, 0, candidate_indices, device)
        
        for i, cand in enumerate(candidate_data):
            cand['gnn_score'] = gnn_scores[i]
            prot_boost = 0
            if cand['log2fc'] is not None and cand['pvalue'] is not None:
                if cand['pvalue'] < 0.05:
                    prot_boost = 0.2 * abs(cand['log2fc'])
            cand['combined_score'] = cand['gnn_score'] + prot_boost
        
        candidate_data.sort(key=lambda x: x['combined_score'], reverse=True)
        
        for cand in candidate_data[:top_k]:
            pathway_enzyme_indices.append(cand['idx'])
        
        pathway_result.append({
            'step': step_num,
            'reaction': reaction,
            'enzyme_class': enzyme_class,
            'top_candidates': candidate_data[:top_k],
            'all_candidates': len(candidate_data),
            'status': 'predicted'
        })
    
    # Find associated TFs
    print("[6/7] Finding associated transcription factors...")
    
    tf_enz_edges = data['TF', 'interacts', 'Enzyme'].edge_index
    interacting_tfs = {}
    
    for enz_idx in pathway_enzyme_indices:
        mask = tf_enz_edges[1] == enz_idx
        tf_indices = tf_enz_edges[0][mask].tolist()
        
        for tf_idx in tf_indices:
            if tf_idx not in interacting_tfs:
                tf_uniprot = idx_to_tf_uniprot.get(tf_idx, f'TF_{tf_idx}')
                prot_info = proteomics_features.get(tf_uniprot, None)
                
                interacting_tfs[tf_idx] = {
                    'idx': tf_idx,
                    'uniprot': tf_uniprot,
                    'name': tf_idx_to_name.get(tf_idx, tf_uniprot),
                    'connected_enzymes': [],
                    'connection_count': 0,
                    'log2fc': prot_info[0] if prot_info else None,
                    'pvalue': prot_info[1] if prot_info else None,
                    'significant': prot_info[1] < 0.05 if prot_info and prot_info[1] else False
                }
            interacting_tfs[tf_idx]['connected_enzymes'].append(enz_idx)
            interacting_tfs[tf_idx]['connection_count'] += 1
    
    associated_tfs = list(interacting_tfs.values())
    associated_tfs.sort(key=lambda x: (x['significant'], x['connection_count'],
                                        abs(x['log2fc']) if x['log2fc'] else 0), reverse=True)
    
    # Generate outputs
    print("[7/7] Generating output files...")
    
    # 1. Enzyme ranking CSV
    enzyme_csv_path = os.path.join(output_dir, 'isoflavonoid_pathway_enzymes.csv')
    generate_pathway_enzymes_csv(pathway_result, proteomics_features, enzyme_csv_path)
    
    # 2. Associated TFs CSV
    tf_csv_path = os.path.join(output_dir, 'pathway_associated_tfs.csv')
    generate_associated_tfs_csv(associated_tfs, proteomics_features, tf_csv_path)
    
    # 3. Proteomics validation CSV
    validation_csv_path = os.path.join(output_dir, 'pathway_proteomics_validation.csv')
    generate_proteomics_validation_csv(pathway_result, proteomics_features, validation_csv_path)
    
    # 4. Pathway summary TXT
    summary_path = os.path.join(output_dir, 'pathway_summary.txt')
    generate_pathway_summary_txt(target_metabolite, pathway_result, associated_tfs,
                                  proteomics_features, template, summary_path)
    
    # 5. Pathway diagram PNG
    diagram_path = os.path.join(output_dir, 'pathway_diagram.png')
    generate_pathway_diagram(pathway_result, proteomics_features, template,
                             target_metabolite, diagram_path)
    
    print(f"\n{'='*70}")
    print("REPORT GENERATION COMPLETE")
    print(f"{'='*70}")
    print(f"\nOutput files saved to: {output_dir}/")
    print(f"  - isoflavonoid_pathway_enzymes.csv")
    print(f"  - pathway_associated_tfs.csv")
    print(f"  - pathway_proteomics_validation.csv")
    print(f"  - pathway_summary.txt")
    print(f"  - pathway_diagram.png")
    
    return {
        'pathway': pathway_result,
        'associated_tfs': associated_tfs,
        'output_dir': output_dir
    }


def main():
    parser = argparse.ArgumentParser(description='Generate Pathway Report')
    parser.add_argument('--target', type=str, default='Daidzein',
                        help='Target metabolite name (default: Daidzein)')
    parser.add_argument('--pathway', type=str, default='isoflavonoid',
                        choices=['isoflavonoid', 'phaseollin'],
                        help='Pathway type (default: isoflavonoid)')
    parser.add_argument('--top-k', type=int, default=3,
                        help='Number of top enzyme isoforms per step (default: 3)')
    parser.add_argument('--output-dir', type=str, default='results/pathway',
                        help='Output directory (default: results/pathway)')
    
    args = parser.parse_args()
    
    run_pathway_analysis(
        target_metabolite=args.target,
        pathway_type=args.pathway,
        top_k=args.top_k,
        output_dir=args.output_dir
    )


if __name__ == '__main__':
    main()
