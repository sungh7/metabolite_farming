#!/usr/bin/env python
"""
Model Evaluation Runner Script

Executes the complete HGT attention analysis pipeline implementing
recommendations from the GNN vs Transformer critical review.

Usage:
    python -m src.model_evaluation.run_evaluation [--output-dir PATH]

Output:
    - Attention analysis reports
    - Cross-attention visualizations
    - Explanation reports
    - Summary JSON files
"""

import torch
import os
import sys
import argparse
import json
from datetime import datetime

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.model import HGT, LinkPredictor
from src.model_evaluation.attention_extractor import (
    generate_attention_report,
    compute_cross_type_attention
)
from src.model_evaluation.cross_attention import (
    generate_cross_attention_report
)
from src.model_evaluation.attention_explainer import (
    generate_explanation_report
)


def load_data_and_models(data_path: str, device: torch.device):
    """Load graph data and initialize models."""
    print(f"Loading data from {data_path}...")
    data = torch.load(data_path)
    data = data.to(device)

    print(f"  Node types: {data.metadata()[0]}")
    print(f"  Edge types: {len(data.metadata()[1])} types")

    # Get feature dimension
    sample_type = list(data.x_dict.keys())[0]
    in_dim = data.x_dict[sample_type].size(1)

    print(f"  Feature dimension: {in_dim}")

    # Initialize models
    print("\nInitializing HGT model...")
    model = HGT(
        metadata=data.metadata(),
        in_channels=in_dim,
        hidden_channels=64,
        out_channels=64,
        num_heads=4,
        num_layers=2
    )
    predictor = LinkPredictor(64)

    model.to(device)
    predictor.to(device)

    return data, model, predictor


def train_model(data, model, predictor, device, epochs: int = 20):
    """Train model for evaluation."""
    print(f"\nTraining model for {epochs} epochs...")

    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=0.01
    )

    model.train()
    predictor.train()

    # Find edge type for training
    train_edge_type = data.edge_types[0]
    for et in data.edge_types:
        if 'TF' in et[0] or 'Enzyme' in et[0]:
            train_edge_type = et
            break

    src_type, _, dst_type = train_edge_type
    print(f"  Training on edge type: {train_edge_type}")

    for epoch in range(epochs):
        optimizer.zero_grad()

        x_dict = model(data.x_dict, data.edge_index_dict)

        pos_edge = data[train_edge_type].edge_index

        # Contrastive loss
        pos_out = predictor(x_dict[src_type], x_dict[dst_type], pos_edge)

        # Negative sampling
        num_neg = pos_edge.size(1)
        neg_dst = torch.randint(0, x_dict[dst_type].size(0), (num_neg,), device=device)
        neg_edge = torch.stack([pos_edge[0], neg_dst])
        neg_out = predictor(x_dict[src_type], x_dict[dst_type], neg_edge)

        loss = -torch.log(torch.sigmoid(pos_out) + 1e-8).mean() - \
               torch.log(1 - torch.sigmoid(neg_out) + 1e-8).mean()

        loss.backward()
        optimizer.step()

        if (epoch + 1) % 5 == 0:
            print(f"    Epoch {epoch+1}/{epochs}, Loss: {loss.item():.4f}")

    return model, predictor


def run_full_evaluation(data_path: str, output_dir: str):
    """
    Run complete model evaluation pipeline.

    Implements recommendations from GNN vs Transformer review:
    1. HGT attention analysis
    2. Cross-attention (Proteomics ↔ Metabolomics)
    3. Attention-weighted explanations
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # Create output directories
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_base = os.path.join(output_dir, f"evaluation_{timestamp}")
    os.makedirs(output_base, exist_ok=True)

    print(f"\nOutput directory: {output_base}")

    # Load data and models
    data, model, predictor = load_data_and_models(data_path, device)

    # Train model
    model, predictor = train_model(data, model, predictor, device, epochs=20)

    # Analysis pipeline
    results = {}

    # 1. HGT Attention Analysis
    print("\n" + "="*60)
    print("STEP 1: HGT Attention Analysis")
    print("="*60)
    attention_dir = os.path.join(output_base, "attention_analysis")
    try:
        attention_summary = generate_attention_report(
            data, model, predictor, device, attention_dir
        )
        results['attention_analysis'] = 'completed'
    except Exception as e:
        print(f"Warning: Attention analysis failed: {e}")
        results['attention_analysis'] = f'failed: {str(e)}'

    # 2. Cross-Attention Analysis
    print("\n" + "="*60)
    print("STEP 2: Cross-Attention Analysis (Proteomics ↔ Metabolomics)")
    print("="*60)
    cross_attn_dir = os.path.join(output_base, "cross_attention")
    try:
        # Try common node type combinations
        for protein_type in ['Enzyme', 'Protein', 'TF']:
            for metab_type in ['Metabolite', 'Enzyme']:
                if protein_type != metab_type:
                    if protein_type in data.x_dict and metab_type in data.x_dict:
                        print(f"  Analyzing: {protein_type} ↔ {metab_type}")
                        cross_attn_scores = generate_cross_attention_report(
                            data, model, device, cross_attn_dir,
                            protein_type=protein_type,
                            metabolite_type=metab_type
                        )
                        results['cross_attention'] = 'completed'
                        break
            else:
                continue
            break
        else:
            print("  No suitable node type pairs found for cross-attention")
            results['cross_attention'] = 'no_suitable_types'
    except Exception as e:
        print(f"Warning: Cross-attention analysis failed: {e}")
        results['cross_attention'] = f'failed: {str(e)}'

    # 3. Explanation Analysis
    print("\n" + "="*60)
    print("STEP 3: Attention-Weighted Path Explanations")
    print("="*60)
    explain_dir = os.path.join(output_base, "explanations")
    try:
        # Try TF -> Enzyme explanations
        for src_type in ['TF', 'Signaling', 'Protein']:
            for dst_type in ['Enzyme', 'Protein', 'Metabolite']:
                if src_type != dst_type:
                    if src_type in data.x_dict and dst_type in data.x_dict:
                        print(f"  Explaining: {src_type} → {dst_type}")
                        explanations = generate_explanation_report(
                            data, model, predictor, device, explain_dir,
                            src_type=src_type, dst_type=dst_type, top_k=10
                        )
                        results['explanations'] = 'completed'
                        break
            else:
                continue
            break
        else:
            print("  No suitable node type pairs found for explanations")
            results['explanations'] = 'no_suitable_types'
    except Exception as e:
        print(f"Warning: Explanation analysis failed: {e}")
        results['explanations'] = f'failed: {str(e)}'

    # Save overall summary
    summary = {
        'timestamp': timestamp,
        'device': str(device),
        'data_path': data_path,
        'output_dir': output_base,
        'node_types': list(data.metadata()[0]),
        'edge_types': [str(et) for et in data.metadata()[1]],
        'results': results,
        'notes': {
            'data_type': 'Static comparison (72h single time point)',
            'interpretation': 'Attention scores indicate association strength, not causality',
            'recommendation': 'For causal claims, perturbation experiments required'
        }
    }

    with open(os.path.join(output_base, 'evaluation_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print("\n" + "="*60)
    print("EVALUATION COMPLETE")
    print("="*60)
    print(f"\nResults saved to: {output_base}/")
    print("\nSubdirectories:")
    print("  - attention_analysis/: HGT attention weights and heatmaps")
    print("  - cross_attention/: Proteomics-Metabolomics cross-attention")
    print("  - explanations/: Attention-weighted path explanations")
    print("\nKey output files:")
    print("  - evaluation_summary.json: Overall summary")
    print("  - */cross_type_attention_heatmap.png: Node type interactions")
    print("  - */top_cross_attention_pairs.csv: Top protein-metabolite pairs")
    print("  - */explanation_summary.csv: Top TF-Enzyme predictions")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Run HGT model evaluation implementing GNN vs Transformer review recommendations'
    )
    parser.add_argument(
        '--data-path',
        type=str,
        default='data/processed/graph.pt',
        help='Path to graph data file'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='results/model_evaluation',
        help='Output directory for results'
    )
    args = parser.parse_args()

    # Check for alternative data paths
    if not os.path.exists(args.data_path):
        alt_paths = [
            'data/processed/strict_bipartite_v2.pt',
            'data/processed/strict_bipartite.pt',
            'data/processed/hetero_graph.pt'
        ]
        for alt in alt_paths:
            if os.path.exists(alt):
                print(f"Using alternative data path: {alt}")
                args.data_path = alt
                break
        else:
            print(f"Error: Data file not found at {args.data_path}")
            print(f"Tried alternatives: {alt_paths}")
            sys.exit(1)

    run_full_evaluation(args.data_path, args.output_dir)


if __name__ == "__main__":
    main()
