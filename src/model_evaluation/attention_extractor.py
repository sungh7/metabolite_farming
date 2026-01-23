"""
HGT Attention Extraction and Visualization Module

Extracts attention weights from Heterogeneous Graph Transformer (HGT)
for interpretable analysis of TF-Enzyme regulatory relationships.

Based on plan review:
- Current data is static (72h single time point)
- Focus on edge-level attention interpretation
- HGT already combines GNN structure with Transformer attention
"""

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
from typing import Dict, Optional
from torch_geometric.nn import HGTConv
import warnings
warnings.filterwarnings('ignore')


class AttentionHGT(nn.Module):
    """
    HGT model with attention weight extraction capability.

    Wraps standard HGT to capture layer-wise attention scores
    for interpretability analysis.
    """

    def __init__(self, metadata, in_channels, hidden_channels, out_channels,
                 num_heads=4, num_layers=2):
        super().__init__()
        self.metadata = metadata
        self.num_heads = num_heads
        self.num_layers = num_layers

        # Node type projections
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = nn.Linear(in_channels, hidden_channels)

        # HGT convolution layers
        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HGTConv(hidden_channels, hidden_channels, metadata, heads=num_heads)
            self.convs.append(conv)

        self.out_lin = nn.Linear(hidden_channels, out_channels)

        # Storage for attention weights
        self.attention_weights = {}
        self._register_hooks()

    def _register_hooks(self):
        """Register forward hooks to capture attention weights."""
        self._hooks = []
        for i, conv in enumerate(self.convs):
            hook = conv.register_forward_hook(
                lambda module, _input, _output, layer_idx=i:
                self._capture_attention(module, layer_idx)
            )
            self._hooks.append(hook)

    def _capture_attention(self, module, layer_idx):
        """Capture attention weights from HGTConv layer."""
        if hasattr(module, 'alpha'):
            self.attention_weights[f'layer_{layer_idx}'] = {
                'alpha': module.alpha.detach().cpu() if module.alpha is not None else None
            }

    def forward(self, x_dict, edge_index_dict, return_attention=False):
        """Forward pass with optional attention extraction."""
        self.attention_weights = {}

        # Initial projection
        x_dict = {
            node_type: torch.relu(self.lin_dict[node_type](x))
            for node_type, x in x_dict.items()
        }

        # HGT convolution layers
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)

        if return_attention:
            return x_dict, self.attention_weights
        return x_dict

    def get_edge_attention(self, data, src_type: str, src_idx: int,
                           dst_type: str, dst_idx: int) -> Dict:
        """Get attention scores for edges involving specific nodes."""
        self.eval()
        with torch.no_grad():
            x_dict, attn = self.forward(data.x_dict, data.edge_index_dict,
                                        return_attention=True)

        edge_attention = {}

        for edge_type in data.edge_types:
            et_src, rel, et_dst = edge_type
            edge_index = data[edge_type].edge_index

            if et_src == src_type:
                mask = edge_index[0] == src_idx
                if mask.any():
                    neighbors = edge_index[1, mask].tolist()
                    edge_attention[f'{src_type}->{et_dst}:{rel}'] = {
                        'neighbors': neighbors,
                        'neighbor_type': et_dst
                    }

            if et_dst == src_type:
                mask = edge_index[1] == src_idx
                if mask.any():
                    neighbors = edge_index[0, mask].tolist()
                    edge_attention[f'{et_src}->{src_type}:{rel}'] = {
                        'neighbors': neighbors,
                        'neighbor_type': et_src
                    }

        return {
            'source': {'type': src_type, 'idx': src_idx},
            'target': {'type': dst_type, 'idx': dst_idx},
            'edge_connections': edge_attention,
            'layer_attention': attn
        }


def compute_cross_type_attention(data, model, device) -> pd.DataFrame:
    """
    Compute cross-type attention matrix showing how different node types
    attend to each other through the HGT layers.

    Returns:
        DataFrame with attention flow between node types
    """
    model.eval()
    model.to(device)
    data = data.to(device)

    node_types = data.metadata()[0]

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)

    # Compute cross-type similarity as proxy for attention
    attention_matrix = np.zeros((len(node_types), len(node_types)))

    for i, src_type in enumerate(node_types):
        for j, dst_type in enumerate(node_types):
            if src_type in x_dict and dst_type in x_dict:
                src_emb = x_dict[src_type]
                dst_emb = x_dict[dst_type]

                src_norm = src_emb / (src_emb.norm(dim=1, keepdim=True) + 1e-8)
                dst_norm = dst_emb / (dst_emb.norm(dim=1, keepdim=True) + 1e-8)

                n_src = min(100, src_norm.size(0))
                n_dst = min(100, dst_norm.size(0))

                src_sample = src_norm[:n_src]
                dst_sample = dst_norm[:n_dst]

                sim = torch.mm(src_sample, dst_sample.t())
                attention_matrix[i, j] = sim.mean().item()

    return pd.DataFrame(
        attention_matrix,
        index=node_types,
        columns=node_types
    )


def compute_edge_importance_scores(data, model, predictor,
                                    src_type: str, dst_type: str,
                                    device, top_k: int = 100) -> pd.DataFrame:
    """
    Compute importance scores for edges of a given type.

    Based on plan: "Attention 해석 → 어떤 TF-Enzyme edge가 중요한가 시각화"
    """
    model.eval()
    predictor.eval()
    model.to(device)
    predictor.to(device)
    data = data.to(device)

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)

    src_emb = x_dict[src_type]
    dst_emb = x_dict[dst_type]

    n_src = src_emb.size(0)
    n_dst = dst_emb.size(0)

    # Sample if too large
    max_pairs = 10000
    if n_src * n_dst > max_pairs:
        src_indices = torch.randperm(n_src)[:min(n_src, 100)]
        dst_indices = torch.randperm(n_dst)[:min(n_dst, 100)]
    else:
        src_indices = torch.arange(n_src)
        dst_indices = torch.arange(n_dst)

    results = []

    for src_idx in src_indices:
        for dst_idx in dst_indices:
            edge = torch.tensor([[src_idx], [dst_idx]], device=device)
            score = predictor(src_emb, dst_emb, edge).sigmoid().item()
            results.append({
                'src_type': src_type,
                'src_idx': src_idx.item(),
                'dst_type': dst_type,
                'dst_idx': dst_idx.item(),
                'score': score
            })

    df = pd.DataFrame(results)
    df = df.sort_values('score', ascending=False).head(top_k)

    return df


def visualize_attention_heatmap(attention_df: pd.DataFrame,
                                 output_path: str,
                                 title: str = "Cross-Type Attention Matrix"):
    """Create heatmap visualization of cross-type attention."""
    fig, ax = plt.subplots(figsize=(10, 8))

    sns.heatmap(
        attention_df,
        annot=True,
        fmt='.3f',
        cmap='YlOrRd',
        ax=ax,
        square=True,
        cbar_kws={'label': 'Average Attention Score'}
    )

    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('Target Node Type', fontsize=12)
    ax.set_ylabel('Source Node Type', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved attention heatmap to {output_path}")


def visualize_tf_enzyme_importance(importance_df: pd.DataFrame,
                                   output_path: str,
                                   node_names: Optional[Dict] = None,
                                   title: str = "TF-Enzyme Edge Importance"):
    """Create visualization of TF-Enzyme edge importance scores."""
    fig, ax = plt.subplots(figsize=(12, 8))

    top_df = importance_df.head(20).copy()

    if node_names:
        top_df['label'] = top_df.apply(
            lambda r: f"{node_names.get(('TF', r['src_idx']), f'TF_{r.src_idx}')} → "
                     f"{node_names.get(('Enzyme', r['dst_idx']), f'Enz_{r.dst_idx}')}",
            axis=1
        )
    else:
        top_df['label'] = top_df.apply(
            lambda r: f"TF_{r['src_idx']} → Enzyme_{r['dst_idx']}",
            axis=1
        )

    cmap = plt.get_cmap('viridis')
    colors = [cmap(x) for x in np.linspace(0.2, 0.8, len(top_df))]

    bars = ax.barh(range(len(top_df)), top_df['score'], color=colors)
    ax.set_yticks(range(len(top_df)))
    ax.set_yticklabels(top_df['label'])
    ax.invert_yaxis()

    ax.set_xlabel('Prediction Score', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    for bar, score in zip(bars, top_df['score']):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                f'{score:.3f}', va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved TF-Enzyme importance plot to {output_path}")


def generate_attention_report(data, model, predictor, device, output_dir: str):
    """
    Generate comprehensive attention analysis report.

    Implements the plan's recommendations for HGT interpretation enhancement.
    """
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("HGT ATTENTION ANALYSIS REPORT")
    print("=" * 60)

    # 1. Cross-type attention matrix
    print("\n1. Computing cross-type attention matrix...")
    attention_df = compute_cross_type_attention(data, model, device)
    attention_df.to_csv(os.path.join(output_dir, 'cross_type_attention.csv'))

    visualize_attention_heatmap(
        attention_df,
        os.path.join(output_dir, 'cross_type_attention_heatmap.png'),
        title='HGT Cross-Type Attention: Node Type Interactions'
    )

    print(f"\nCross-Type Attention Matrix:\n{attention_df.round(3)}")

    # 2. TF-Enzyme importance
    print("\n2. Computing TF-Enzyme edge importance...")
    tf_enz_importance = compute_edge_importance_scores(
        data, model, predictor, 'TF', 'Enzyme', device, top_k=100
    )
    tf_enz_importance.to_csv(
        os.path.join(output_dir, 'tf_enzyme_importance.csv'),
        index=False
    )

    visualize_tf_enzyme_importance(
        tf_enz_importance,
        os.path.join(output_dir, 'tf_enzyme_importance.png'),
        title='Top TF-Enzyme Regulatory Predictions'
    )

    # 3. Node type statistics
    print("\n3. Computing node type statistics...")
    node_stats = {}
    for node_type in data.metadata()[0]:
        if node_type in data.x_dict:
            x = data.x_dict[node_type]
            node_stats[node_type] = {
                'count': x.size(0),
                'feature_dim': x.size(1),
                'mean_norm': x.norm(dim=1).mean().item()
            }

    stats_df = pd.DataFrame(node_stats).T
    stats_df.to_csv(os.path.join(output_dir, 'node_type_statistics.csv'))

    print(f"\nNode Type Statistics:\n{stats_df.round(3)}")

    # 4. Edge type statistics
    print("\n4. Computing edge type statistics...")
    edge_stats = {}
    for edge_type in data.metadata()[1]:
        edge_index = data[edge_type].edge_index
        edge_stats[str(edge_type)] = {
            'num_edges': edge_index.size(1),
            'src_type': edge_type[0],
            'relation': edge_type[1],
            'dst_type': edge_type[2]
        }

    edge_df = pd.DataFrame(edge_stats).T
    edge_df.to_csv(os.path.join(output_dir, 'edge_type_statistics.csv'))

    print(f"\nEdge Type Statistics:\n{edge_df}")

    # 5. Summary report
    summary = {
        'analysis_type': 'HGT Attention Analysis',
        'note': 'Static comparison (72h single time point)',
        'model_architecture': 'Heterogeneous Graph Transformer (HGT)',
        'node_types': list(data.metadata()[0]),
        'edge_types': [str(et) for et in data.metadata()[1]],
        'cross_type_attention': attention_df.to_dict(),
        'top_tf_enzyme_pairs': tf_enz_importance.head(10).to_dict('records'),
        'node_statistics': node_stats
    }

    with open(os.path.join(output_dir, 'attention_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n{'='*60}")
    print(f"Report saved to {output_dir}/")
    print(f"{'='*60}")

    return summary
