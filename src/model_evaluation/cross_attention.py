"""
Cross-Attention Visualization for Multi-omics Integration

Implements cross-attention analysis between Proteomics and Metabolomics data
as recommended in the GNN vs Transformer review plan.

Key insight from plan review:
- Current data is static (72h single time point)
- Cross-attention shows evidence strength, not causality
- Serves as "Proteomics evidence ↔ Metabolomics result" direct linking
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
from typing import Dict, Optional


class CrossAttentionLayer(nn.Module):
    """
    Cross-attention layer for multi-omics data fusion.

    Computes bidirectional attention between two omics modalities
    (e.g., Proteomics and Metabolomics).
    """

    def __init__(self, embed_dim: int, num_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.head_dim = embed_dim // num_heads

        assert embed_dim % num_heads == 0, "embed_dim must be divisible by num_heads"

        # Query, Key, Value projections for modality A -> B attention
        self.q_proj_a = nn.Linear(embed_dim, embed_dim)
        self.k_proj_b = nn.Linear(embed_dim, embed_dim)
        self.v_proj_b = nn.Linear(embed_dim, embed_dim)

        # Query, Key, Value projections for modality B -> A attention
        self.q_proj_b = nn.Linear(embed_dim, embed_dim)
        self.k_proj_a = nn.Linear(embed_dim, embed_dim)
        self.v_proj_a = nn.Linear(embed_dim, embed_dim)

        self.dropout = nn.Dropout(dropout)
        self.scale = self.head_dim ** -0.5

        self.out_proj_a = nn.Linear(embed_dim, embed_dim)
        self.out_proj_b = nn.Linear(embed_dim, embed_dim)

        self.attn_weights_a_to_b = None
        self.attn_weights_b_to_a = None

    def forward(self, x_a: torch.Tensor, x_b: torch.Tensor,
                return_attention: bool = False):
        """Compute bidirectional cross-attention."""
        if x_a.dim() == 2:
            x_a = x_a.unsqueeze(0)
        if x_b.dim() == 2:
            x_b = x_b.unsqueeze(0)

        batch_size = x_a.size(0)
        seq_len_a = x_a.size(1)
        seq_len_b = x_b.size(1)

        # A attends to B
        q_a = self.q_proj_a(x_a).view(batch_size, seq_len_a, self.num_heads, self.head_dim).transpose(1, 2)
        k_b = self.k_proj_b(x_b).view(batch_size, seq_len_b, self.num_heads, self.head_dim).transpose(1, 2)
        v_b = self.v_proj_b(x_b).view(batch_size, seq_len_b, self.num_heads, self.head_dim).transpose(1, 2)

        attn_a_to_b = torch.matmul(q_a, k_b.transpose(-2, -1)) * self.scale
        attn_a_to_b = F.softmax(attn_a_to_b, dim=-1)
        self.attn_weights_a_to_b = attn_a_to_b.detach()
        attn_a_to_b = self.dropout(attn_a_to_b)

        out_a = torch.matmul(attn_a_to_b, v_b)
        out_a = out_a.transpose(1, 2).contiguous().view(batch_size, seq_len_a, self.embed_dim)
        out_a = self.out_proj_a(out_a)

        # B attends to A
        q_b = self.q_proj_b(x_b).view(batch_size, seq_len_b, self.num_heads, self.head_dim).transpose(1, 2)
        k_a = self.k_proj_a(x_a).view(batch_size, seq_len_a, self.num_heads, self.head_dim).transpose(1, 2)
        v_a = self.v_proj_a(x_a).view(batch_size, seq_len_a, self.num_heads, self.head_dim).transpose(1, 2)

        attn_b_to_a = torch.matmul(q_b, k_a.transpose(-2, -1)) * self.scale
        attn_b_to_a = F.softmax(attn_b_to_a, dim=-1)
        self.attn_weights_b_to_a = attn_b_to_a.detach()
        attn_b_to_a = self.dropout(attn_b_to_a)

        out_b = torch.matmul(attn_b_to_a, v_a)
        out_b = out_b.transpose(1, 2).contiguous().view(batch_size, seq_len_b, self.embed_dim)
        out_b = self.out_proj_b(out_b)

        out_a = out_a.squeeze(0)
        out_b = out_b.squeeze(0)

        if return_attention:
            return out_a, out_b, self.attn_weights_a_to_b, self.attn_weights_b_to_a
        return out_a, out_b


class MultiOmicsCrossAttention(nn.Module):
    """
    Multi-omics integration using cross-attention.

    Designed for Proteomics-Metabolomics fusion in the context of
    ethylene-induced isoflavonoid biosynthesis.
    """

    def __init__(self, embed_dim: int, num_heads: int = 4, num_layers: int = 2):
        super().__init__()
        self.layers = nn.ModuleList([
            CrossAttentionLayer(embed_dim, num_heads)
            for _ in range(num_layers)
        ])
        self.layer_norm_prot = nn.LayerNorm(embed_dim)
        self.layer_norm_metab = nn.LayerNorm(embed_dim)

    def forward(self, prot_emb: torch.Tensor, metab_emb: torch.Tensor,
                return_all_attention: bool = False):
        """Forward pass through cross-attention layers."""
        all_attention = []

        for layer in self.layers:
            prot_out, metab_out, attn_p2m, attn_m2p = layer(
                prot_emb, metab_emb, return_attention=True
            )

            prot_emb = self.layer_norm_prot(prot_emb + prot_out)
            metab_emb = self.layer_norm_metab(metab_emb + metab_out)

            if return_all_attention:
                all_attention.append({
                    'protein_to_metabolite': attn_p2m,
                    'metabolite_to_protein': attn_m2p
                })

        if return_all_attention:
            return prot_emb, metab_emb, all_attention
        return prot_emb, metab_emb


def compute_cross_attention_scores(data, model, device,
                                    protein_type: str = 'Enzyme',
                                    metabolite_type: str = 'Metabolite') -> Dict:
    """
    Compute cross-attention scores between proteins and metabolites.

    Uses HGT embeddings as input to cross-attention analysis.
    """
    model.eval()
    model.to(device)
    data = data.to(device)

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)

    prot_emb = x_dict.get(protein_type)
    metab_emb = x_dict.get(metabolite_type)

    if prot_emb is None or metab_emb is None:
        return {'error': f'Missing embeddings for {protein_type} or {metabolite_type}'}

    embed_dim = prot_emb.size(1)
    cross_attn = MultiOmicsCrossAttention(embed_dim, num_heads=4, num_layers=1)
    cross_attn.to(device)
    cross_attn.eval()

    with torch.no_grad():
        _, _, all_attn = cross_attn(
            prot_emb, metab_emb, return_all_attention=True
        )

    attn_p2m = all_attn[0]['protein_to_metabolite'].squeeze(0).mean(dim=0)
    attn_m2p = all_attn[0]['metabolite_to_protein'].squeeze(0).mean(dim=0)

    n_prot = prot_emb.size(0)
    n_metab = metab_emb.size(0)

    attn_flat = attn_p2m.flatten()
    top_k = min(100, attn_flat.numel())
    top_values, top_indices = torch.topk(attn_flat, top_k)

    top_pairs = []
    for idx, val in zip(top_indices.cpu().numpy(), top_values.cpu().numpy()):
        prot_idx = idx // n_metab
        metab_idx = idx % n_metab
        top_pairs.append({
            'protein_idx': int(prot_idx),
            'metabolite_idx': int(metab_idx),
            'attention_score': float(val)
        })

    return {
        'protein_to_metabolite_attention': attn_p2m.cpu().numpy(),
        'metabolite_to_protein_attention': attn_m2p.cpu().numpy(),
        'top_pairs': top_pairs,
        'n_proteins': n_prot,
        'n_metabolites': n_metab,
    }


def visualize_cross_attention(attention_scores: Dict, output_dir: str,
                               protein_names: Optional[Dict] = None,
                               metabolite_names: Optional[Dict] = None,
                               max_display: int = 30):
    """Create visualizations of cross-attention analysis."""
    os.makedirs(output_dir, exist_ok=True)

    attn_p2m = attention_scores['protein_to_metabolite_attention']
    attn_m2p = attention_scores['metabolite_to_protein_attention']
    top_pairs = attention_scores['top_pairs']

    # 1. Cross-attention heatmaps
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    n_prot = min(max_display, attn_p2m.shape[0])
    n_metab = min(max_display, attn_p2m.shape[1])

    prot_importance = attn_p2m.sum(axis=1)
    metab_importance = attn_p2m.sum(axis=0)

    top_prot_idx = np.argsort(prot_importance)[-n_prot:]
    top_metab_idx = np.argsort(metab_importance)[-n_metab:]

    attn_subset = attn_p2m[np.ix_(top_prot_idx, top_metab_idx)]

    if protein_names:
        prot_labels = [protein_names.get(i, f'P_{i}') for i in top_prot_idx]
    else:
        prot_labels = [f'P_{i}' for i in top_prot_idx]

    if metabolite_names:
        metab_labels = [metabolite_names.get(i, f'M_{i}') for i in top_metab_idx]
    else:
        metab_labels = [f'M_{i}' for i in top_metab_idx]

    sns.heatmap(attn_subset, ax=axes[0], cmap='YlOrRd',
                xticklabels=metab_labels, yticklabels=prot_labels)
    axes[0].set_title('Protein → Metabolite Attention', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Metabolite', fontsize=10)
    axes[0].set_ylabel('Protein/Enzyme', fontsize=10)
    axes[0].tick_params(axis='both', labelsize=7)

    attn_m2p_subset = attn_m2p[np.ix_(top_metab_idx, top_prot_idx)]
    sns.heatmap(attn_m2p_subset, ax=axes[1], cmap='YlGnBu',
                xticklabels=prot_labels, yticklabels=metab_labels)
    axes[1].set_title('Metabolite → Protein Attention', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Protein/Enzyme', fontsize=10)
    axes[1].set_ylabel('Metabolite', fontsize=10)
    axes[1].tick_params(axis='both', labelsize=7)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cross_attention_heatmap.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # 2. Top pairs bar chart
    fig, ax = plt.subplots(figsize=(12, 8))

    top_20 = top_pairs[:20]
    labels = []
    scores = []

    for pair in top_20:
        p_name = protein_names.get(pair['protein_idx'], f"P_{pair['protein_idx']}") if protein_names else f"P_{pair['protein_idx']}"
        m_name = metabolite_names.get(pair['metabolite_idx'], f"M_{pair['metabolite_idx']}") if metabolite_names else f"M_{pair['metabolite_idx']}"
        labels.append(f"{p_name} → {m_name}")
        scores.append(pair['attention_score'])

    cmap = plt.get_cmap('viridis')
    colors = [cmap(x) for x in np.linspace(0.3, 0.9, len(labels))]

    ax.barh(range(len(labels)), scores, color=colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel('Cross-Attention Score', fontsize=12)
    ax.set_title('Top Protein-Metabolite Cross-Attention Pairs', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'top_cross_attention_pairs.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # 3. Attention distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(attn_p2m.flatten(), bins=50, color='coral', alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('Attention Score', fontsize=11)
    axes[0].set_ylabel('Frequency', fontsize=11)
    axes[0].set_title('Protein → Metabolite Attention Distribution', fontsize=12)
    axes[0].axvline(np.mean(attn_p2m), color='red', linestyle='--', label=f'Mean: {np.mean(attn_p2m):.4f}')
    axes[0].legend()

    axes[1].hist(attn_m2p.flatten(), bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    axes[1].set_xlabel('Attention Score', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_title('Metabolite → Protein Attention Distribution', fontsize=12)
    axes[1].axvline(np.mean(attn_m2p), color='blue', linestyle='--', label=f'Mean: {np.mean(attn_m2p):.4f}')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'attention_distribution.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Cross-attention visualizations saved to {output_dir}/")


def generate_cross_attention_report(data, model, device, output_dir: str,
                                     protein_type: str = 'Enzyme',
                                     metabolite_type: str = 'Metabolite'):
    """
    Generate comprehensive cross-attention analysis report.

    Implements: "Cross-attention layer → Proteomics 증거 ↔ Metabolomics 결과 직접 연결"
    """
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("CROSS-ATTENTION ANALYSIS: Proteomics ↔ Metabolomics")
    print("=" * 60)
    print("\nNote: Current data is static (72h single time point).")
    print("Cross-attention shows evidence strength, not causality.\n")

    print("Computing cross-attention scores...")
    scores = compute_cross_attention_scores(
        data, model, device, protein_type, metabolite_type
    )

    if 'error' in scores:
        print(f"Error: {scores['error']}")
        return scores

    print("Generating visualizations...")
    visualize_cross_attention(scores, output_dir)

    top_pairs_df = pd.DataFrame(scores['top_pairs'])
    top_pairs_df.to_csv(os.path.join(output_dir, 'top_cross_attention_pairs.csv'), index=False)

    summary = {
        'analysis': 'Cross-Attention Analysis',
        'data_type': 'Static (72h single time point)',
        'interpretation': 'Association strength, not causality',
        'n_proteins': scores['n_proteins'],
        'n_metabolites': scores['n_metabolites'],
        'attention_stats': {
            'p2m_mean': float(np.mean(scores['protein_to_metabolite_attention'])),
            'p2m_std': float(np.std(scores['protein_to_metabolite_attention'])),
            'p2m_max': float(np.max(scores['protein_to_metabolite_attention'])),
            'm2p_mean': float(np.mean(scores['metabolite_to_protein_attention'])),
            'm2p_std': float(np.std(scores['metabolite_to_protein_attention'])),
            'm2p_max': float(np.max(scores['metabolite_to_protein_attention']))
        },
        'top_10_pairs': scores['top_pairs'][:10]
    }

    with open(os.path.join(output_dir, 'cross_attention_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nSummary:")
    print(f"  - Proteins analyzed: {scores['n_proteins']}")
    print(f"  - Metabolites analyzed: {scores['n_metabolites']}")
    print(f"  - Mean P→M attention: {summary['attention_stats']['p2m_mean']:.4f}")
    print(f"  - Mean M→P attention: {summary['attention_stats']['m2p_mean']:.4f}")

    print(f"\nTop 5 Protein-Metabolite pairs by attention:")
    for i, pair in enumerate(scores['top_pairs'][:5], 1):
        print(f"  {i}. Protein[{pair['protein_idx']}] → Metabolite[{pair['metabolite_idx']}]: {pair['attention_score']:.4f}")

    print(f"\nResults saved to {output_dir}/")

    return scores
