"""
Model Evaluation Package for GNN vs Transformer Analysis

This package implements the recommendations from the critical review of
GNN vs Transformer vs Causal Model approaches for isoflavonoid biosynthesis research.

Key findings from review:
1. Current data is static (72h single time point) - no time-series
2. HGT (Heterogeneous Graph Transformer) already combines GNN + Transformer attention
3. Focus on attention interpretation rather than temporal modeling
4. Express findings as "association" not "causality" without perturbation data

Modules:
- attention_extractor: HGT attention weight extraction and visualization
- cross_attention: Multi-omics cross-attention analysis (Proteomics ↔ Metabolomics)
- attention_explainer: Attention-weighted path explanations
- review_report: Analysis documentation generator
"""

from .attention_extractor import (
    AttentionHGT,
    compute_cross_type_attention,
    compute_edge_importance_scores,
    visualize_attention_heatmap,
    generate_attention_report
)

from .cross_attention import (
    CrossAttentionLayer,
    MultiOmicsCrossAttention,
    compute_cross_attention_scores,
    generate_cross_attention_report
)

from .attention_explainer import (
    AttentionWeightedExplainer,
    explain_top_predictions,
    visualize_attention_paths
)

__all__ = [
    'AttentionHGT',
    'compute_cross_type_attention',
    'compute_edge_importance_scores',
    'visualize_attention_heatmap',
    'generate_attention_report',
    'CrossAttentionLayer',
    'MultiOmicsCrossAttention',
    'compute_cross_attention_scores',
    'generate_cross_attention_report',
    'AttentionWeightedExplainer',
    'explain_top_predictions',
    'visualize_attention_paths'
]
