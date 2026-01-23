"""
Attention-Weighted Path Explainer

Combines BFS path finding with attention weights to provide
interpretable explanations of TF-Enzyme regulatory predictions.

Based on plan review:
- "Attention 해석 → 어떤 TF-Enzyme edge가 중요한가 시각화"
- Focus on association strength, not causality claims
"""

import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import os
import json
from collections import deque
from typing import Dict, List, Optional


class AttentionWeightedExplainer:
    """
    Explains HGT link predictions using attention-weighted paths.

    Combines structural paths with embedding similarity and
    prediction scores for comprehensive explanations.
    """

    def __init__(self, model, predictor, data, device):
        """
        Initialize explainer with trained model and data.

        Args:
            model: Trained HGT model
            predictor: Link predictor module
            data: HeteroData graph
            device: torch device
        """
        self.model = model
        self.predictor = predictor
        self.data = data
        self.device = device

        self.model.eval()
        self.predictor.eval()

        # Precompute embeddings
        with torch.no_grad():
            self.x_dict = self.model(
                self.data.x_dict.copy(),
                self.data.edge_index_dict
            )

        # Build adjacency for path finding
        self._build_adjacency()

    def _build_adjacency(self):
        """Build adjacency list for BFS path finding."""
        self.adj = {}
        self.edge_type_map = {}

        for edge_type in self.data.edge_types:
            st, rel, dt = edge_type
            edge_index = self.data[edge_type].edge_index

            for i in range(edge_index.size(1)):
                src = edge_index[0, i].item()
                dst = edge_index[1, i].item()

                u = f"{st}_{src}"
                v = f"{dt}_{dst}"

                if u not in self.adj:
                    self.adj[u] = set()
                if v not in self.adj:
                    self.adj[v] = set()

                self.adj[u].add(v)
                self.adj[v].add(u)

                self.edge_type_map[(u, v)] = rel
                self.edge_type_map[(v, u)] = rel

    def find_paths(self, start_node: str, end_node: str,
                   max_depth: int = 3) -> List[List[str]]:
        """
        Find all paths between start and end nodes using BFS.

        Args:
            start_node: Starting node key (e.g., "TF_5")
            end_node: Target node key (e.g., "Enzyme_10")
            max_depth: Maximum path length

        Returns:
            List of paths, where each path is a list of node keys
        """
        queue = deque([(start_node, [start_node])])
        paths = []

        while queue:
            curr, path = queue.popleft()

            if len(path) > max_depth + 1:
                continue

            if curr == end_node and len(path) > 1:
                paths.append(path)
                continue

            if len(path) == max_depth + 1:
                continue

            for neighbor in self.adj.get(curr, []):
                if neighbor not in path:
                    queue.append((neighbor, path + [neighbor]))

        return paths

    def compute_path_attention(self, path: List[str]) -> Dict:
        """
        Compute attention-weighted importance for a path.

        Uses embedding similarity along path edges as attention proxy.

        Args:
            path: List of node keys forming a path

        Returns:
            Dictionary with path attention scores
        """
        if len(path) < 2:
            return {'total_attention': 0.0, 'edge_scores': []}

        edge_scores = []

        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]

            u_type, u_idx = u.split('_')
            v_type, v_idx = v.split('_')
            u_idx, v_idx = int(u_idx), int(v_idx)

            # Compute embedding similarity as attention proxy
            if u_type in self.x_dict and v_type in self.x_dict:
                u_emb = self.x_dict[u_type][u_idx]
                v_emb = self.x_dict[v_type][v_idx]

                u_norm = u_emb / (u_emb.norm() + 1e-8)
                v_norm = v_emb / (v_emb.norm() + 1e-8)

                similarity = torch.dot(u_norm, v_norm).item()
            else:
                similarity = 0.0

            rel = self.edge_type_map.get((u, v), 'unknown')

            edge_scores.append({
                'source': u,
                'target': v,
                'relation': rel,
                'attention': similarity
            })

        total_attention = np.mean([e['attention'] for e in edge_scores])

        return {
            'path': path,
            'total_attention': total_attention,
            'edge_scores': edge_scores,
            'path_length': len(path) - 1
        }

    def explain_prediction(self, src_type: str, src_idx: int,
                           dst_type: str, dst_idx: int,
                           max_paths: int = 10,
                           max_depth: int = 3) -> Dict:
        """
        Generate explanation for a specific TF-Enzyme prediction.

        Args:
            src_type: Source node type (e.g., 'TF')
            src_idx: Source node index
            dst_type: Target node type (e.g., 'Enzyme')
            dst_idx: Target node index
            max_paths: Maximum paths to return
            max_depth: Maximum path length

        Returns:
            Explanation dictionary with prediction score and paths
        """
        start = f"{src_type}_{src_idx}"
        end = f"{dst_type}_{dst_idx}"

        # Compute prediction score
        if src_type in self.x_dict and dst_type in self.x_dict:
            edge = torch.tensor([[src_idx], [dst_idx]], device=self.device)
            with torch.no_grad():
                score = self.predictor(
                    self.x_dict[src_type],
                    self.x_dict[dst_type],
                    edge
                ).sigmoid().item()
        else:
            score = 0.0

        # Find paths
        paths = self.find_paths(start, end, max_depth)

        # Score paths by attention
        scored_paths = []
        for path in paths:
            path_info = self.compute_path_attention(path)
            scored_paths.append(path_info)

        # Sort by attention
        scored_paths.sort(key=lambda x: -x['total_attention'])
        top_paths = scored_paths[:max_paths]

        return {
            'source': {'type': src_type, 'idx': src_idx},
            'target': {'type': dst_type, 'idx': dst_idx},
            'prediction_score': score,
            'num_paths_found': len(paths),
            'top_paths': top_paths,
            'interpretation_note': 'Attention scores indicate association strength, not causality'
        }


def explain_top_predictions(data, model, predictor, device,
                            src_type: str = 'TF',
                            dst_type: str = 'Enzyme',
                            top_k: int = 10) -> List[Dict]:
    """
    Generate explanations for top-k predicted TF-Enzyme pairs.

    Args:
        data: HeteroData graph
        model: Trained HGT model
        predictor: Link predictor
        device: torch device
        src_type: Source node type
        dst_type: Target node type
        top_k: Number of top predictions to explain

    Returns:
        List of explanation dictionaries
    """
    model.eval()
    predictor.eval()

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)

    if src_type not in x_dict or dst_type not in x_dict:
        return []

    src_emb = x_dict[src_type]
    dst_emb = x_dict[dst_type]

    n_src = min(100, src_emb.size(0))
    n_dst = min(100, dst_emb.size(0))

    # Score all pairs (sampled)
    predictions = []
    for i in range(n_src):
        for j in range(n_dst):
            edge = torch.tensor([[i], [j]], device=device)
            score = predictor(src_emb, dst_emb, edge).sigmoid().item()
            predictions.append({
                'src_idx': i,
                'dst_idx': j,
                'score': score
            })

    predictions.sort(key=lambda x: -x['score'])
    top_predictions = predictions[:top_k]

    # Generate explanations
    explainer = AttentionWeightedExplainer(model, predictor, data, device)

    explanations = []
    for pred in top_predictions:
        exp = explainer.explain_prediction(
            src_type, pred['src_idx'],
            dst_type, pred['dst_idx'],
            max_paths=5, max_depth=3
        )
        explanations.append(exp)

    return explanations


def visualize_attention_paths(explanation: Dict, output_path: str,
                              node_names: Optional[Dict] = None):
    """
    Visualize attention-weighted paths for an explanation.

    Args:
        explanation: Output from explain_prediction
        output_path: Path to save figure
        node_names: Optional mapping of node keys to names
    """
    G = nx.DiGraph()

    # Add nodes and edges from top paths
    for path_info in explanation['top_paths'][:5]:
        path = path_info['path']
        for node in path:
            if node not in G:
                node_type = node.split('_')[0]
                label = node_names.get(node, node) if node_names else node
                G.add_node(node, label=label, node_type=node_type)

        for edge in path_info['edge_scores']:
            G.add_edge(
                edge['source'],
                edge['target'],
                weight=edge['attention'],
                relation=edge['relation']
            )

    if G.number_of_nodes() == 0:
        print("No paths to visualize")
        return

    # Draw graph
    _, ax = plt.subplots(figsize=(14, 10))

    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Color by node type
    color_map = {'TF': '#a2d2ff', 'Enzyme': '#b5e48c',
                 'Signaling': '#fec89a', 'Protein': '#e5e5e5',
                 'Metabolite': '#ffccd5'}

    node_list = list(G.nodes())
    node_colors = [color_map.get(G.nodes[n].get('node_type', 'Protein'), '#e5e5e5')
                   for n in node_list]

    # Draw
    nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_color=node_colors,
                           node_size=2000, alpha=0.9, ax=ax)

    labels = {n: G.nodes[n].get('label', n) for n in node_list}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)

    # Edge weights
    edge_list = list(G.edges(data=True))
    edge_weights = [e[2].get('weight', 0.5) for e in edge_list]
    cmap = plt.get_cmap('Reds')
    edge_colors_arr = [cmap(w) for w in edge_weights]

    nx.draw_networkx_edges(G, pos, edgelist=[(e[0], e[1]) for e in edge_list],
                           edge_color=edge_colors_arr, width=2,
                           arrows=True, arrowsize=15, ax=ax,
                           connectionstyle="arc3,rad=0.1")

    # Title
    src = explanation['source']
    dst = explanation['target']
    score = explanation['prediction_score']
    ax.set_title(
        f"Attention Paths: {src['type']}[{src['idx']}] → {dst['type']}[{dst['idx']}]\n"
        f"Prediction Score: {score:.4f}",
        fontsize=12, fontweight='bold'
    )

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved path visualization to {output_path}")


def generate_explanation_report(data, model, predictor, device,
                                 output_dir: str,
                                 src_type: str = 'TF',
                                 dst_type: str = 'Enzyme',
                                 top_k: int = 10):
    """
    Generate comprehensive explanation report for top predictions.

    Args:
        data: HeteroData graph
        model: Trained HGT model
        predictor: Link predictor
        device: torch device
        output_dir: Output directory
        src_type: Source node type
        dst_type: Target node type
        top_k: Number of predictions to explain
    """
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print(f"EXPLANATION REPORT: {src_type} → {dst_type} Predictions")
    print("=" * 60)

    print("\nGenerating explanations for top predictions...")
    explanations = explain_top_predictions(
        data, model, predictor, device,
        src_type, dst_type, top_k
    )

    if not explanations:
        print("No explanations generated")
        return

    # Save explanations to JSON
    with open(os.path.join(output_dir, 'explanations.json'), 'w') as f:
        json.dump(explanations, f, indent=2, default=str)

    # Create summary table
    summary_rows = []
    for i, exp in enumerate(explanations, 1):
        src = exp['source']
        dst = exp['target']
        summary_rows.append({
            'rank': i,
            'source_type': src['type'],
            'source_idx': src['idx'],
            'target_type': dst['type'],
            'target_idx': dst['idx'],
            'prediction_score': exp['prediction_score'],
            'num_paths': exp['num_paths_found'],
            'top_path_attention': exp['top_paths'][0]['total_attention'] if exp['top_paths'] else 0
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(output_dir, 'explanation_summary.csv'), index=False)

    print(f"\nTop {top_k} {src_type}→{dst_type} Predictions:")
    print(summary_df.to_string(index=False))

    # Visualize top 3 predictions
    print("\nGenerating path visualizations...")
    for i, exp in enumerate(explanations[:3], 1):
        visualize_attention_paths(
            exp,
            os.path.join(output_dir, f'attention_paths_rank{i}.png')
        )

    # Create analysis notes
    notes = {
        'analysis_type': 'Attention-Weighted Path Explanation',
        'data_characteristic': 'Static comparison (72h single time point)',
        'interpretation_guideline': [
            'Attention scores indicate ASSOCIATION strength, not causality',
            'High prediction scores suggest strong multi-omics evidence',
            'Path analysis shows intermediate nodes connecting source and target',
            'For causal claims, perturbation experiments are required'
        ],
        'recommendations': [
            'Use docking/MD simulations for structural validation',
            'Search for time-series data to support temporal ordering',
            'Consider ChIP-seq/DAP-seq data for TF-target directness'
        ]
    }

    with open(os.path.join(output_dir, 'analysis_notes.json'), 'w') as f:
        json.dump(notes, f, indent=2)

    print(f"\nResults saved to {output_dir}/")
    return explanations
