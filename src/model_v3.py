"""
Model v3 for Ethylene-Isoflavonoid GNN Project

학술적으로 방어 가능한 모델 설계 (Academically Defensible Design)

Key Changes from v2:
1. Learnable embeddings as nn.Embedding (model parameters)
2. Omics features concat'd with learnable embeddings
3. Layer-wise edge type filtering (rxn_neighbor/TF only in Layer 1)
4. rxn_neighbor dropout for information leakage control
5. TF domain one-hot features support

Usage:
    from src.model_v3 import HGTv3, LinkPredictor
    model = HGTv3(data.metadata(), num_nodes_dict, config)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import HGTConv, Linear
from typing import Dict, List, Optional, Tuple

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import V3_CONFIG, TRAINING_CONFIG


class HGTv3(nn.Module):
    """
    Heterogeneous Graph Transformer v3 with academically defensible design.

    Key Features:
    1. Learnable node embeddings (nn.Embedding) + omics features
    2. Layer-wise edge type filtering
    3. rxn_neighbor dropout for information leakage control

    Args:
        metadata: Graph metadata (node_types, edge_types)
        num_nodes_dict: Dict mapping node_type → num_nodes
        hidden_channels: Hidden dimension (default: 64)
        out_channels: Output dimension (default: 64)
        num_heads: Number of attention heads (default: 4)
        num_layers: Number of HGT layers (default: 2)
        rxn_neighbor_dropout: Dropout rate for rxn_neighbor edges
        rxn_neighbor_dropout_mode: 'per_run_fixed' or 'per_epoch'
    """

    def __init__(
        self,
        metadata: Tuple[List[str], List[Tuple[str, str, str]]],
        num_nodes_dict: Dict[str, int],
        hidden_channels: int = None,
        out_channels: int = None,
        num_heads: int = None,
        num_layers: int = 2,
        rxn_neighbor_dropout: float = None,
        rxn_neighbor_dropout_mode: str = None,
        met_omics_dim: int = None,
        enz_omics_dim: int = None,
        tf_domain_dim: int = None,
    ):
        super().__init__()

        # Use config defaults
        hidden_channels = hidden_channels or TRAINING_CONFIG['hidden_channels']
        out_channels = out_channels or TRAINING_CONFIG['out_channels']
        num_heads = num_heads or TRAINING_CONFIG['num_heads']
        rxn_neighbor_dropout = rxn_neighbor_dropout if rxn_neighbor_dropout is not None \
            else V3_CONFIG['rxn_neighbor_dropout']
        rxn_neighbor_dropout_mode = rxn_neighbor_dropout_mode or \
            V3_CONFIG['rxn_neighbor_dropout_mode']

        # Feature dimensions
        met_omics_dim = met_omics_dim or V3_CONFIG['met_omics_dim']
        enz_omics_dim = enz_omics_dim or V3_CONFIG['enz_omics_dim']
        tf_domain_dim = tf_domain_dim or V3_CONFIG['tf_domain_dim']

        self.hidden_channels = hidden_channels
        self.out_channels = out_channels
        self.num_layers = num_layers
        self.rxn_neighbor_dropout = rxn_neighbor_dropout
        self.rxn_neighbor_dropout_mode = rxn_neighbor_dropout_mode

        # Store metadata
        self.node_types = metadata[0]
        self.edge_types = metadata[1]

        # Layer-wise edge types from config
        self.layer1_edge_types = set(V3_CONFIG['layer1_edge_types'])
        self.layer2_edge_types = set(V3_CONFIG['layer2_edge_types'])

        # Learnable embeddings for each node type
        # These are model parameters (updated during training)
        self.embeddings = nn.ModuleDict()
        self.omics_dims = {}

        for node_type in self.node_types:
            num_nodes = num_nodes_dict.get(node_type, 0)
            if num_nodes == 0:
                continue

            # Determine omics dimension for this node type
            if node_type == 'Metabolite':
                omics_dim = met_omics_dim
            elif node_type == 'Enzyme':
                omics_dim = enz_omics_dim
            elif node_type == 'TF':
                omics_dim = tf_domain_dim
            else:
                omics_dim = 0

            self.omics_dims[node_type] = omics_dim
            learnable_dim = hidden_channels - omics_dim

            # Create learnable embedding
            self.embeddings[node_type] = nn.Embedding(num_nodes, learnable_dim)

        # Input projection layer for each node type
        self.lin_dict = nn.ModuleDict()
        for node_type in self.node_types:
            if node_type in self.embeddings:
                self.lin_dict[node_type] = Linear(hidden_channels, hidden_channels)

        # HGT convolution layers
        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HGTConv(hidden_channels, hidden_channels, metadata, heads=num_heads)
            self.convs.append(conv)

        # Output projection
        self.out_lin = Linear(hidden_channels, out_channels)

        # Per-run fixed mask for rxn_neighbor dropout (set during first forward)
        self.rxn_neighbor_mask: Optional[torch.Tensor] = None
        self._rxn_neighbor_num_edges: Optional[int] = None

    def _get_device(self, data) -> torch.device:
        """Get device from data object."""
        for node_type in data.node_types:
            if hasattr(data[node_type], 'x') and data[node_type].x is not None:
                return data[node_type].x.device
            if hasattr(data[node_type], 'node_ids'):
                return data[node_type].node_ids.device
            if hasattr(data[node_type], 'omics_x'):
                return data[node_type].omics_x.device
        return torch.device('cpu')

    def _get_node_features(
        self,
        data,
        node_type: str,
    ) -> torch.Tensor:
        """
        Get node features by concatenating learnable embeddings with omics features.

        Linear concat (no MLP projection) - omics features are few, MLP would over-parameterize.

        Args:
            data: HeteroData object
            node_type: Node type string

        Returns:
            Node feature tensor [num_nodes, hidden_channels]
        """
        num_nodes = data[node_type].num_nodes

        # Check if we have a learnable embedding for this type
        if node_type not in self.embeddings:
            # Fallback: use existing x if available, else random
            if hasattr(data[node_type], 'x') and data[node_type].x is not None:
                x = data[node_type].x
                if x.size(1) != self.hidden_channels:
                    # Project to hidden_channels
                    device = x.device
                    proj = nn.Linear(x.size(1), self.hidden_channels).to(device)
                    return proj(x)
                return x
            else:
                # Get device from any available tensor
                device = self._get_device(data)
                return torch.randn(num_nodes, self.hidden_channels, device=device)

        # Get node_ids for embedding lookup
        if hasattr(data[node_type], 'node_ids'):
            node_ids = data[node_type].node_ids
        else:
            device = self._get_device(data)
            node_ids = torch.arange(num_nodes, device=device)

        learnable = self.embeddings[node_type](node_ids)

        # Get omics features
        omics_dim = self.omics_dims.get(node_type, 0)
        if omics_dim > 0:
            if node_type == 'Metabolite' and hasattr(data[node_type], 'omics_x'):
                omics = data[node_type].omics_x
            elif node_type == 'Enzyme' and hasattr(data[node_type], 'omics_x'):
                omics = data[node_type].omics_x
            elif node_type == 'TF' and hasattr(data[node_type], 'domain_x'):
                omics = data[node_type].domain_x
            else:
                # No omics available, use zeros
                omics = torch.zeros(learnable.size(0), omics_dim, device=learnable.device)

            # Concat: [learnable, omics]
            x = torch.cat([learnable, omics], dim=-1)
        else:
            x = learnable

        return x

    def _filter_edges_for_layer(
        self,
        edge_index_dict: Dict[Tuple[str, str, str], torch.Tensor],
        layer_idx: int,
    ) -> Dict[Tuple[str, str, str], torch.Tensor]:
        """
        Filter edge types based on layer index.

        Layer 1: All edge types (local context)
        Layer 2+: Strong relationships only (PPI, catalyzes_R)

        Rationale: rxn_neighbor/TF are weak associations; multi-hop
        propagation causes noise and representation pollution.

        Args:
            edge_index_dict: Original edge indices
            layer_idx: Current layer index (0-based)

        Returns:
            Filtered edge index dict
        """
        if layer_idx == 0:
            allowed_types = self.layer1_edge_types
        else:
            allowed_types = self.layer2_edge_types

        filtered = {}
        for edge_type, edge_index in edge_index_dict.items():
            if edge_type in allowed_types:
                filtered[edge_type] = edge_index

        return filtered

    def _apply_rxn_neighbor_dropout(
        self,
        edge_index_dict: Dict[Tuple[str, str, str], torch.Tensor],
    ) -> Dict[Tuple[str, str, str], torch.Tensor]:
        """
        Apply dropout to rxn_neighbor edges during training.

        This controls information leakage from metabolite community structure.

        Args:
            edge_index_dict: Original edge indices

        Returns:
            Edge index dict with rxn_neighbor edges masked
        """
        if not self.training or self.rxn_neighbor_dropout <= 0:
            return edge_index_dict

        rxn_key = ('Metabolite', 'rxn_neighbor', 'Metabolite')
        if rxn_key not in edge_index_dict:
            return edge_index_dict

        edge_index = edge_index_dict[rxn_key]
        num_edges = edge_index.size(1)

        if self.rxn_neighbor_dropout_mode == 'per_run_fixed':
            # Fixed mask for reproducibility (set once per run)
            if self.rxn_neighbor_mask is None or self._rxn_neighbor_num_edges != num_edges:
                mask = torch.rand(num_edges, device=edge_index.device) > self.rxn_neighbor_dropout
                self.rxn_neighbor_mask = mask
                self._rxn_neighbor_num_edges = num_edges

            mask = self.rxn_neighbor_mask
        else:
            # Per-epoch dropout (regularization effect)
            mask = torch.rand(num_edges, device=edge_index.device) > self.rxn_neighbor_dropout

        # Apply mask
        result = dict(edge_index_dict)
        result[rxn_key] = edge_index[:, mask]
        return result

    def forward(
        self,
        data,
        return_all_layers: bool = False,
    ) -> Dict[str, torch.Tensor]:
        """
        Forward pass through HGT layers.

        Args:
            data: HeteroData object with node features and edge indices
            return_all_layers: If True, return embeddings from all layers

        Returns:
            Dict mapping node_type → embeddings
        """
        # Build initial node features (learnable + omics)
        x_dict = {}
        for node_type in self.node_types:
            if node_type in self.embeddings:
                x_dict[node_type] = self._get_node_features(data, node_type)

        # Apply input projection with ReLU
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
            if node_type in self.lin_dict
        }

        # Get edge indices
        edge_index_dict = data.edge_index_dict

        # Apply rxn_neighbor dropout
        edge_index_dict = self._apply_rxn_neighbor_dropout(edge_index_dict)

        # Store layer outputs if needed
        layer_outputs = [x_dict] if return_all_layers else None

        # Message passing through HGT layers
        for layer_idx, conv in enumerate(self.convs):
            # Filter edges for this layer
            filtered_edges = self._filter_edges_for_layer(edge_index_dict, layer_idx)

            # Apply HGT convolution
            x_dict = conv(x_dict, filtered_edges)

            if return_all_layers:
                layer_outputs.append(x_dict)

        if return_all_layers:
            return layer_outputs

        return x_dict

    def get_embeddings(self, data) -> Dict[str, torch.Tensor]:
        """Get final node embeddings."""
        return self.forward(data)


class LinkPredictor(nn.Module):
    """
    Link predictor using dot product scoring.

    For Enzyme→Metabolite link prediction:
        score(e, m) = <h_e, h_m>
    """

    def __init__(self, in_channels: int = None):
        super().__init__()
        self.in_channels = in_channels or TRAINING_CONFIG['out_channels']

    def forward(
        self,
        x_src: torch.Tensor,
        x_dst: torch.Tensor,
        edge_label_index: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute link prediction scores.

        Args:
            x_src: Source node embeddings [num_src, dim]
            x_dst: Target node embeddings [num_dst, dim]
            edge_label_index: Edge indices [2, num_edges]

        Returns:
            Scores [num_edges]
        """
        row, col = edge_label_index
        src_feats = x_src[row]
        dst_feats = x_dst[col]
        return (src_feats * dst_feats).sum(dim=-1)

    def predict_all(
        self,
        x_src: torch.Tensor,
        x_dst: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute all pairwise scores.

        Args:
            x_src: Source embeddings [num_src, dim]
            x_dst: Target embeddings [num_dst, dim]

        Returns:
            Scores [num_src, num_dst]
        """
        return torch.mm(x_src, x_dst.t())


class HGTv3Ablation(HGTv3):
    """
    HGT v3 with ablation support.

    Allows disabling specific edge types for ablation studies.
    """

    def __init__(
        self,
        *args,
        use_rxn_neighbor: bool = True,
        use_tf_edges: bool = True,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.use_rxn_neighbor = use_rxn_neighbor
        self.use_tf_edges = use_tf_edges

    def _filter_edges_for_layer(
        self,
        edge_index_dict: Dict[Tuple[str, str, str], torch.Tensor],
        layer_idx: int,
    ) -> Dict[Tuple[str, str, str], torch.Tensor]:
        """Filter edges with ablation support."""
        filtered = super()._filter_edges_for_layer(edge_index_dict, layer_idx)

        # Ablation: Remove rxn_neighbor
        if not self.use_rxn_neighbor:
            rxn_key = ('Metabolite', 'rxn_neighbor', 'Metabolite')
            if rxn_key in filtered:
                del filtered[rxn_key]

        # Ablation: Remove TF edges
        if not self.use_tf_edges:
            tf_keys = [
                ('TF', 'associates', 'Enzyme'),
                ('Enzyme', 'rev_associates', 'TF'),
            ]
            for key in tf_keys:
                if key in filtered:
                    del filtered[key]

        return filtered


def create_model_v3(
    data,
    num_layers: int = 2,
    use_rxn_neighbor: bool = True,
    use_tf_edges: bool = True,
) -> Tuple[nn.Module, nn.Module]:
    """
    Factory function to create HGT v3 model and link predictor.

    Args:
        data: HeteroData object
        num_layers: Number of GNN layers
        use_rxn_neighbor: Whether to use rxn_neighbor edges
        use_tf_edges: Whether to use TF edges

    Returns:
        Tuple of (model, predictor)
    """
    # Get num_nodes for each type
    num_nodes_dict = {}
    for node_type in data.node_types:
        num_nodes_dict[node_type] = data[node_type].num_nodes

    # Create model
    if use_rxn_neighbor and use_tf_edges:
        model = HGTv3(
            metadata=data.metadata(),
            num_nodes_dict=num_nodes_dict,
            num_layers=num_layers,
        )
    else:
        model = HGTv3Ablation(
            metadata=data.metadata(),
            num_nodes_dict=num_nodes_dict,
            num_layers=num_layers,
            use_rxn_neighbor=use_rxn_neighbor,
            use_tf_edges=use_tf_edges,
        )

    predictor = LinkPredictor()

    return model, predictor
