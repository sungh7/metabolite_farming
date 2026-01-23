"""
Unified Data Splitting for Link Prediction

This module provides consistent data splitting strategies
for train/validation/test splits across all model scripts.
"""

import torch
import numpy as np
from typing import Literal, Optional, Tuple
from sklearn.model_selection import train_test_split
from torch_geometric.data import HeteroData

from src.config import SPLIT_CONFIG


class DataSplitter:
    """
    Unified data splitter for heterogeneous graphs.

    Supports multiple split strategies:
    - 'node_split': Split nodes, edges follow (OGB/PyG standard)
    - 'edge_split': Randomly split edges directly
    - 'random_link_split': PyG's RandomLinkSplit transform

    Example:
        >>> splitter = DataSplitter(strategy='node_split')
        >>> train_data, val_data, test_data = splitter.split(data)
    """

    def __init__(
        self,
        strategy: Literal['node_split', 'edge_split', 'random_link_split'] = None,
        train_ratio: float = None,
        val_ratio: float = None,
        test_ratio: float = None,
        disjoint: bool = None,
        seed: int = 42,
    ):
        """
        Initialize the data splitter.

        Args:
            strategy: Split strategy
            train_ratio: Fraction for training (default: 0.70)
            val_ratio: Fraction for validation (default: 0.15)
            test_ratio: Fraction for testing (default: 0.15)
            disjoint: Whether to use disjoint node splits
            seed: Random seed for reproducibility
        """
        config = SPLIT_CONFIG

        self.strategy = strategy or config['strategy']
        self.train_ratio = train_ratio or config['train_ratio']
        self.val_ratio = val_ratio or config['val_ratio']
        self.test_ratio = test_ratio or config['test_ratio']
        self.disjoint = disjoint if disjoint is not None else config['disjoint']
        self.seed = seed

        # Validate ratios
        total = self.train_ratio + self.val_ratio + self.test_ratio
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f"Split ratios must sum to 1.0, got {total}")

    def split(
        self,
        data: HeteroData,
        edge_type: Tuple[str, str, str] = ('Enzyme', 'catalyzes', 'Metabolite'),
        source_node_type: str = 'Enzyme',
    ) -> Tuple[HeteroData, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Split the data according to the configured strategy.

        Args:
            data: HeteroData graph
            edge_type: Edge type to split on
            source_node_type: Node type to use for node-based splitting

        Returns:
            Tuple of (train_data, train_edges, val_edges, test_edges)
        """
        if self.strategy == 'node_split':
            return self._node_split(data, edge_type, source_node_type)
        elif self.strategy == 'edge_split':
            return self._edge_split(data, edge_type)
        elif self.strategy == 'random_link_split':
            return self._random_link_split(data, edge_type)
        else:
            raise ValueError(f"Unknown strategy: {self.strategy}")

    def _node_split(
        self,
        data: HeteroData,
        edge_type: Tuple[str, str, str],
        source_node_type: str,
    ) -> Tuple[HeteroData, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Node-disjoint split: divide nodes, edges follow.

        This is the recommended strategy for inductive evaluation,
        where we want to assess performance on unseen nodes.
        """
        device = data[source_node_type].x.device
        num_nodes = data[source_node_type].num_nodes

        # Create node indices split
        train_idx, val_idx, test_idx = self._create_node_splits(num_nodes)

        # Create masks
        train_mask = torch.zeros(num_nodes, dtype=torch.bool, device=device)
        val_mask = torch.zeros(num_nodes, dtype=torch.bool, device=device)
        test_mask = torch.zeros(num_nodes, dtype=torch.bool, device=device)

        train_mask[train_idx] = True
        val_mask[val_idx] = True
        test_mask[test_idx] = True

        # Split edges based on source node membership
        edge_index = data[edge_type].edge_index
        src = edge_index[0]

        train_edge_mask = train_mask[src]
        val_edge_mask = val_mask[src]
        test_edge_mask = test_mask[src]

        train_edges = edge_index[:, train_edge_mask]
        val_edges = edge_index[:, val_edge_mask]
        test_edges = edge_index[:, test_edge_mask]

        # Create training graph (only train edges)
        train_data = data.clone()
        train_data[edge_type].edge_index = train_edges

        # Handle edge weights if present
        if hasattr(data[edge_type], 'edge_weight'):
            train_data[edge_type].edge_weight = data[edge_type].edge_weight[train_edge_mask]

        # Handle reverse edges
        src_type, rel, dst_type = edge_type
        rev_edge_type = (dst_type, f'rev_{rel}', src_type)

        if rev_edge_type in data.edge_types:
            rev_index = data[rev_edge_type].edge_index
            rev_mask = train_mask[rev_index[1]]  # Reverse: dst is now at index 1
            train_data[rev_edge_type].edge_index = rev_index[:, rev_mask]

            if hasattr(data[rev_edge_type], 'edge_weight'):
                train_data[rev_edge_type].edge_weight = data[rev_edge_type].edge_weight[rev_mask]

        # Store node masks for evaluation
        train_data._train_mask = train_mask
        train_data._val_mask = val_mask
        train_data._test_mask = test_mask

        return train_data, train_edges, val_edges, test_edges

    def _edge_split(
        self,
        data: HeteroData,
        edge_type: Tuple[str, str, str],
    ) -> Tuple[HeteroData, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Edge-based split: randomly divide edges.
        """
        edge_index = data[edge_type].edge_index
        num_edges = edge_index.size(1)

        # Permute edge indices
        np.random.seed(self.seed)
        perm = np.random.permutation(num_edges)

        # Calculate split points
        train_end = int(self.train_ratio * num_edges)
        val_end = train_end + int(self.val_ratio * num_edges)

        train_perm = perm[:train_end]
        val_perm = perm[train_end:val_end]
        test_perm = perm[val_end:]

        train_edges = edge_index[:, train_perm]
        val_edges = edge_index[:, val_perm]
        test_edges = edge_index[:, test_perm]

        # Create training graph
        train_data = data.clone()
        train_data[edge_type].edge_index = train_edges

        return train_data, train_edges, val_edges, test_edges

    def _random_link_split(
        self,
        data: HeteroData,
        edge_type: Tuple[str, str, str],
    ) -> Tuple[HeteroData, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        PyG's RandomLinkSplit transform.
        """
        import torch_geometric.transforms as T

        transform = T.RandomLinkSplit(
            num_val=self.val_ratio,
            num_test=self.test_ratio,
            is_undirected=False,
            edge_types=[edge_type],
            add_negative_train_samples=False,
        )

        train_data, val_data, test_data = transform(data)

        # Extract edges
        train_edges = train_data[edge_type].edge_label_index
        val_edges = val_data[edge_type].edge_label_index
        test_edges = test_data[edge_type].edge_label_index

        return train_data, train_edges, val_edges, test_edges

    def _create_node_splits(
        self,
        num_nodes: int,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Create train/val/test node index splits.
        """
        indices = np.arange(num_nodes)

        # First split: train vs (val+test)
        train_idx, temp_idx = train_test_split(
            indices,
            train_size=self.train_ratio,
            random_state=self.seed
        )

        # Second split: val vs test
        val_ratio_adjusted = self.val_ratio / (self.val_ratio + self.test_ratio)
        val_idx, test_idx = train_test_split(
            temp_idx,
            train_size=val_ratio_adjusted,
            random_state=self.seed
        )

        return train_idx, val_idx, test_idx

    def get_node_masks(
        self,
        train_data: HeteroData,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Get node masks from split data (if available).
        """
        if hasattr(train_data, '_train_mask'):
            return (
                train_data._train_mask,
                train_data._val_mask,
                train_data._test_mask
            )
        raise AttributeError("Node masks not available. Use node_split strategy.")


def create_splits(
    num_nodes: int,
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convenience function for creating node splits.

    This is a drop-in replacement for the function in eval_standard.py.

    Args:
        num_nodes: Total number of nodes
        train_ratio: Fraction for training
        val_ratio: Fraction for validation
        test_ratio: Fraction for testing
        seed: Random seed

    Returns:
        Tuple of (train_indices, val_indices, test_indices)
    """
    splitter = DataSplitter(
        strategy='node_split',
        train_ratio=train_ratio,
        val_ratio=val_ratio,
        test_ratio=test_ratio,
        seed=seed,
    )
    return splitter._create_node_splits(num_nodes)
