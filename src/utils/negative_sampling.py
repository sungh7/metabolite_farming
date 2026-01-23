"""
Unified Negative Sampling for Link Prediction

This module provides consistent negative sampling strategies
for training and evaluation across all model scripts.
"""

import torch
import random
from typing import Literal, Optional, Dict, Set
from src.config import NEGATIVE_SAMPLING_CONFIG


class NegativeSampler:
    """
    Unified negative sampler for link prediction tasks.

    Supports four strategies:
    - 'random': Uniformly random negative samples
    - 'hard': Hard negatives by sampling nearby metabolites (±offset)
    - 'mixed': Combination of random and hard negatives
    - 'ec_class': EC-class aware sampling (avoids same EC class metabolites)

    Example:
        >>> sampler = NegativeSampler(strategy='hard')
        >>> neg_src, neg_dst = sampler.sample(pos_edges, num_metabolites, device)

        # For EC-class strategy, provide mappings:
        >>> sampler = NegativeSampler(
        ...     strategy='ec_class',
        ...     ec_to_indices=ec_to_indices,
        ...     met_to_ecs=met_to_ecs
        ... )
    """

    def __init__(
        self,
        strategy: Literal['random', 'hard', 'mixed', 'ec_class'] = None,
        hard_offset_min: int = None,
        hard_offset_max: int = None,
        hard_ratio: float = None,
        neg_ratio: float = None,
        ec_to_indices: Optional[Dict[str, list]] = None,
        met_to_ecs: Optional[Dict[int, Set[str]]] = None,
    ):
        """
        Initialize the negative sampler.

        Args:
            strategy: Sampling strategy ('random', 'hard', 'mixed', or 'ec_class')
            hard_offset_min: Minimum offset for hard negatives (default: 1)
            hard_offset_max: Maximum offset for hard negatives (default: 5)
            hard_ratio: Ratio of hard negatives when using 'mixed' strategy
            neg_ratio: Ratio of negative to positive samples
            ec_to_indices: EC class to metabolite indices mapping (for ec_class strategy)
            met_to_ecs: Metabolite index to EC classes mapping (for ec_class strategy)
        """
        config = NEGATIVE_SAMPLING_CONFIG

        self.strategy = strategy or config['strategy']
        self.hard_offset_min = hard_offset_min or config['hard_offset_min']
        self.hard_offset_max = hard_offset_max or config['hard_offset_max']
        self.hard_ratio = hard_ratio or config['hard_ratio']
        self.neg_ratio = neg_ratio or config['neg_ratio']

        # EC-class specific mappings
        self.ec_to_indices = ec_to_indices or {}
        self.met_to_ecs = met_to_ecs or {}
        self.ec_class_fallback = config.get('ec_class_fallback', 'random')

    def sample(
        self,
        pos_edge_index: torch.Tensor,
        num_target_nodes: int,
        device: torch.device,
        num_samples: Optional[int] = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Sample negative edges.

        Args:
            pos_edge_index: Positive edge indices [2, num_edges]
            num_target_nodes: Total number of target nodes (metabolites)
            device: Device to create tensors on
            num_samples: Number of negative samples (default: same as positive)

        Returns:
            Tuple of (negative source indices, negative target indices)
        """
        num_pos = pos_edge_index.size(1)
        num_neg = num_samples or int(num_pos * self.neg_ratio)

        if self.strategy == 'random':
            return self._sample_random(pos_edge_index, num_target_nodes, num_neg, device)
        elif self.strategy == 'hard':
            return self._sample_hard(pos_edge_index, num_target_nodes, num_neg, device)
        elif self.strategy == 'mixed':
            return self._sample_mixed(pos_edge_index, num_target_nodes, num_neg, device)
        elif self.strategy == 'ec_class':
            return self._sample_ec_class(pos_edge_index, num_target_nodes, num_neg, device)
        else:
            raise ValueError(f"Unknown strategy: {self.strategy}")

    def _sample_random(
        self,
        pos_edge_index: torch.Tensor,
        num_target_nodes: int,
        num_neg: int,
        device: torch.device,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Sample uniformly random negative edges.

        For each positive edge (src, dst), sample a random target node.
        """
        neg_src = pos_edge_index[0][:num_neg]
        neg_dst = torch.randint(0, num_target_nodes, (num_neg,), device=device)
        return neg_src, neg_dst

    def _sample_hard(
        self,
        pos_edge_index: torch.Tensor,
        num_target_nodes: int,
        num_neg: int,
        device: torch.device,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Sample hard negative edges.

        For each positive edge (src, dst), sample a nearby target node
        by adding a small random offset to the target index. This simulates
        "same pathway, different reaction" without explicit pathway data.
        """
        neg_src = pos_edge_index[0][:num_neg]
        pos_dst = pos_edge_index[1][:num_neg]

        # Random offset in range [hard_offset_min, hard_offset_max]
        offset_magnitude = torch.randint(
            self.hard_offset_min,
            self.hard_offset_max + 1,
            (num_neg,),
            device=device
        )

        # Random sign (±1)
        sign = 2 * torch.randint(0, 2, (num_neg,), device=device) - 1

        # Apply offset with wrap-around
        neg_dst = (pos_dst + offset_magnitude * sign) % num_target_nodes

        return neg_src, neg_dst

    def _sample_mixed(
        self,
        pos_edge_index: torch.Tensor,
        num_target_nodes: int,
        num_neg: int,
        device: torch.device,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Sample mixed negative edges (combination of hard and random).
        """
        num_hard = int(num_neg * self.hard_ratio)
        num_random = num_neg - num_hard

        # Sample hard negatives
        hard_src, hard_dst = self._sample_hard(
            pos_edge_index, num_target_nodes, num_hard, device
        )

        # Sample random negatives
        random_src, random_dst = self._sample_random(
            pos_edge_index[:, num_hard:], num_target_nodes, num_random, device
        )

        # Concatenate
        neg_src = torch.cat([hard_src, random_src])
        neg_dst = torch.cat([hard_dst, random_dst])

        return neg_src, neg_dst

    def _sample_ec_class(
        self,
        pos_edge_index: torch.Tensor,
        num_target_nodes: int,
        num_neg: int,
        device: torch.device,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Sample EC-class aware negative edges.

        For each positive edge (src, dst), sample a metabolite that does NOT
        belong to the same EC class as the positive target. This avoids
        false negatives where the sampled 'negative' is actually connected
        to similar enzymatic reactions.

        Falls back to random sampling if EC mapping is not available for a
        metabolite or if all metabolites share the same EC class.
        """
        neg_src = pos_edge_index[0][:num_neg]
        pos_dst = pos_edge_index[1][:num_neg]

        # If no EC mapping available, fall back to random or hard sampling
        if not self.met_to_ecs or not self.ec_to_indices:
            if self.ec_class_fallback == 'hard':
                return self._sample_hard(pos_edge_index, num_target_nodes, num_neg, device)
            else:
                return self._sample_random(pos_edge_index, num_target_nodes, num_neg, device)

        all_mets = set(range(num_target_nodes))
        neg_dst_list = []

        for dst_idx in pos_dst.tolist():
            # Get EC classes for current metabolite
            current_ecs = self.met_to_ecs.get(dst_idx, set())

            if not current_ecs:
                # No EC info for this metabolite, sample randomly
                neg_dst_list.append(random.randint(0, num_target_nodes - 1))
                continue

            # Collect all metabolites that share any EC class with current
            excluded = set()
            for ec in current_ecs:
                excluded.update(self.ec_to_indices.get(ec, []))

            # Also exclude the positive target itself
            excluded.add(dst_idx)

            # Get candidates (metabolites not in same EC class)
            candidates = list(all_mets - excluded)

            if candidates:
                neg_dst_list.append(random.choice(candidates))
            else:
                # All metabolites share EC class, fall back to random
                # (but at least avoid the exact same metabolite)
                fallback = random.randint(0, num_target_nodes - 1)
                while fallback == dst_idx and num_target_nodes > 1:
                    fallback = random.randint(0, num_target_nodes - 1)
                neg_dst_list.append(fallback)

        neg_dst = torch.tensor(neg_dst_list, device=device, dtype=torch.long)
        return neg_src, neg_dst


def sample_negatives(
    pos_edge_index: torch.Tensor,
    num_target_nodes: int,
    device: torch.device,
    strategy: str = 'hard',
    **kwargs
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Convenience function for one-off negative sampling.

    Args:
        pos_edge_index: Positive edge indices [2, num_edges]
        num_target_nodes: Total number of target nodes
        device: Device to create tensors on
        strategy: Sampling strategy
        **kwargs: Additional arguments passed to NegativeSampler

    Returns:
        Tuple of (negative source indices, negative target indices)
    """
    sampler = NegativeSampler(strategy=strategy, **kwargs)
    return sampler.sample(pos_edge_index, num_target_nodes, device)
